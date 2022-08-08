//! Taken from https://crates.io/crates/wavefront but making it compile on stable

use core::{
    fmt,
    num::NonZeroUsize,
    ops::{Deref, DerefMut},
};
use std::{
    collections::HashMap,
    error,
    fs::File,
    io::{self, Read, Write},
    path::Path,
};

/// A number used to index into vertex attribute arrays.
pub type Index = usize;

/// An error that may be encountered while attempting to parse an OBJ.
#[derive(Debug)]
#[non_exhaustive]
pub enum Error {
    /// An error occurred when attempting to perform I/O.
    Io(io::Error),
    /// Expected a term on the given line but no term was found instead.
    ExpectedTerm(usize),
    /// Expected an index on the given line but something else was found.
    ExpectedIdx(usize),
    /// Expected a name but something else was found instead.
    ExpectedName(usize),
    /// An invalid index was encountered.
    InvalidIndex(isize),
}

impl From<io::Error> for Error {
    fn from(e: io::Error) -> Self {
        Error::Io(e)
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Error::Io(e) => write!(f, "{}", e),
            Error::ExpectedTerm(line) => write!(f, "Expected term on line {}", line),
            Error::ExpectedIdx(line) => write!(f, "Expected index on line {}", line),
            Error::ExpectedName(line) => {
                write!(f, "Expected object or group name on line {}", line)
            }
            Error::InvalidIndex(idx) => write!(f, "Invalid index '{}'", idx),
        }
    }
}

impl error::Error for Error {}

/// A struct representing the contents of a parsed OBJ file.
#[derive(Clone)]
pub struct Obj {
    buffers: Buffers,
    objects: HashMap<String, HashMap<String, Vec<VertexRange>>>,
}

impl Obj {
    /// Read an OBJ from a file.
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, Error> {
        Self::from_reader(io::BufReader::new(File::open(path)?))
    }

    /// Read an OBJ from a reader (something implementing [`std::io::Read`]).
    pub fn from_reader<R: Read>(mut reader: R) -> Result<Self, Error> {
        let mut buf = String::new();
        reader.read_to_string(&mut buf)?;
        Self::from_lines(buf.lines())
    }

    /// Read an OBJ from an iterator over its lines.
    pub fn from_lines<I: Iterator<Item = L>, L: AsRef<str>>(lines: I) -> Result<Self, Error> {
        let mut positions = Vec::new();
        let mut uvs = Vec::new();
        let mut normals = Vec::new();
        let mut vertices = Vec::new();
        let mut objects = HashMap::new();

        let mut object = (None, HashMap::new());

        let mut default_group = Vec::new();
        let mut groups = HashMap::<_, Vec<VertexRange>>::new();
        let mut selected_groups = Vec::new();

        for (i, line) in lines.enumerate() {
            let line = line.as_ref();
            let line_num = i + 1;
            let mut terms = line.split_ascii_whitespace();
            match terms.next() {
                Some("v") => {
                    let mut nums = terms.map_while(|t| t.parse().ok());
                    positions.push([
                        nums.next().unwrap_or(0.0),
                        nums.next().unwrap_or(0.0),
                        nums.next().unwrap_or(0.0),
                    ]);
                }
                Some("vt") => {
                    let mut nums = terms.map_while(|t| t.parse().ok());
                    uvs.push([
                        nums.next().unwrap_or(0.0),
                        nums.next().unwrap_or(0.0),
                        nums.next().unwrap_or(0.0),
                    ]);
                }
                Some("vn") => {
                    let mut nums = terms.map_while(|t| t.parse().ok());
                    normals.push([
                        nums.next().unwrap_or(0.0),
                        nums.next().unwrap_or(0.0),
                        nums.next().unwrap_or(0.0),
                    ]);
                }
                Some("f") => {
                    let parse_vert = |lengths: [usize; 3], v: &str| {
                        v.split('/')
                            .enumerate()
                            .take(3)
                            .map(|(i, idx)| match idx.trim() {
                                "" => Ok(None),
                                s => s
                                    .parse::<isize>()
                                    .map_err(|_| Error::ExpectedIdx(line_num))
                                    .and_then(|idx| {
                                        Ok(Some(if idx >= 0 {
                                            NonZeroUsize::new(idx as usize)
                                                .ok_or_else(|| Error::InvalidIndex(idx))?
                                        } else {
                                            lengths[i]
                                                .checked_sub((-idx - 1) as usize)
                                                .map(|idx| NonZeroUsize::new(idx).unwrap())
                                                .ok_or_else(|| Error::InvalidIndex(idx))?
                                        }))
                                    }),
                            })
                            .collect::<Result<Vec<_>, Error>>()
                    };

                    let lengths = [positions.len(), uvs.len(), normals.len()];
                    let poly_start = vertices.len();

                    for term in terms {
                        let v = parse_vert(lengths, term)?;

                        vertices.push((
                            // Position
                            v.get(0)
                                .copied()
                                .flatten()
                                .ok_or_else(|| Error::ExpectedIdx(line_num))?,
                            // Uv
                            v.get(1).copied().flatten(),
                            // Normal
                            v.get(2).copied().flatten(),
                        ));
                    }

                    let poly = VertexRange {
                        start: poly_start,
                        end: vertices.len(),
                    };

                    if selected_groups.len() == 0 {
                        default_group.push(poly);
                    } else {
                        selected_groups
                            .iter()
                            .for_each(|g| groups.get_mut(g).unwrap().push(poly));
                    }
                }
                Some("g") => {
                    selected_groups = terms
                        .map_while(|t| Some(t).filter(|t| util::name_is_valid(t)))
                        .map(|g| {
                            groups.entry(g.to_string()).or_default();
                            g.to_string()
                        })
                        .collect();
                }
                Some("o") => {
                    // Clean up old object
                    object.1 = core::mem::take(&mut groups);
                    if default_group.len() > 0 {
                        object
                            .1
                            .insert(String::new(), core::mem::take(&mut default_group));
                    }
                    selected_groups.clear();
                    if object.1.len() > 0 {
                        objects
                            .insert(object.0.unwrap_or_default(), core::mem::take(&mut object.1));
                    }

                    // Create new object
                    let name = terms
                        .map_while(|t| Some(t).filter(|t| util::name_is_valid(t)))
                        .next()
                        .ok_or_else(|| Error::ExpectedName(line_num))?
                        .to_string();
                    object.0 = Some(name);
                }
                _ => {}
            }
        }

        // Clean up old object
        object.1 = core::mem::take(&mut groups);
        if default_group.len() > 0 {
            object
                .1
                .insert(String::new(), core::mem::take(&mut default_group));
        }
        selected_groups.clear();
        if object.1.len() > 0 {
            objects.insert(object.0.unwrap_or_default(), core::mem::take(&mut object.1));
        }

        // Validate indices
        for (pos, uv, norm) in &vertices {
            if pos.get() - 1 >= positions.len() {
                return Err(Error::InvalidIndex(pos.get() as isize));
            }
            if let Some(uv) = *uv {
                if uv.get() - 1 >= uvs.len() {
                    return Err(Error::InvalidIndex(uv.get() as isize));
                }
            }
            if let Some(norm) = *norm {
                if norm.get() - 1 >= normals.len() {
                    return Err(Error::InvalidIndex(norm.get() as isize));
                }
            }
        }

        Ok(Self {
            buffers: Buffers {
                positions,
                uvs,
                normals,
                vertices,
            },
            objects,
        })
    }

    /// Write this [`Obj`] to a writer (something implementing [`std::io::Write`]) in OBJ format.
    pub fn write<W: Write>(&self, mut writer: W) -> Result<(), Error> {
        Ok(write!(writer, "{}", self)?)
    }

    /// Write this [`Obj`] to a file in OBJ format.
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), Error> {
        Ok(self.write(File::create(path)?)?)
    }

    /// Returns a specific [`Object`] by name.
    ///
    /// Note that if a name is not specified in the OBJ file, the name defaults to an empty string.
    pub fn object(&self, name: &str) -> Option<Object> {
        self.objects.get(name).map(|groups| Object {
            buffers: &self.buffers,
            groups,
        })
    }

    /// Returns an iterator over the [`Object`]s in this [`Obj`].
    pub fn objects(&self) -> impl ExactSizeIterator<Item = (&str, Object)> + Clone + '_ {
        self.objects.iter().map(move |(name, groups)| {
            (
                name.as_str(),
                Object {
                    buffers: &self.buffers,
                    groups,
                },
            )
        })
    }

    /// Returns an iterator over the [`Group`]s in this [`Obj`].
    pub fn groups(&self) -> impl Iterator<Item = (&String, Group)> + Clone + '_ {
        self.objects().map(|(_, object)| object.groups()).flatten()
    }

    /// Returns an iterator over the [`Polygon`]s in this [`Obj`].
    pub fn polygons(&self) -> impl Iterator<Item = Polygon> + Clone + '_ {
        self.groups().map(|(_, group)| group.polygons()).flatten()
    }

    /// Returns an iterator over the triangles in this [`Obj`].
    ///
    /// See [`Polygon::triangles`] for more information.
    pub fn triangles(&self) -> impl Iterator<Item = [Vertex; 3]> + Clone + '_ {
        self.polygons().map(|poly| poly.triangles()).flatten()
    }

    /// Returns an iterator over the vertices in this [`Obj`].
    pub fn vertices(&self) -> impl ExactSizeIterator<Item = Vertex> + Clone + '_ {
        self.buffers.vertices.iter().map(move |indices| Vertex {
            buffers: &self.buffers,
            indices: *indices,
        })
    }
}

impl Deref for Obj {
    type Target = Buffers;

    fn deref(&self) -> &Buffers {
        &self.buffers
    }
}

impl DerefMut for Obj {
    fn deref_mut(&mut self) -> &mut Buffers {
        &mut self.buffers
    }
}

impl fmt::Debug for Obj {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.debug_struct("Obj")
            .field("positions", &self.buffers.positions.len())
            .field("uvs", &self.buffers.uvs.len())
            .field("normals", &self.buffers.normals.len())
            .field("vertices", &self.buffers.vertices.len())
            .field("objects", &self.objects().collect::<Vec<_>>())
            .finish()
    }
}

impl fmt::Display for Obj {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for [x, y, z] in &self.buffers.positions {
            writeln!(f, "v {} {} {}", x, y, z)?;
        }
        for [u, v, w] in &self.buffers.uvs {
            writeln!(f, "vt {} {} {}", u, v, w)?;
        }
        for [x, y, z] in &self.buffers.normals {
            writeln!(f, "vn {} {} {}", x, y, z)?;
        }
        for (name, groups) in self.objects.iter() {
            if name.len() > 0 {
                writeln!(f, "o {}", name)?;
            }
            for (name, polys) in groups.iter() {
                if name.len() > 0 {
                    writeln!(f, "g {}", name)?;
                }
                for range in polys {
                    self.buffers.lookup(*range).display(f)?;
                    writeln!(f)?;
                }
            }
        }
        Ok(())
    }
}

/// An object defined in an OBJ.
#[derive(Copy, Clone)]
pub struct Object<'a> {
    buffers: &'a Buffers,
    groups: &'a HashMap<String, Vec<VertexRange>>,
}

impl<'a> Object<'a> {
    /// Returns a specific [`Group`] by name.
    ///
    /// Note that if a name is not specified in the OBJ file, the name defaults to an empty string.
    pub fn group(&self, name: &str) -> Option<Group<'a>> {
        self.groups.get(name).map(|polygons| Group {
            buffers: self.buffers,
            polygons: polygons.as_slice(),
        })
    }

    /// Returns an iterator over the [`Group`]s in this [`Object`].
    pub fn groups(&self) -> impl ExactSizeIterator<Item = (&'a String, Group<'a>)> + Clone + 'a {
        let buffers = self.buffers;
        self.groups.iter().map(move |(name, polygons)| {
            (
                name,
                Group {
                    buffers,
                    polygons: polygons.as_slice(),
                },
            )
        })
    }

    /// Returns an iterator over the [`Polygon`]s in this [`Object`].
    pub fn polygons(&self) -> impl Iterator<Item = Polygon<'a>> + Clone + 'a {
        self.groups().map(|(_, group)| group.polygons()).flatten()
    }

    /// Returns an iterator over the triangles in this [`Object`].
    ///
    /// See [`Polygon::triangles`] for more information.
    pub fn triangles(&self) -> impl Iterator<Item = [Vertex<'a>; 3]> + Clone + 'a {
        self.polygons().map(|poly| poly.triangles()).flatten()
    }
}

impl<'a> fmt::Debug for Object<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.debug_struct("Object")
            .field("groups", &self.groups().collect::<Vec<_>>())
            .finish()
    }
}

/// A group defined in an OBJ.
#[derive(Copy, Clone)]
pub struct Group<'a> {
    buffers: &'a Buffers,
    polygons: &'a [VertexRange],
}

impl<'a> Group<'a> {
    /// Returns a specific [`Polygon`] by index.
    pub fn polygon(&self, index: Index) -> Option<Polygon<'a>> {
        self.polygons
            .get(index)
            .map(|range| self.buffers.lookup(*range))
    }

    /// Returns an iterator over the [`Polygon`]s in this [`Group`].
    pub fn polygons(&self) -> impl ExactSizeIterator<Item = Polygon<'a>> + Clone + 'a {
        let buffers = self.buffers;
        self.polygons
            .iter()
            .map(move |range| buffers.lookup(*range))
    }

    /// Returns an iterator over the triangles in this [`Group`].
    ///
    /// See [`Polygon::triangles`] for more information.
    pub fn triangles(&self) -> impl Iterator<Item = [Vertex<'a>; 3]> + Clone + 'a {
        self.polygons().map(|poly| poly.triangles()).flatten()
    }
}

impl<'a> fmt::Debug for Group<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.debug_struct("Group")
            .field("polygons", &self.polygons().collect::<Vec<_>>())
            .finish()
    }
}

/// A polygon defined in an OBJ.
#[derive(Copy, Clone)]
pub struct Polygon<'a> {
    buffers: &'a Buffers,
    vertices: &'a [VertexIndices],
}

impl<'a> Polygon<'a> {
    /// Returns a specific [`Vertex`] by index.
    pub fn vertex(&self, index: usize) -> Option<Vertex<'a>> {
        self.vertices.get(index).map(|indices| Vertex {
            buffers: self.buffers,
            indices: *indices,
        })
    }

    /// Returns an iterator over the [`Vertex`]s in this [`Polygon`].
    pub fn vertices(&self) -> impl ExactSizeIterator<Item = Vertex<'a>> + Clone + 'a {
        let buffers = self.buffers;
        self.vertices.iter().map(move |indices| Vertex {
            buffers,
            indices: *indices,
        })
    }

    /// Returns an iterator over triangles in this [`Polygon`] by splitting up the polygon into smaller pieces.
    ///
    ///
    /// This function is useful when your application supports only triangles as input but the OBJ contains unusual
    /// polygons such as quads.
    ///
    /// The triangles produced will be arranged in a fan and will follow the winding order of the original
    /// polygon.
    ///
    /// This function assumes that:
    ///
    /// - The polygon is concave
    /// - The vertices of the polygon all lie in the same plane
    pub fn triangles(&self) -> impl ExactSizeIterator<Item = [Vertex<'a>; 3]> + Clone + 'a {
        let this = *self;
        (0..this.vertices.len().saturating_sub(2)).map(move |i| {
            [
                this.vertex(0).unwrap(),
                this.vertex(i + 1).unwrap(),
                this.vertex(i + 2).unwrap(),
            ]
        })
    }

    fn display(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "f")?;
        for v in self.vertices() {
            write!(f, " ")?;
            v.display(f)?;
        }
        Ok(())
    }
}

impl<'a> fmt::Debug for Polygon<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.debug_struct("Polygon")
            .field("vertices", &self.vertices().collect::<Vec<_>>())
            .finish()
    }
}

/// A vertex defined in an OBJ.
#[derive(Copy, Clone)]
pub struct Vertex<'a> {
    buffers: &'a Buffers,
    indices: VertexIndices,
}

impl<'a> Vertex<'a> {
    /// Returns the index of the vertex's position in the slice given by [`Buffers::positions`].
    ///
    /// Note that, unlike OBJ files themselves, this is zero-indexed.
    pub fn position_index(&self) -> Index {
        self.indices.0.get() - 1
    }

    /// Returns the position of this vertex.
    pub fn position(&self) -> [f32; 3] {
        self.buffers.positions[self.position_index()]
    }

    /// Returns the index of the vertex's texture coordinate, if it has one, in the slice given by [`Buffers::uvs`].
    ///
    /// Note that, unlike OBJ files themselves, this is zero-indexed.
    pub fn uv_index(&self) -> Option<Index> {
        self.indices.1.map(|idx| idx.get() - 1)
    }

    /// Returns the texture coordinate of this vertex, if it has one.
    pub fn uv(&self) -> Option<[f32; 3]> {
        Some(self.buffers.uvs[self.uv_index()?])
    }

    /// Returns the index of the vertex's normal, if it has one, in the slice given by [`Buffers::normals`].
    ///
    /// Note that, unlike OBJ files themselves, this is zero-indexed.
    pub fn normal_index(&self) -> Option<Index> {
        self.indices.2.map(|idx| idx.get() - 1)
    }

    /// Returns the normal of this vertex, if it has one.
    pub fn normal(&self) -> Option<[f32; 3]> {
        Some(self.buffers.normals[self.normal_index()?])
    }

    fn display(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, " {}/", self.indices.0)?;
        if let Some(uv) = self.indices.1 {
            write!(f, "{}", uv)?;
        }
        write!(f, "/")?;
        if let Some(norm) = self.indices.2 {
            write!(f, "{}", norm)?;
        }
        Ok(())
    }
}

impl<'a> fmt::Debug for Vertex<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.debug_struct("Vertex")
            .field("position", &self.position())
            .field("uv", &self.uv())
            .field("normal", &self.normal())
            .finish()
    }
}

/// Utilities relating to the OBJ format.
pub mod util {
    /// Determine whether a name (of either an object or a group) is valid.
    pub fn name_is_valid(name: &str) -> bool {
        name.chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '.' || c == '_')
    }
}

type VertexIndices = (NonZeroUsize, Option<NonZeroUsize>, Option<NonZeroUsize>);

/// A struct of buffers that may be indexed by [`Obj`] vertices.
#[derive(Clone, Default)]
pub struct Buffers {
    positions: Vec<[f32; 3]>,
    uvs: Vec<[f32; 3]>,
    normals: Vec<[f32; 3]>,
    vertices: Vec<VertexIndices>,
}

impl Buffers {
    fn lookup(&self, range: VertexRange) -> Polygon {
        Polygon {
            buffers: self,
            vertices: &self.vertices[range.start..range.end],
        }
    }

    /// Returns a reference to the position attributes contained within this [`Obj`].
    pub fn positions(&self) -> &[[f32; 3]] {
        &self.positions
    }

    /// Returns a reference to the texture coordinate attributes contained within this [`Obj`].
    pub fn uvs(&self) -> &[[f32; 3]] {
        &self.uvs
    }

    /// Returns a reference to the normal attributes contained within this [`Obj`].
    pub fn normals(&self) -> &[[f32; 3]] {
        &self.normals
    }

    /// Add a new position attribute to this [`Obj`], returning its index.
    pub fn add_position(&mut self, position: [f32; 3]) -> usize {
        let idx = self.positions.len();
        self.positions.push(position);
        idx
    }

    /// Add a new texture coordinate attribute to this [`Obj`], returning its index.
    pub fn add_uv(&mut self, uv: [f32; 3]) -> usize {
        let idx = self.uvs.len();
        self.uvs.push(uv);
        idx
    }

    /// Add a new normal attribute to this [`Obj`], returning its index.
    pub fn add_normal(&mut self, normal: [f32; 3]) -> usize {
        let idx = self.normals.len();
        self.normals.push(normal);
        idx
    }
}

#[derive(Copy, Clone)]
struct VertexRange {
    start: usize,
    end: usize,
}
