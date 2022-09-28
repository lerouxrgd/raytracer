use std::cmp;
use std::f32::consts::PI;

use lazy_static::lazy_static;
use ordered_float::OrderedFloat;
use parking_lot::RwLock;
use serde::{Deserialize, Serialize};
use slotmap::{new_key_type, SlotMap};

use crate::canvas::Canvas;
use crate::matrices::Matrix;
use crate::shapes::Shape;
use crate::tuples::{Color, Point, Vector};

#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Pattern {
    Striped(Striped),
    Gradient(Gradient),
    Ring(Ring),
    Checker(Checker),
    XyzRgb(XyzRgb),
    TextureMap(TextureMap),
    CubeMap(CubeMap),
}

impl Pattern {
    pub fn with_transform<T: Into<Matrix<4, 4>>>(mut self, transform: T) -> Self {
        self.set_transform(transform.into());
        self
    }

    pub fn get_transform(&self) -> Matrix<4, 4> {
        match self {
            &Self::Striped(Striped { transform, .. })
            | &Self::Gradient(Gradient { transform, .. })
            | &Self::Ring(Ring { transform, .. })
            | &Self::Checker(Checker { transform, .. })
            | &Self::XyzRgb(XyzRgb { transform, .. })
            | &Self::TextureMap(TextureMap { transform, .. })
            | &Self::CubeMap(CubeMap { transform, .. }) => transform,
        }
    }

    pub fn set_transform<T: Into<Matrix<4, 4>>>(&mut self, t: T) {
        match self {
            Self::Striped(Striped {
                ref mut transform, ..
            })
            | Self::Gradient(Gradient {
                ref mut transform, ..
            })
            | Self::Ring(Ring {
                ref mut transform, ..
            })
            | Self::Checker(Checker {
                ref mut transform, ..
            })
            | Self::XyzRgb(XyzRgb {
                ref mut transform, ..
            })
            | Self::TextureMap(TextureMap {
                ref mut transform, ..
            })
            | Self::CubeMap(CubeMap {
                ref mut transform, ..
            }) => *transform = t.into(),
        }
    }

    pub fn pattern_at_shape(&self, shape: Shape, world_point: Point) -> Color {
        let object_point = shape.world_to_object(world_point);
        let pattern_point = self.get_transform().inverse().unwrap() * object_point;
        match self {
            Self::Striped(s) => s.pattern_at(pattern_point),
            Self::Gradient(g) => g.pattern_at(pattern_point),
            Self::Ring(r) => r.pattern_at(pattern_point),
            Self::Checker(c) => c.pattern_at(pattern_point),
            Self::XyzRgb(xyz) => xyz.pattern_at(pattern_point),
            Self::TextureMap(tm) => tm.pattern_at(pattern_point),
            Self::CubeMap(cm) => cm.pattern_at(pattern_point),
        }
    }
}

/// Stripes based on x value
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Striped {
    pub a: Color,
    pub b: Color,
    pub transform: Matrix<4, 4>,
}

impl Striped {
    pub fn new(a: Color, b: Color) -> Self {
        Self {
            a,
            b,
            transform: Matrix::identity(),
        }
    }

    pub fn with_transform<T: Into<Matrix<4, 4>>>(mut self, transform: T) -> Self {
        self.transform = transform.into();
        self
    }

    pub fn pattern_at(&self, pattern_point: Point) -> Color {
        if pattern_point.x().floor() % 2. == 0. {
            self.a
        } else {
            self.b
        }
    }
}

impl From<Striped> for Pattern {
    fn from(s: Striped) -> Self {
        Pattern::Striped(s)
    }
}

impl From<Striped> for Option<Pattern> {
    fn from(s: Striped) -> Self {
        Some(Pattern::Striped(s))
    }
}

/// Linearly blend two colors
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Gradient {
    pub a: Color,
    pub b: Color,
    pub transform: Matrix<4, 4>,
}

impl Gradient {
    pub fn new(a: Color, b: Color) -> Self {
        Self {
            a,
            b,
            transform: Matrix::identity(),
        }
    }

    pub fn with_transform<T: Into<Matrix<4, 4>>>(mut self, transform: T) -> Self {
        self.transform = transform.into();
        self
    }

    pub fn pattern_at(&self, pattern_point: Point) -> Color {
        let distance = self.b - self.a;
        let fraction = pattern_point.x() - pattern_point.x().floor();
        self.a + distance * fraction
    }
}

impl From<Gradient> for Pattern {
    fn from(s: Gradient) -> Self {
        Pattern::Gradient(s)
    }
}

impl From<Gradient> for Option<Pattern> {
    fn from(s: Gradient) -> Self {
        Some(Pattern::Gradient(s))
    }
}

/// An xz ring
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ring {
    pub a: Color,
    pub b: Color,
    pub transform: Matrix<4, 4>,
}

impl Ring {
    pub fn new(a: Color, b: Color) -> Self {
        Self {
            a,
            b,
            transform: Matrix::identity(),
        }
    }

    pub fn with_transform<T: Into<Matrix<4, 4>>>(mut self, transform: T) -> Self {
        self.transform = transform.into();
        self
    }

    pub fn pattern_at(&self, pattern_point: Point) -> Color {
        if (pattern_point.x().powi(2) + pattern_point.z().powi(2))
            .sqrt()
            .floor()
            % 2.
            == 0.
        {
            self.a
        } else {
            self.b
        }
    }
}

impl From<Ring> for Pattern {
    fn from(s: Ring) -> Self {
        Pattern::Ring(s)
    }
}

impl From<Ring> for Option<Pattern> {
    fn from(s: Ring) -> Self {
        Some(Pattern::Ring(s))
    }
}

/// Repeating pattern of squares (adjacent squares have different colors)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Checker {
    pub a: Color,
    pub b: Color,
    pub transform: Matrix<4, 4>,
}

impl Checker {
    pub fn new(a: Color, b: Color) -> Self {
        Self {
            a,
            b,
            transform: Matrix::identity(),
        }
    }

    pub fn with_transform<T: Into<Matrix<4, 4>>>(mut self, transform: T) -> Self {
        self.transform = transform.into();
        self
    }

    pub fn pattern_at(&self, pattern_point: Point) -> Color {
        if (pattern_point.x().floor() + pattern_point.y().floor() + pattern_point.z().floor()) % 2.
            == 0.
        {
            self.a
        } else {
            self.b
        }
    }
}

impl From<Checker> for Pattern {
    fn from(s: Checker) -> Self {
        Pattern::Checker(s)
    }
}

impl From<Checker> for Option<Pattern> {
    fn from(s: Checker) -> Self {
        Some(Pattern::Checker(s))
    }
}

/// Coordinates are the color components
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct XyzRgb {
    pub transform: Matrix<4, 4>,
}

impl XyzRgb {
    pub fn new() -> Self {
        Self {
            transform: Matrix::identity(),
        }
    }

    pub fn with_transform<T: Into<Matrix<4, 4>>>(mut self, transform: T) -> Self {
        self.transform = transform.into();
        self
    }

    pub fn pattern_at(&self, pattern_point: Point) -> Color {
        Color::new(pattern_point.x(), pattern_point.y(), pattern_point.z())
    }
}

impl Default for XyzRgb {
    fn default() -> Self {
        Self::new()
    }
}

impl From<XyzRgb> for Pattern {
    fn from(xyz: XyzRgb) -> Self {
        Pattern::XyzRgb(xyz)
    }
}

impl From<XyzRgb> for Option<Pattern> {
    fn from(xyz: XyzRgb) -> Self {
        Some(Pattern::XyzRgb(xyz))
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TextureMap {
    pub transform: Matrix<4, 4>,
    pub uv_pattern: UvPattern,
    pub uv_mapping: UvMapping,
}

impl TextureMap {
    pub fn new(uv_pattern: UvPattern, uv_mapping: UvMapping) -> Self {
        Self {
            transform: Matrix::identity(),
            uv_pattern,
            uv_mapping,
        }
    }

    pub fn with_transform<T: Into<Matrix<4, 4>>>(mut self, transform: T) -> Self {
        self.transform = transform.into();
        self
    }

    pub fn pattern_at(&self, pattern_point: Point) -> Color {
        let (u, v) = self.uv_mapping.map_point(pattern_point);
        self.uv_pattern.uv_pattern_at(u, v)
    }
}

impl From<TextureMap> for Pattern {
    fn from(tm: TextureMap) -> Self {
        Pattern::TextureMap(tm)
    }
}

impl From<TextureMap> for Option<Pattern> {
    fn from(tm: TextureMap) -> Self {
        Some(Pattern::TextureMap(tm))
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CubeMap {
    pub transform: Matrix<4, 4>,
    pub left: UvPattern,
    pub front: UvPattern,
    pub right: UvPattern,
    pub back: UvPattern,
    pub up: UvPattern,
    pub down: UvPattern,
}

impl CubeMap {
    pub fn new(
        left: UvPattern,
        front: UvPattern,
        right: UvPattern,
        back: UvPattern,
        up: UvPattern,
        down: UvPattern,
    ) -> Self {
        Self {
            transform: Matrix::identity(),
            left,
            front,
            right,
            back,
            up,
            down,
        }
    }

    pub fn with_transform<T: Into<Matrix<4, 4>>>(mut self, transform: T) -> Self {
        self.transform = transform.into();
        self
    }

    pub fn pattern_at(&self, pattern_point: Point) -> Color {
        let face = CubeFace::from(pattern_point);
        let (u, v) = face.uv(pattern_point);
        match face {
            CubeFace::Left => self.left.uv_pattern_at(u, v),
            CubeFace::Front => self.front.uv_pattern_at(u, v),
            CubeFace::Right => self.right.uv_pattern_at(u, v),
            CubeFace::Back => self.back.uv_pattern_at(u, v),
            CubeFace::Up => self.up.uv_pattern_at(u, v),
            CubeFace::Down => self.down.uv_pattern_at(u, v),
        }
    }
}

impl From<CubeMap> for Pattern {
    fn from(cm: CubeMap) -> Self {
        Pattern::CubeMap(cm)
    }
}

impl From<CubeMap> for Option<Pattern> {
    fn from(cm: CubeMap) -> Self {
        Some(Pattern::CubeMap(cm))
    }
}

////////////////////////////////////////////////////////////////////////////////////////

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum UvPattern {
    Checker(UvChecker),
    AlignCheck(UvAlignCheck),
    Image(UvImage),
}

impl UvPattern {
    pub fn uv_pattern_at(&self, u: f32, v: f32) -> Color {
        match self {
            Self::Checker(p) => p.uv_pattern_at(u, v),
            Self::AlignCheck(p) => p.uv_pattern_at(u, v),
            Self::Image(p) => p.uv_pattern_at(u, v),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UvChecker {
    pub width: f32,
    pub height: f32,
    pub a: Color,
    pub b: Color,
}

impl UvChecker {
    pub fn new(width: f32, height: f32, a: Color, b: Color) -> Self {
        Self {
            width,
            height,
            a,
            b,
        }
    }

    pub fn uv_pattern_at(&self, u: f32, v: f32) -> Color {
        let u2 = (u * self.width).floor();
        let v2 = (v * self.height).floor();
        if (u2 + v2) % 2. == 0. {
            self.a
        } else {
            self.b
        }
    }
}

impl From<UvChecker> for UvPattern {
    fn from(uv_checker: UvChecker) -> Self {
        UvPattern::Checker(uv_checker)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UvAlignCheck {
    main: Color,
    ul: Color,
    ur: Color,
    bl: Color,
    br: Color,
}

impl From<UvAlignCheck> for UvPattern {
    fn from(uv_align_check: UvAlignCheck) -> Self {
        UvPattern::AlignCheck(uv_align_check)
    }
}

impl UvAlignCheck {
    pub fn new(main: Color, ul: Color, ur: Color, bl: Color, br: Color) -> Self {
        Self {
            main,
            ul,
            ur,
            bl,
            br,
        }
    }

    pub fn uv_pattern_at(&self, u: f32, v: f32) -> Color {
        if v > 0.8 {
            if u < 0.2 {
                return self.ul;
            } else if u > 0.8 {
                return self.ur;
            }
        } else if v < 0.2 {
            if u < 0.2 {
                return self.bl;
            } else if u > 0.8 {
                return self.br;
            }
        }
        self.main
    }
}

enum CubeFace {
    Left,
    Right,
    Front,
    Back,
    Up,
    Down,
}

impl From<Point> for CubeFace {
    fn from(p: Point) -> Self {
        let abs_x = OrderedFloat(p.x().abs());
        let abs_y = OrderedFloat(p.y().abs());
        let abs_z = OrderedFloat(p.z().abs());
        let coord: f32 = cmp::max(cmp::max(abs_x, abs_y), abs_z).into();
        if coord == p.x() {
            Self::Right
        } else if coord == -p.x() {
            Self::Left
        } else if coord == p.y() {
            Self::Up
        } else if coord == -p.y() {
            Self::Down
        } else if coord == p.z() {
            Self::Front
        } else {
            Self::Back
        }
    }
}

impl CubeFace {
    pub fn uv(&self, p: Point) -> (f32, f32) {
        match self {
            Self::Front => {
                let u = ((p.x() + 1.) % 2.) / 2.;
                let v = ((p.y() + 1.) % 2.) / 2.;
                (u, v)
            }
            Self::Back => {
                let u = ((1. - p.x()) % 2.) / 2.;
                let v = ((p.y() + 1.) % 2.) / 2.;
                (u, v)
            }
            Self::Up => {
                let u = ((p.x() + 1.) % 2.) / 2.;
                let v = ((1. - p.z()) % 2.) / 2.;
                (u, v)
            }
            Self::Down => {
                let u = ((p.x() + 1.) % 2.) / 2.;
                let v = ((p.z() + 1.) % 2.) / 2.;
                (u, v)
            }
            Self::Left => {
                let u = ((p.z() + 1.) % 2.) / 2.;
                let v = ((p.y() + 1.) % 2.) / 2.;
                (u, v)
            }
            Self::Right => {
                let u = ((1. - p.z()) % 2.) / 2.;
                let v = ((p.y() + 1.) % 2.) / 2.;
                (u, v)
            }
        }
    }
}

lazy_static! {
    static ref UV_IMAGES: RwLock<SlotMap<UvImage, Canvas>> = RwLock::new(SlotMap::with_key());
}

new_key_type! {
    pub struct UvImage;
}

impl UvImage {
    pub fn new(canvas: Canvas) -> Self {
        UV_IMAGES.write().insert(canvas)
    }

    pub fn uv_pattern_at(&self, u: f32, v: f32) -> Color {
        // flip v over so it matches the image layout, with y at the top
        let v = 1. - v;

        let images = UV_IMAGES.read();
        let canvas = images.get(*self).unwrap();
        let x = u * (canvas.width() as f32 - 1.);
        let y = v * (canvas.height() as f32 - 1.);

        canvas.pixel_at(x.round() as u32, y.round() as u32)
    }
}

impl From<UvImage> for UvPattern {
    fn from(uv_image: UvImage) -> Self {
        UvPattern::Image(uv_image)
    }
}

////////////////////////////////////////////////////////////////////////////////////////

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum UvMapping {
    Spherical,
    Planar,
    Cylindrical,
}

impl UvMapping {
    pub fn map_point(&self, p: Point) -> (f32, f32) {
        match self {
            Self::Spherical => spherical_map(p),
            Self::Planar => planar_map(p),
            Self::Cylindrical => cylindrical_map(p),
        }
    }
}

fn spherical_map(p: Point) -> (f32, f32) {
    // compute the azimuthal angle
    // -π < theta <= π
    // angle increases clockwise as viewed from above,
    // which is opposite of what we want, but we'll fix it later.
    let theta = p.x().atan2(p.z());

    // v is the vector pointing from the sphere's origin (the world origin)
    // to the point, which will also happen to be exactly equal to the sphere's
    // radius.
    let v = Vector::new(p.x(), p.y(), p.z());
    let radius = v.magnitude();

    // compute the polar angle
    // 0 <= phi <= π
    let phi = (p.y() / radius).acos();

    // -0.5 < raw_u <= 0.5
    let raw_u = theta / (2. * PI);

    // 0 <= u < 1
    // here's also where we fix the direction of u. Subtract it from 1,
    // so that it increases counterclockwise as viewed from above.
    let u = 1. - (raw_u + 0.5);

    // we want v to be 0 at the south pole of the sphere,
    // and 1 at the north pole, so we have to "flip it over"
    // by subtracting it from 1.
    let v = 1. - phi / PI;

    (u, v)
}

fn planar_map(p: Point) -> (f32, f32) {
    (p.x().rem_euclid(1.), p.z().rem_euclid(1.))
}

fn cylindrical_map(p: Point) -> (f32, f32) {
    // compute the azimuthal angle, same as with spherical_map()
    let theta = p.x().atan2(p.z());
    let raw_u = theta / (2. * PI);
    let u = 1. - (raw_u + 0.5);
    //let v go from 0 to 1 between whole units of y
    let v = p.y().rem_euclid(1.);
    (u, v)
}

////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;
    use crate::shapes::Sphere;
    use crate::transformations::{scaling, translation};

    #[test]
    fn striped_basics() {
        let p = Striped::new(Color::white(), Color::black());
        assert!(p.pattern_at(Point::new(0., 0., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(0., 1., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(0., 2., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(0., 0., 1.)) == Color::white());
        assert!(p.pattern_at(Point::new(0., 0., 2.)) == Color::white());
        assert!(p.pattern_at(Point::new(0.9, 0., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(1., 0., 0.)) == Color::black());
        assert!(p.pattern_at(Point::new(-0.1, 0., 0.)) == Color::black());
        assert!(p.pattern_at(Point::new(-1., 0., 0.)) == Color::black());
        assert!(p.pattern_at(Point::new(-1.1, 0., 0.)) == Color::white());
    }

    #[test]
    fn striped_transform() {
        let shape = Sphere::default().with_transform(scaling(2., 2., 2.));
        let p: Pattern = Striped::new(Color::white(), Color::black()).into();
        let c = p.pattern_at_shape(shape.into(), Point::new(1.5, 0., 0.));
        assert!(c == Color::white());

        let shape = Sphere::default();
        let p: Pattern = Striped::new(Color::white(), Color::black())
            .with_transform(scaling(2., 2., 2.))
            .into();
        let c = p.pattern_at_shape(shape.into(), Point::new(1.5, 0., 0.));
        assert!(c == Color::white());

        let shape = Sphere::default().with_transform(scaling(2., 2., 2.));
        let p: Pattern = Striped::new(Color::white(), Color::black())
            .with_transform(translation(0.5, 0., 0.))
            .into();
        let c = p.pattern_at_shape(shape.into(), Point::new(2.5, 0., 0.));
        assert!(c == Color::white());
    }

    #[test]
    fn gradient_basics() {
        let p = Gradient::new(Color::white(), Color::black());
        assert!(p.pattern_at(Point::new(0., 0., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(0.25, 0., 0.)) == Color::new(0.75, 0.75, 0.75));
        assert!(p.pattern_at(Point::new(0.5, 0., 0.)) == Color::new(0.5, 0.5, 0.5));
        assert!(p.pattern_at(Point::new(0.75, 0., 0.)) == Color::new(0.25, 0.25, 0.25));
    }

    #[test]
    fn ring_basics() {
        let p = Ring::new(Color::white(), Color::black());
        assert!(p.pattern_at(Point::new(0., 0., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(1., 0., 0.)) == Color::black());
        assert!(p.pattern_at(Point::new(0., 0., 1.)) == Color::black());
        assert!(p.pattern_at(Point::new(0.708, 0., 0.708)) == Color::black());
    }

    #[test]
    fn checker_basics() {
        let p = Checker::new(Color::white(), Color::black());
        assert!(p.pattern_at(Point::new(0., 0., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(0.99, 0., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(1.01, 0., 0.)) == Color::black());
        assert!(p.pattern_at(Point::new(0., 0., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(0., 0.99, 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(0., 1.01, 0.)) == Color::black());
        assert!(p.pattern_at(Point::new(0., 0., 0.99)) == Color::white());
        assert!(p.pattern_at(Point::new(0., 0., 1.01)) == Color::black());
    }

    #[test]
    fn uv_checkers_basics() {
        let pattern = UvChecker::new(2., 2., Color::black(), Color::white());
        let cases = vec![
            (0.0, 0.0, Color::black()),
            (0.5, 0.0, Color::white()),
            (0.0, 0.5, Color::white()),
            (0.5, 0.5, Color::black()),
            (1.0, 1.0, Color::black()),
        ];
        for (u, v, res) in cases {
            assert!(pattern.uv_pattern_at(u, v) == res);
        }
    }

    #[test]
    fn spherical_mapping() {
        let cases = vec![
            (Point::new(0., 0., -1.), 0.0, 0.5),
            (Point::new(1., 0., 0.), 0.25, 0.5),
            (Point::new(0., 0., 1.), 0.5, 0.5),
            (Point::new(-1., 0., 0.), 0.75, 0.5),
            (Point::new(0., 1., 0.), 0.5, 1.0),
            (Point::new(0., -1., 0.), 0.5, 0.0),
            (
                Point::new(f32::sqrt(2.) / 2., f32::sqrt(2.) / 2., 0.),
                0.25,
                0.75,
            ),
        ];
        for (p, u, v) in cases {
            assert!(spherical_map(p) == (u, v));
        }
    }

    #[test]
    fn texture_map_spherical() {
        let cases = vec![
            (Point::new(0.4315, 0.4670, 0.7719), Color::white()),
            (Point::new(-0.9654, 0.2552, -0.0534), Color::black()),
            (Point::new(0.1039, 0.7090, 0.6975), Color::white()),
            (Point::new(-0.4986, -0.7856, -0.3663), Color::black()),
            (Point::new(-0.0317, -0.9395, 0.3411), Color::black()),
            (Point::new(0.4809, -0.7721, 0.4154), Color::black()),
            (Point::new(0.0285, -0.9612, -0.2745), Color::black()),
            (Point::new(-0.5734, -0.2162, -0.7903), Color::white()),
            (Point::new(0.7688, -0.1470, 0.6223), Color::black()),
            (Point::new(-0.7652, 0.2175, 0.6060), Color::black()),
        ];
        let tm = TextureMap::new(
            UvChecker::new(16., 8., Color::black(), Color::white()).into(),
            UvMapping::Spherical,
        );
        for (p, color) in cases {
            assert!(tm.pattern_at(p) == color);
        }
    }

    #[test]
    fn planar_mapping() {
        let cases = vec![
            (Point::new(0.25, 0., 0.5), 0.25, 0.5),
            (Point::new(0.25, 0., -0.25), 0.25, 0.75),
            (Point::new(0.25, 0.5, -0.25), 0.25, 0.75),
            (Point::new(1.25, 0., 0.5), 0.25, 0.5),
            (Point::new(0.25, 0., -1.75), 0.25, 0.25),
            (Point::new(1., 0., -1.), 0.0, 0.0),
            (Point::new(0., 0., 0.), 0.0, 0.0),
        ];
        for (p, u, v) in cases {
            assert!(planar_map(p) == (u, v));
        }
    }

    #[test]
    fn cylindrical_mapping() {
        let cases = vec![
            (Point::new(0., 0., -1.), 0.0, 0.0),
            (Point::new(0., 0.5, -1.), 0.0, 0.5),
            (Point::new(0., 1., -1.), 0.0, 0.0),
            (Point::new(0.70711, 0.5, -0.70711), 0.125, 0.5),
            (Point::new(1., 0.5, 0.), 0.25, 0.5),
            (Point::new(0.70711, 0.5, 0.70711), 0.375, 0.5),
            (Point::new(0., -0.25, 1.), 0.5, 0.75),
            (Point::new(-0.70711, 0.5, 0.70711), 0.625, 0.5),
            (Point::new(-1., 1.25, 0.), 0.75, 0.25),
            (Point::new(-0.70711, 0.5, -0.70711), 0.875, 0.5),
        ];
        for (p, u, v) in cases {
            assert!(cylindrical_map(p) == (u, v));
        }
    }

    #[test]
    fn uv_align_check_basics() {
        let main = Color::white();
        let ul = Color::new(1., 0., 0.);
        let ur = Color::new(1., 1., 0.);
        let bl = Color::new(0., 1., 0.);
        let br = Color::new(0., 1., 1.);
        let pattern = UvAlignCheck::new(main, ul, ur, bl, br);
        let cases = vec![
            (0.5, 0.5, main),
            (0.1, 0.9, ul),
            (0.9, 0.9, ur),
            (0.1, 0.1, bl),
            (0.9, 0.1, br),
        ];
        for (u, v, res) in cases {
            assert!(pattern.uv_pattern_at(u, v) == res);
        }
    }

    #[test]
    fn uv_cube_faces() {
        let cases = vec![
            (Point::new(-0.5, 0.5, 1.), 0.25, 0.75),
            (Point::new(0.5, -0.5, 1.), 0.75, 0.25),
        ];
        for (p, u, v) in cases {
            assert!(CubeFace::Front.uv(p) == (u, v));
        }

        let cases = vec![
            (Point::new(0.5, 0.5, -1.), 0.25, 0.75),
            (Point::new(-0.5, -0.5, -1.), 0.75, 0.25),
        ];
        for (p, u, v) in cases {
            assert!(CubeFace::Back.uv(p) == (u, v));
        }

        let cases = vec![
            (Point::new(-1., 0.5, -0.5), 0.25, 0.75),
            (Point::new(-1., -0.5, 0.5), 0.75, 0.25),
        ];
        for (p, u, v) in cases {
            assert!(CubeFace::Left.uv(p) == (u, v));
        }

        let cases = vec![
            (Point::new(1., 0.5, 0.5), 0.25, 0.75),
            (Point::new(1., -0.5, -0.5), 0.75, 0.25),
        ];
        for (p, u, v) in cases {
            assert!(CubeFace::Right.uv(p) == (u, v));
        }

        let cases = vec![
            (Point::new(-0.5, 1., -0.5), 0.25, 0.75),
            (Point::new(0.5, 1., 0.5), 0.75, 0.25),
        ];
        for (p, u, v) in cases {
            assert!(CubeFace::Up.uv(p) == (u, v));
        }

        let cases = vec![
            (Point::new(-0.5, -1., 0.5), 0.25, 0.75),
            (Point::new(0.5, -1., -0.5), 0.75, 0.25),
        ];
        for (p, u, v) in cases {
            assert!(CubeFace::Down.uv(p) == (u, v));
        }
    }

    #[test]
    fn uv_image_basics() {
        let ppm_data = b"P3
10 10
10
0 0 0  1 1 1  2 2 2  3 3 3  4 4 4  5 5 5  6 6 6  7 7 7  8 8 8  9 9 9
1 1 1  2 2 2  3 3 3  4 4 4  5 5 5  6 6 6  7 7 7  8 8 8  9 9 9  0 0 0
2 2 2  3 3 3  4 4 4  5 5 5  6 6 6  7 7 7  8 8 8  9 9 9  0 0 0  1 1 1
3 3 3  4 4 4  5 5 5  6 6 6  7 7 7  8 8 8  9 9 9  0 0 0  1 1 1  2 2 2
4 4 4  5 5 5  6 6 6  7 7 7  8 8 8  9 9 9  0 0 0  1 1 1  2 2 2  3 3 3
5 5 5  6 6 6  7 7 7  8 8 8  9 9 9  0 0 0  1 1 1  2 2 2  3 3 3  4 4 4
6 6 6  7 7 7  8 8 8  9 9 9  0 0 0  1 1 1  2 2 2  3 3 3  4 4 4  5 5 5
7 7 7  8 8 8  9 9 9  0 0 0  1 1 1  2 2 2  3 3 3  4 4 4  5 5 5  6 6 6
8 8 8  9 9 9  0 0 0  1 1 1  2 2 2  3 3 3  4 4 4  5 5 5  6 6 6  7 7 7
9 9 9  0 0 0  1 1 1  2 2 2  3 3 3  4 4 4  5 5 5  6 6 6  7 7 7  8 8 8";

        let reader = Cursor::new(ppm_data);
        let canvas = Canvas::from_ppm(reader).unwrap();
        let pattern = UvImage::new(canvas);
        let cases = vec![
            (0., 0., Color::new(0.9, 0.9, 0.9)),
            (0.3, 0., Color::new(0.2, 0.2, 0.2)),
            (0.6, 0.3, Color::new(0.1, 0.1, 0.1)),
            (1., 1., Color::new(0.9, 0.9, 0.9)),
        ];
        for (u, v, res) in cases {
            assert!(pattern.uv_pattern_at(u, v) == res);
        }
    }
}
