#![allow(clippy::large_enum_variant)]

use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufReader, Write};
use std::path::PathBuf;

use serde::{de, Deserialize, Deserializer};

use crate::camera::Camera;
use crate::csg::{Csg, CsgChild, CsgOp};
use crate::groups::Group;
use crate::lights::{AreaLight, Light, PointLight};
use crate::materials::Material;
use crate::obj::parse_obj;
use crate::patterns::{Checker, Gradient, Pattern, Ring, Striped, XyzRgb};
use crate::shapes::{Cone, Cube, Cylinder, Plane, Shape, SmoothTriangle, Sphere, Triangle};
use crate::transformations::{view_transform, Transform};
use crate::world::World;

#[derive(Debug, Clone)]
pub struct Scene {
    instructions: Vec<Instruction>,
}

impl<'de> Deserialize<'de> for Scene {
    fn deserialize<D>(deserializer: D) -> Result<Scene, D::Error>
    where
        D: Deserializer<'de>,
    {
        let instructions = Vec::<Instruction>::deserialize(deserializer)?;
        if !instructions
            .iter()
            .any(|i| matches!(i, Instruction::Add(Add::AddGear(AddGear::Camera { .. }))))
        {
            return Err(de::Error::custom("Missing camera"));
        }
        if !instructions.iter().any(|i| {
            matches!(
                i,
                Instruction::Add(Add::AddGear(AddGear::PointLight { .. }))
            ) || matches!(i, Instruction::Add(Add::AddGear(AddGear::AreaLight { .. })))
        }) {
            return Err(de::Error::custom("Missing light"));
        }
        Ok(Scene { instructions })
    }
}

impl Scene {
    pub fn render(&self, out: impl Write, obj_files: &[PathBuf]) -> Result<(), Box<dyn Error>> {
        let mut camera: Option<Camera> = None;
        let mut light: Option<Light> = None;
        let mut define_transforms = HashMap::<String, Vec<TransformSpec>>::new();
        let mut define_materials = HashMap::<String, Vec<MaterialSpec>>::new();
        let mut shapes: Vec<Shape> = Vec::new();
        let mut groups: Vec<Group> = Vec::new();
        let mut csgs: Vec<Csg> = Vec::new();

        for instruction in &self.instructions {
            match instruction {
                &Instruction::Add(Add::AddGear(AddGear::Camera {
                    width,
                    height,
                    field_of_view,
                    from,
                    to,
                    up,
                })) => {
                    camera =
                        Some(
                            Camera::new(width, height, field_of_view)
                                .with_transform(view_transform(from.into(), to.into(), up.into())),
                        )
                }

                &Instruction::Add(Add::AddGear(AddGear::PointLight { at, intensity })) => {
                    light = Some(PointLight::new(at.into(), intensity.into()).into())
                }

                &Instruction::Add(Add::AddGear(AddGear::AreaLight {
                    corner,
                    uvec,
                    vvec,
                    usteps,
                    vsteps,
                    intensity,
                })) => {
                    light = Some(
                        AreaLight::new(
                            corner.into(),
                            uvec.into(),
                            usteps,
                            vvec.into(),
                            vsteps,
                            intensity.into(),
                            #[cfg(test)]
                            vec![0.5],
                        )
                        .into(),
                    )
                }

                Instruction::Define(Define {
                    define,
                    extend,
                    transform,
                    material,
                }) => {
                    if let Some(specs) = transform {
                        let mut final_specs = vec![];
                        if let Some(definitions) = extend {
                            for definition in definitions {
                                if let Some(base_specs) = define_transforms.get(definition) {
                                    final_specs.extend(base_specs.clone());
                                }
                            }
                        }
                        final_specs.extend(specs.clone());
                        define_transforms.insert(define.clone(), final_specs);
                    }
                    if let Some(spec) = material {
                        let mut final_specs = vec![];
                        if let Some(extensions) = extend {
                            for ext in extensions {
                                if let Some(base_specs) = define_materials.get(ext) {
                                    final_specs.extend(base_specs.clone());
                                }
                            }
                        }
                        final_specs.push(spec.clone());
                        define_materials.insert(define.clone(), final_specs);
                    }
                }

                Instruction::Add(Add::AddShape(a @ AddShape::Sphere { .. }))
                | Instruction::Add(Add::AddShape(a @ AddShape::Plane { .. }))
                | Instruction::Add(Add::AddShape(a @ AddShape::Cube { .. }))
                | Instruction::Add(Add::AddShape(a @ AddShape::Cylinder { .. }))
                | Instruction::Add(Add::AddShape(a @ AddShape::Cone { .. }))
                | Instruction::Add(Add::AddShape(a @ AddShape::Triangle { .. }))
                | Instruction::Add(Add::AddShape(a @ AddShape::SmoothTriangle { .. })) => {
                    shapes.push(a.make_shape(&define_transforms, &define_materials));
                }

                Instruction::Add(Add::AddShape(a @ AddShape::Group { .. })) => {
                    groups.push(a.make_group(&define_transforms, &define_materials));
                }

                Instruction::Add(Add::AddCsg(a @ AddCsg::Csg { .. })) => {
                    csgs.push(a.make_csg(&define_transforms, &define_materials));
                }
            }
        }

        for obj_file in obj_files {
            let file_name = obj_file
                .file_name()
                .map(|name| name.to_string_lossy().to_string())
                .ok_or("Invalid obj file")?;

            let mut final_transform = Transform::default();
            if let Some(specs) = define_transforms.get(&file_name) {
                for t in specs {
                    final_transform = t.update(final_transform);
                }
            }

            let reader = BufReader::new(File::open(obj_file)?);
            let obj_group = parse_obj(reader)?.with_transform(final_transform);

            groups.push(obj_group);
        }

        let world = World {
            shapes,
            groups,
            csgs,
            light: light.unwrap(),
        };
        let canvas = camera.unwrap().render(&world);

        canvas.to_ppm(out);
        Ok(())
    }
}

#[derive(Debug, Clone, Deserialize)]
#[serde(untagged)]
pub enum Instruction {
    Add(Add),
    Define(Define),
}

#[derive(Debug, Clone, Deserialize)]
#[serde(untagged)]
pub enum Add {
    AddGear(AddGear),
    AddShape(AddShape),
    AddCsg(AddCsg),
}

#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "add")]
#[serde(rename_all = "kebab-case")]
pub enum AddGear {
    PointLight {
        at: [f32; 3],
        intensity: [f32; 3],
    },
    AreaLight {
        corner: [f32; 3],
        uvec: [f32; 3],
        vvec: [f32; 3],
        usteps: usize,
        vsteps: usize,
        intensity: [f32; 3],
    },
    #[serde(rename_all = "kebab-case")]
    Camera {
        width: u32,
        height: u32,
        #[serde(deserialize_with = "deser_meval")]
        field_of_view: f32,
        from: [f32; 3],
        up: [f32; 3],
        to: [f32; 3],
    },
}

fn deser_meval<'de, D>(deserializer: D) -> Result<f32, D::Error>
where
    D: Deserializer<'de>,
{
    use serde_yaml::Value;

    let expr = Value::deserialize(deserializer)?;
    match expr {
        Value::String(expr) => meval::eval_str(&expr)
            .map(|v| v as f32)
            .map_err(de::Error::custom),
        Value::Number(num) => Ok(to_f32(&num)),
        _ => Err(de::Error::custom(format!(
            "Invalid math expression: {expr:?}"
        ))),
    }
}

fn to_f32(num: &serde_yaml::Number) -> f32 {
    if let Some(num) = num.as_f64() {
        num as f32
    } else if let Some(num) = num.as_i64() {
        num as f32
    } else if let Some(num) = num.as_u64() {
        num as f32
    } else {
        unreachable!()
    }
}

#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "add")]
#[serde(rename_all = "kebab-case")]
pub enum AddShape {
    Sphere {
        extend: Option<Vec<String>>,
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
        material: Option<MaterialSpec>,
        shadow: Option<bool>,
    },
    Plane {
        extend: Option<Vec<String>>,
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
        material: Option<MaterialSpec>,
        shadow: Option<bool>,
    },
    Cube {
        extend: Option<Vec<String>>,
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
        material: Option<MaterialSpec>,
        shadow: Option<bool>,
    },
    Cylinder {
        extend: Option<Vec<String>>,
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
        material: Option<MaterialSpec>,
        shadow: Option<bool>,
        min: Option<f32>,
        max: Option<f32>,
        closed: Option<bool>,
    },
    Cone {
        extend: Option<Vec<String>>,
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
        material: Option<MaterialSpec>,
        shadow: Option<bool>,
        min: Option<f32>,
        max: Option<f32>,
        closed: Option<bool>,
    },
    Triangle {
        extend: Option<Vec<String>>,
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
        material: Option<MaterialSpec>,
        shadow: Option<bool>,
        p1: [f32; 3],
        p2: [f32; 3],
        p3: [f32; 3],
    },
    SmoothTriangle {
        extend: Option<Vec<String>>,
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
        material: Option<MaterialSpec>,
        shadow: Option<bool>,
        p1: [f32; 3],
        p2: [f32; 3],
        p3: [f32; 3],
        n1: [f32; 3],
        n2: [f32; 3],
        n3: [f32; 3],
    },
    Group {
        extend: Option<Vec<String>>,
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
        shapes: Vec<AddShape>,
    },
}

impl AddShape {
    fn make_shape(
        &self,
        define_transforms: &HashMap<String, Vec<TransformSpec>>,
        define_materials: &HashMap<String, Vec<MaterialSpec>>,
    ) -> Shape {
        match self {
            Self::Sphere {
                extend,
                transform,
                material,
                shadow,
            } => {
                let (transform, material) = make_transform_material(
                    transform.as_ref(),
                    material.as_ref(),
                    extend.as_ref(),
                    define_transforms,
                    define_materials,
                );

                Sphere::default()
                    .with_transform(transform)
                    .with_material(material)
                    .with_shadow(shadow.unwrap_or(true))
                    .into()
            }

            Self::Plane {
                extend,
                transform,
                material,
                shadow,
            } => {
                let (transform, material) = make_transform_material(
                    transform.as_ref(),
                    material.as_ref(),
                    extend.as_ref(),
                    define_transforms,
                    define_materials,
                );

                Plane::default()
                    .with_transform(transform)
                    .with_material(material)
                    .with_shadow(shadow.unwrap_or(true))
                    .into()
            }

            Self::Cube {
                extend,
                transform,
                material,
                shadow,
            } => {
                let (transform, material) = make_transform_material(
                    transform.as_ref(),
                    material.as_ref(),
                    extend.as_ref(),
                    define_transforms,
                    define_materials,
                );

                Cube::default()
                    .with_transform(transform)
                    .with_material(material)
                    .with_shadow(shadow.unwrap_or(true))
                    .into()
            }

            Self::Cylinder {
                extend,
                transform,
                material,
                shadow,
                min,
                max,
                closed,
            } => {
                let (transform, material) = make_transform_material(
                    transform.as_ref(),
                    material.as_ref(),
                    extend.as_ref(),
                    define_transforms,
                    define_materials,
                );

                Cylinder::default()
                    .with_transform(transform)
                    .with_material(material)
                    .with_shadow(shadow.unwrap_or(true))
                    .min(min.unwrap_or(f32::NEG_INFINITY))
                    .max(max.unwrap_or(f32::INFINITY))
                    .closed(closed.unwrap_or(false))
                    .into()
            }

            Self::Cone {
                extend,
                transform,
                material,
                shadow,
                min,
                max,
                closed,
            } => {
                let (transform, material) = make_transform_material(
                    transform.as_ref(),
                    material.as_ref(),
                    extend.as_ref(),
                    define_transforms,
                    define_materials,
                );

                Cone::default()
                    .with_transform(transform)
                    .with_material(material)
                    .with_shadow(shadow.unwrap_or(true))
                    .min(min.unwrap_or(f32::NEG_INFINITY))
                    .max(max.unwrap_or(f32::INFINITY))
                    .closed(closed.unwrap_or(false))
                    .into()
            }

            Self::Triangle {
                extend,
                transform,
                material,
                shadow,
                p1,
                p2,
                p3,
            } => {
                let (transform, material) = make_transform_material(
                    transform.as_ref(),
                    material.as_ref(),
                    extend.as_ref(),
                    define_transforms,
                    define_materials,
                );

                Triangle::new((*p1).into(), (*p2).into(), (*p3).into())
                    .with_transform(transform)
                    .with_material(material)
                    .with_shadow(shadow.unwrap_or(true))
                    .into()
            }

            Self::SmoothTriangle {
                extend,
                transform,
                material,
                shadow,
                p1,
                p2,
                p3,
                n1,
                n2,
                n3,
            } => {
                let (transform, material) = make_transform_material(
                    transform.as_ref(),
                    material.as_ref(),
                    extend.as_ref(),
                    define_transforms,
                    define_materials,
                );

                SmoothTriangle::new(
                    (*p1).into(),
                    (*p2).into(),
                    (*p3).into(),
                    (*n1).into(),
                    (*n2).into(),
                    (*n3).into(),
                )
                .with_transform(transform)
                .with_material(material)
                .with_shadow(shadow.unwrap_or(true))
                .into()
            }

            Self::Group { .. } => unreachable!(),
        }
    }

    fn make_group(
        &self,
        define_transforms: &HashMap<String, Vec<TransformSpec>>,
        define_materials: &HashMap<String, Vec<MaterialSpec>>,
    ) -> Group {
        match self {
            Self::Group {
                extend,
                transform,
                shapes,
            } => {
                let (transform, _) = make_transform_material(
                    transform.as_ref(),
                    None,
                    extend.as_ref(),
                    define_transforms,
                    define_materials,
                );
                let mut group = Group::default().with_transform(transform);
                for shape in shapes {
                    match shape {
                        a @ AddShape::Group { .. } => {
                            let g = a.make_group(define_transforms, define_materials);
                            group.add_child(g);
                        }
                        a => {
                            let s = a.make_shape(define_transforms, define_materials);
                            group.add_shape(s);
                        }
                    }
                }
                group
            }
            _ => unreachable!(),
        }
    }
}

fn make_transform_material(
    transform_specs: Option<&Vec<TransformSpec>>,
    material_spec: Option<&MaterialSpec>,
    extend: Option<&Vec<String>>,
    define_transforms: &HashMap<String, Vec<TransformSpec>>,
    define_materials: &HashMap<String, Vec<MaterialSpec>>,
) -> (Transform, Material) {
    let mut final_transform = Transform::default();
    let mut final_material = Material::default();
    if let Some(definitions) = extend {
        for definition in definitions {
            if let Some(specs) = define_transforms.get(definition) {
                for t in specs {
                    final_transform = t.update(final_transform);
                }
            }
            if let Some(specs) = define_materials.get(definition) {
                for m in specs {
                    final_material = m.update(final_material);
                }
            }
        }
    }
    if let Some(specs) = transform_specs {
        for t in specs {
            final_transform = t.update(final_transform);
        }
    }
    if let Some(spec) = material_spec {
        final_material = spec.update(final_material);
    }
    (final_transform, final_material)
}

#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "add")]
#[serde(rename_all = "kebab-case")]
pub enum AddCsg {
    Csg { args: [CsgElem; 2], op: CsgOp },
}

impl AddCsg {
    pub fn make_csg(
        &self,
        define_transforms: &HashMap<String, Vec<TransformSpec>>,
        define_materials: &HashMap<String, Vec<MaterialSpec>>,
    ) -> Csg {
        match self {
            Self::Csg { args, op } => {
                let left = args[0].make_csg_child(define_transforms, define_materials);
                let right = args[1].make_csg_child(define_transforms, define_materials);
                Csg::new(*op, left, right)
            }
        }
    }
}

#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "kind")]
#[serde(rename_all = "kebab-case")]
pub enum CsgElem {
    Csg {
        op: CsgOp,
        args: [Box<CsgElem>; 2],
    },
    Sphere {
        extend: Option<Vec<String>>,
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
        material: Option<MaterialSpec>,
        shadow: Option<bool>,
    },
    Plane {
        extend: Option<Vec<String>>,
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
        material: Option<MaterialSpec>,
        shadow: Option<bool>,
    },
    Cube {
        extend: Option<Vec<String>>,
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
        material: Option<MaterialSpec>,
        shadow: Option<bool>,
    },
    Cylinder {
        extend: Option<Vec<String>>,
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
        material: Option<MaterialSpec>,
        shadow: Option<bool>,
        min: Option<f32>,
        max: Option<f32>,
        closed: Option<bool>,
    },
    Cone {
        extend: Option<Vec<String>>,
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
        material: Option<MaterialSpec>,
        shadow: Option<bool>,
        min: Option<f32>,
        max: Option<f32>,
        closed: Option<bool>,
    },
    Triangle {
        extend: Option<Vec<String>>,
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
        material: Option<MaterialSpec>,
        shadow: Option<bool>,
        p1: [f32; 3],
        p2: [f32; 3],
        p3: [f32; 3],
    },
    SmoothTriangle {
        extend: Option<Vec<String>>,
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
        material: Option<MaterialSpec>,
        shadow: Option<bool>,
        p1: [f32; 3],
        p2: [f32; 3],
        p3: [f32; 3],
        n1: [f32; 3],
        n2: [f32; 3],
        n3: [f32; 3],
    },
}

impl CsgElem {
    pub fn make_csg_child(
        &self,
        define_transforms: &HashMap<String, Vec<TransformSpec>>,
        define_materials: &HashMap<String, Vec<MaterialSpec>>,
    ) -> CsgChild {
        match &self {
            Self::Sphere { .. }
            | Self::Plane { .. }
            | Self::Cube { .. }
            | Self::Cylinder { .. }
            | Self::Cone { .. }
            | Self::Triangle { .. }
            | Self::SmoothTriangle { .. } => {
                let spec = AddShape::from(self.clone());
                let shape = spec.make_shape(define_transforms, define_materials);
                CsgChild::Shape(Box::new(shape))
            }
            Self::Csg { op, args } => {
                let left = args[0].make_csg_child(define_transforms, define_materials);
                let right = args[1].make_csg_child(define_transforms, define_materials);
                let csg = Csg::new(*op, left, right);
                CsgChild::Csg(Box::new(csg))
            }
        }
    }
}

impl From<CsgElem> for AddShape {
    fn from(csg: CsgElem) -> Self {
        match csg {
            CsgElem::Sphere {
                extend,
                transform,
                material,
                shadow,
            } => Self::Sphere {
                extend,
                transform,
                material,
                shadow,
            },

            CsgElem::Plane {
                extend,
                transform,
                material,
                shadow,
            } => Self::Plane {
                extend,
                transform,
                material,
                shadow,
            },

            CsgElem::Cube {
                extend,
                transform,
                material,
                shadow,
            } => Self::Cube {
                extend,
                transform,
                material,
                shadow,
            },

            CsgElem::Cylinder {
                extend,
                transform,
                material,
                shadow,
                min,
                max,
                closed,
            } => Self::Cylinder {
                extend,
                transform,
                material,
                shadow,
                min,
                max,
                closed,
            },

            CsgElem::Cone {
                extend,
                transform,
                material,
                shadow,
                min,
                max,
                closed,
            } => Self::Cone {
                extend,
                transform,
                material,
                shadow,
                min,
                max,
                closed,
            },

            CsgElem::Triangle {
                extend,
                transform,
                material,
                shadow,
                p1,
                p2,
                p3,
            } => Self::Triangle {
                extend,
                transform,
                material,
                shadow,
                p1,
                p2,
                p3,
            },

            CsgElem::SmoothTriangle {
                extend,
                transform,
                material,
                shadow,
                p1,
                p2,
                p3,
                n1,
                n2,
                n3,
            } => Self::SmoothTriangle {
                extend,
                transform,
                material,
                shadow,
                p1,
                p2,
                p3,
                n1,
                n2,
                n3,
            },

            CsgElem::Csg { .. } => unreachable!(),
        }
    }
}

#[derive(Debug, Clone, Deserialize)]
pub struct Define {
    pub define: String,
    #[serde(default)]
    pub extend: Option<Vec<String>>,
    #[serde(default, deserialize_with = "deser_transform")]
    pub transform: Option<Vec<TransformSpec>>,
    #[serde(default)]
    pub material: Option<MaterialSpec>,
}

#[derive(Debug, Clone)]
pub enum TransformSpec {
    Translate {
        x: f32,
        y: f32,
        z: f32,
    },
    Scale {
        x: f32,
        y: f32,
        z: f32,
    },
    RotateX {
        angle: f32,
    },
    RotateY {
        angle: f32,
    },
    RotateZ {
        angle: f32,
    },
    Shear {
        x_y: f32,
        x_z: f32,
        y_x: f32,
        y_z: f32,
        z_x: f32,
        z_y: f32,
    },
}

impl TransformSpec {
    pub fn update(&self, transform: Transform) -> Transform {
        match *self {
            Self::Translate { x, y, z } => transform.translation(x, y, z),
            Self::Scale { x, y, z } => transform.scaling(x, y, z),
            Self::RotateX { angle } => transform.rotation_x(angle),
            Self::RotateY { angle } => transform.rotation_y(angle),
            Self::RotateZ { angle } => transform.rotation_z(angle),
            Self::Shear {
                x_y,
                x_z,
                y_x,
                y_z,
                z_x,
                z_y,
            } => transform.shearing(x_y, x_z, y_x, y_z, z_x, z_y),
        }
    }
}

fn deser_transform<'de, D>(deserializer: D) -> Result<Option<Vec<TransformSpec>>, D::Error>
where
    D: Deserializer<'de>,
{
    use serde_yaml::Value;

    #[derive(Debug, Deserialize)]
    #[serde(rename_all = "kebab-case")]
    enum Op {
        Translate,
        Scale,
        RotateX,
        RotateY,
        RotateZ,
        Shear,
    }

    let exprs = Option::<Vec<Vec<Value>>>::deserialize(deserializer)?;
    if exprs.is_none() || (exprs.is_some() && exprs.as_ref().unwrap().is_empty()) {
        return Ok(None);
    }
    let exprs = exprs.unwrap();

    let mut tranforms = vec![];
    for expr in exprs {
        let op = expr
            .get(0)
            .ok_or_else(|| de::Error::custom("Missing transform operator"))
            .and_then(|op| match op {
                Value::String(op) => serde_plain::from_str::<Op>(op),
                _ => Err(de::Error::custom(format!(
                    "Invalid transform operator: {op:?}"
                ))),
            })
            .map_err(de::Error::custom)?;
        let operands = &expr[1..];
        match op {
            Op::Translate => match operands {
                [Value::Number(x), Value::Number(y), Value::Number(z)] => {
                    tranforms.push(TransformSpec::Translate {
                        x: to_f32(x),
                        y: to_f32(y),
                        z: to_f32(z),
                    });
                }
                _ => {
                    return Err(de::Error::custom(format!(
                        "Invalid transform translate operands: {operands:?}",
                    )))
                }
            },

            Op::Scale => match operands {
                [Value::Number(x), Value::Number(y), Value::Number(z)] => {
                    tranforms.push(TransformSpec::Scale {
                        x: to_f32(x),
                        y: to_f32(y),
                        z: to_f32(z),
                    });
                }
                _ => {
                    return Err(de::Error::custom(format!(
                        "Invalid transform scale operands: {operands:?}",
                    )))
                }
            },

            Op::RotateX => match operands {
                [Value::String(angle)] => match meval::eval_str(angle).map(|angle| angle as f32) {
                    Ok(angle) => {
                        tranforms.push(TransformSpec::RotateX { angle });
                    }
                    _ => {
                        return Err(de::Error::custom(format!(
                            "Invalid transform rotate-x operand value: {operands:?}",
                        )))
                    }
                },
                [Value::Number(angle)] => {
                    tranforms.push(TransformSpec::RotateX {
                        angle: to_f32(angle),
                    });
                }
                _ => {
                    return Err(de::Error::custom(format!(
                        "Invalid transform rotate-x operand: {operands:?}",
                    )))
                }
            },

            Op::RotateY => match operands {
                [Value::String(angle)] => match meval::eval_str(angle).map(|angle| angle as f32) {
                    Ok(angle) => {
                        tranforms.push(TransformSpec::RotateY { angle });
                    }
                    _ => {
                        return Err(de::Error::custom(format!(
                            "Invalid transform rotate-y operand value: {operands:?}",
                        )))
                    }
                },
                [Value::Number(angle)] => {
                    tranforms.push(TransformSpec::RotateY {
                        angle: to_f32(angle),
                    });
                }
                _ => {
                    return Err(de::Error::custom(format!(
                        "Invalid transform rotate-y operand: {operands:?}",
                    )))
                }
            },

            Op::RotateZ => match operands {
                [Value::String(angle)] => match meval::eval_str(angle).map(|angle| angle as f32) {
                    Ok(angle) => {
                        tranforms.push(TransformSpec::RotateZ { angle });
                    }
                    _ => {
                        return Err(de::Error::custom(format!(
                            "Invalid transform rotate-z operand value: {operands:?}",
                        )))
                    }
                },
                [Value::Number(angle)] => {
                    tranforms.push(TransformSpec::RotateZ {
                        angle: to_f32(angle),
                    });
                }

                _ => {
                    return Err(de::Error::custom(format!(
                        "Invalid transform rotate-z operand: {operands:?}",
                    )))
                }
            },

            Op::Shear => match operands {
                [Value::Number(x_y), Value::Number(x_z), Value::Number(y_x), Value::Number(y_z), Value::Number(z_x), Value::Number(z_y)] =>
                {
                    tranforms.push(TransformSpec::Shear {
                        x_y: to_f32(x_y),
                        x_z: to_f32(x_z),
                        y_x: to_f32(y_x),
                        y_z: to_f32(y_z),
                        z_x: to_f32(z_x),
                        z_y: to_f32(z_y),
                    });
                }
                _ => {
                    return Err(de::Error::custom(format!(
                        "Invalid transform shear operands: {operands:?}",
                    )))
                }
            },
        }
    }

    Ok(Some(tranforms))
}

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub struct MaterialSpec {
    pub pattern: Option<PatternSpec>,
    pub color: Option<[f32; 3]>,
    pub ambient: Option<f32>,
    pub diffuse: Option<f32>,
    pub specular: Option<f32>,
    pub shininess: Option<f32>,
    pub reflective: Option<f32>,
    pub transparency: Option<f32>,
    pub refractive_index: Option<f32>,
}

impl MaterialSpec {
    pub fn update(&self, mut material: Material) -> Material {
        if let Some(pattern) = &self.pattern {
            material = material.pattern(pattern)
        }
        if let Some(color) = self.color {
            material = material.color(color.into());
        }
        if let Some(ambient) = self.ambient {
            material = material.ambient(ambient);
        }
        if let Some(diffuse) = self.diffuse {
            material = material.diffuse(diffuse);
        }
        if let Some(specular) = self.specular {
            material = material.specular(specular);
        }
        if let Some(shininess) = self.shininess {
            material = material.shininess(shininess);
        }
        if let Some(reflective) = self.reflective {
            material = material.reflective(reflective);
        }
        if let Some(transparency) = self.transparency {
            material = material.transparency(transparency);
        }
        if let Some(refractive_index) = self.refractive_index {
            material = material.refractive_index(refractive_index);
        }
        material
    }
}

#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "kind")]
#[serde(rename_all = "kebab-case")]
pub enum PatternSpec {
    Striped {
        colors: [[f32; 3]; 2],
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
    },
    Gradient {
        colors: [[f32; 3]; 2],
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
    },
    Ring {
        colors: [[f32; 3]; 2],
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
    },
    Checker {
        colors: [[f32; 3]; 2],
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
    },
    XyzRgb {
        #[serde(default, deserialize_with = "deser_transform")]
        transform: Option<Vec<TransformSpec>>,
    },
}

impl From<&PatternSpec> for Pattern {
    fn from(source: &PatternSpec) -> Self {
        let mut base_transform = Transform::default();
        match source {
            PatternSpec::Striped { transform, .. }
            | PatternSpec::Gradient { transform, .. }
            | PatternSpec::Ring { transform, .. }
            | PatternSpec::Checker { transform, .. }
            | PatternSpec::XyzRgb { transform, .. } => {
                if let Some(specs) = transform {
                    for op in specs {
                        base_transform = op.update(base_transform);
                    }
                }
            }
        }
        match source {
            PatternSpec::Striped { colors, .. } => Striped::new(colors[0].into(), colors[1].into())
                .with_transform(base_transform)
                .into(),
            PatternSpec::Gradient { colors, .. } => {
                Gradient::new(colors[0].into(), colors[1].into())
                    .with_transform(base_transform)
                    .into()
            }
            PatternSpec::Ring { colors, .. } => Ring::new(colors[0].into(), colors[1].into())
                .with_transform(base_transform)
                .into(),
            PatternSpec::Checker { colors, .. } => Checker::new(colors[0].into(), colors[1].into())
                .with_transform(base_transform)
                .into(),
            PatternSpec::XyzRgb { .. } => XyzRgb::new().with_transform(base_transform).into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deser_scene() {
        let scene = "
- add: camera
  width: 1280
  height: 720
  field-of-view: pi/3
  from: [0.0, 1.5, 5.0]
  to: [0, 0, 0]
  up: [0, 1, 0]
- add: point-light
  at: [-10, 10, -10]
  intensity: [1, 1, 1]
- define: my-def
  extend: [some-other-def]
  transform:
    - [rotate-y, pi/4]
    - [scale, 0.5, 0.5, 0.5]
    - [translate, 0, 1, 0]
  material:
    pattern:
      kind: checker
      colors:
        - [0, 0, 0]
        - [1, 1, 1]
    color: [0.1, 0.2, 0.3]
- add: cube
  extend: []
  transform: null
  material: null
  shadow: false
- add: group
  # extend: [my-def]
  transform: null
  shapes:
    - add: cube
      extend: [my-def]
";

        let scene = serde_yaml::from_str::<Scene>(&scene);
        assert!(scene.is_ok());
    }
}
