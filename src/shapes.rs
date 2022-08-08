mod cone;
mod cube;
mod cylinder;
mod plane;
mod sphere;
mod triangle;

use slotmap::Key;

use crate::groups::Group;
use crate::intersections::Intersection;
use crate::materials::Material;
use crate::matrices::Matrix;
use crate::rays::Ray;
use crate::tuples::{Point, Tuple, Vector};

pub use cone::Cone;
pub use cube::Cube;
pub use cylinder::Cylinder;
pub use plane::Plane;
pub use sphere::Sphere;
pub use triangle::Triangle;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Shape {
    Sphere(Sphere),
    Plane(Plane),
    Cube(Cube),
    Cylinder(Cylinder),
    Cone(Cone),
    Triangle(Triangle),
}

impl Shape {
    pub fn sphere() -> Self {
        Sphere::default().into()
    }

    pub fn plane() -> Self {
        Plane::default().into()
    }

    pub fn cube() -> Self {
        Cube::default().into()
    }

    pub fn cylinder() -> Self {
        Cylinder::default().into()
    }

    pub fn cone() -> Self {
        Cone::default().into()
    }

    pub fn triangle(p1: Point, p2: Point, p3: Point) -> Self {
        Triangle::new(p1, p2, p3).into()
    }

    pub fn with_transform<T: Into<Matrix<4, 4>>>(mut self, transform: T) -> Self {
        self.set_transform(transform.into());
        self
    }

    pub fn get_transform(&self) -> Matrix<4, 4> {
        match self {
            &Self::Sphere(Sphere { transform, .. })
            | &Self::Plane(Plane { transform, .. })
            | &Self::Cylinder(Cylinder { transform, .. })
            | &Self::Cone(Cone { transform, .. })
            | &Self::Triangle(Triangle { transform, .. })
            | &Self::Cube(Cube { transform, .. }) => transform,
        }
    }

    pub fn set_transform<T: Into<Matrix<4, 4>>>(&mut self, t: T) {
        match self {
            Self::Sphere(Sphere {
                ref mut transform, ..
            })
            | Self::Plane(Plane {
                ref mut transform, ..
            })
            | Self::Cylinder(Cylinder {
                ref mut transform, ..
            })
            | Self::Cone(Cone {
                ref mut transform, ..
            })
            | Self::Triangle(Triangle {
                ref mut transform, ..
            })
            | Self::Cube(Cube {
                ref mut transform, ..
            }) => *transform = t.into(),
        }
    }

    pub fn with_material(mut self, material: Material) -> Self {
        self.set_material(material);
        self
    }

    pub fn get_material(&self) -> Material {
        match self {
            &Self::Sphere(Sphere { material, .. })
            | &Self::Plane(Plane { material, .. })
            | &Self::Cylinder(Cylinder { material, .. })
            | &Self::Cone(Cone { material, .. })
            | &Self::Triangle(Triangle { material, .. })
            | &Self::Cube(Cube { material, .. }) => material,
        }
    }

    pub fn set_material(&mut self, m: Material) {
        match self {
            Self::Sphere(Sphere {
                ref mut material, ..
            })
            | Self::Plane(Plane {
                ref mut material, ..
            })
            | Self::Cylinder(Cylinder {
                ref mut material, ..
            })
            | Self::Cone(Cone {
                ref mut material, ..
            })
            | Self::Triangle(Triangle {
                ref mut material, ..
            })
            | Self::Cube(Cube {
                ref mut material, ..
            }) => *material = m,
        }
    }

    #[cfg(test)]
    pub fn material_mut(&mut self) -> &mut Material {
        match self {
            Self::Sphere(Sphere {
                ref mut material, ..
            })
            | Self::Plane(Plane {
                ref mut material, ..
            })
            | Self::Cylinder(Cylinder {
                ref mut material, ..
            })
            | Self::Cone(Cone {
                ref mut material, ..
            })
            | Self::Triangle(Triangle {
                ref mut material, ..
            })
            | Self::Cube(Cube {
                ref mut material, ..
            }) => material,
        }
    }

    pub fn normal_at(&self, world_point: Point) -> Vector {
        let local_point = self.world_to_object(world_point);
        let local_normal = match self {
            Self::Sphere(s) => s.local_normal_at(local_point),
            Self::Plane(p) => p.local_normal_at(local_point),
            Self::Cylinder(c) => c.local_normal_at(local_point),
            Self::Cone(c) => c.local_normal_at(local_point),
            Self::Triangle(t) => t.local_normal_at(local_point),
            Self::Cube(c) => c.local_normal_at(local_point),
        };
        self.normal_to_world(local_normal)
    }

    pub fn intersect(&self, intersections: &mut Vec<Intersection>, world_ray: Ray) {
        let inverse = self.get_transform().inverse().unwrap();
        let local_ray = world_ray.transform(inverse);
        match self {
            Shape::Sphere(s) => {
                if let Some(xs) = s.local_intersect(local_ray) {
                    intersections.extend(xs)
                }
            }
            Shape::Plane(p) => {
                if let Some(xs) = p.local_intersect(local_ray) {
                    intersections.push(xs)
                }
            }
            Shape::Cylinder(c) => {
                c.local_intersect(local_ray)
                    .into_iter()
                    .flatten()
                    .for_each(|xs| intersections.push(xs));
            }
            Shape::Cone(c) => {
                c.local_intersect(local_ray)
                    .into_iter()
                    .flatten()
                    .for_each(|xs| intersections.push(xs));
            }
            Shape::Triangle(t) => {
                if let Some(xs) = t.local_intersect(local_ray) {
                    intersections.push(xs)
                }
            }
            Shape::Cube(s) => {
                if let Some(xs) = s.local_intersect(local_ray) {
                    intersections.extend(xs)
                }
            }
        }
    }

    pub fn parent(&self) -> Group {
        match self {
            &Self::Sphere(Sphere { parent, .. })
            | &Self::Plane(Plane { parent, .. })
            | &Self::Cylinder(Cylinder { parent, .. })
            | &Self::Cone(Cone { parent, .. })
            | &Self::Triangle(Triangle { parent, .. })
            | &Self::Cube(Cube { parent, .. }) => parent,
        }
    }

    pub(crate) fn parent_mut(&mut self) -> &mut Group {
        match self {
            Self::Sphere(Sphere { ref mut parent, .. })
            | Self::Plane(Plane { ref mut parent, .. })
            | Self::Cylinder(Cylinder { ref mut parent, .. })
            | Self::Cone(Cone { ref mut parent, .. })
            | Self::Triangle(Triangle { ref mut parent, .. })
            | Self::Cube(Cube { ref mut parent, .. }) => parent,
        }
    }

    pub fn world_to_object(&self, point: Point) -> Point {
        let point = if !self.parent().is_null() {
            self.parent().world_to_object(point)
        } else {
            point
        };
        self.get_transform().inverse().unwrap() * point
    }

    pub fn normal_to_world(&self, normal: Vector) -> Vector {
        let inverse = self.get_transform().inverse().unwrap();
        let normal = Tuple::from(inverse.transpose() * normal);
        let normal = Vector::new(normal[0], normal[1], normal[2]);
        let normal = normal.normalize();

        if !self.parent().is_null() {
            self.parent().normal_to_world(normal)
        } else {
            normal
        }
    }
}
