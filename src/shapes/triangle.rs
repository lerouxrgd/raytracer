use derivative::Derivative;
use slotmap::Key;

use crate::groups::Group;
use crate::intersections::Intersection;
use crate::materials::Material;
use crate::matrices::Matrix;
use crate::rays::Ray;
use crate::shapes::Shape;
use crate::tuples::{Point, Vector};

#[derive(Debug, Clone, Copy, Derivative)]
#[derivative(PartialEq)]
pub struct Triangle {
    pub(super) transform: Matrix<4, 4>,
    pub(super) material: Material,
    pub(super) shadow: bool,
    pub(super) p1: Point,
    pub(super) p2: Point,
    pub(super) p3: Point,
    pub(super) e1: Vector,
    pub(super) e2: Vector,
    pub(super) normal: Vector,
    #[derivative(PartialEq = "ignore")]
    pub(super) parent: Group,
}

impl Triangle {
    pub const EPSILON: f32 = 1e-4;

    pub fn new(p1: Point, p2: Point, p3: Point) -> Self {
        let e1 = p2 - p1;
        let e2 = p3 - p1;
        let normal = e2.cross(e1).normalize();
        Self {
            transform: Matrix::identity(),
            material: Material::default(),
            shadow: true,
            p1,
            p2,
            p3,
            e1,
            e2,
            normal,
            parent: Group::null(),
        }
    }

    pub fn with_transform<T: Into<Matrix<4, 4>>>(mut self, transform: T) -> Self {
        self.transform = transform.into();
        self
    }

    pub fn get_transform(&self) -> Matrix<4, 4> {
        self.transform
    }

    pub fn set_transform<T: Into<Matrix<4, 4>>>(&mut self, t: T) {
        self.transform = t.into()
    }

    pub fn with_material(mut self, material: Material) -> Self {
        self.material = material;
        self
    }

    pub fn get_material(&self) -> Material {
        self.material
    }

    pub fn set_material(&mut self, m: Material) {
        self.material = m
    }

    pub fn with_shadow(mut self, shadow: bool) -> Self {
        self.shadow = shadow;
        self
    }

    pub fn get_shadow(&self) -> bool {
        self.shadow
    }

    pub fn set_shadow(&mut self, s: bool) {
        self.shadow = s
    }

    pub fn local_normal_at(&self, _local_point: Point) -> Vector {
        self.normal
    }

    pub fn local_intersect(&self, local_ray: Ray) -> Option<Intersection> {
        let dir_cross_e2 = local_ray.direction.cross(self.e2);
        let det = self.e1.dot(dir_cross_e2);
        if det.abs() < Self::EPSILON {
            return None;
        }

        let f = 1. / det;
        let p1_to_origin = local_ray.origin - self.p1;
        let u = f * p1_to_origin.dot(dir_cross_e2);
        if !(0. ..=1.).contains(&u) {
            return None;
        }

        let origin_cross_e1 = p1_to_origin.cross(self.e1);
        let v = f * local_ray.direction.dot(origin_cross_e1);
        if v < 0. || (u + v) > 1. {
            return None;
        }

        let t = f * self.e2.dot(origin_cross_e1);
        Some(Intersection::new(t, self.into()))
    }
}

impl From<Triangle> for Shape {
    fn from(t: Triangle) -> Self {
        Self::Triangle(t)
    }
}

impl From<&Triangle> for Shape {
    fn from(t: &Triangle) -> Self {
        Self::Triangle(*t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tuples::Vector;

    #[test]
    fn triangle_basics() {
        let t = Triangle::new(
            Point::new(0., 1., 0.),
            Point::new(-1., 0., 0.),
            Point::new(1., 0., 0.),
        );
        assert!(t.e1 == Vector::new(-1., -1., 0.));
        assert!(t.e2 == Vector::new(1., -1., 0.));
        assert!(t.normal == Vector::new(0., 0., -1.));

        let r = Ray::new(Point::new(0., -1., -2.), Vector::new(0., 1., 0.));
        assert!(t.local_intersect(r).is_none());

        let r = Ray::new(Point::new(1., 1., -2.), Vector::new(0., 0., 1.));
        assert!(t.local_intersect(r).is_none());
        let r = Ray::new(Point::new(-1., 1., -2.), Vector::new(0., 0., 1.));
        assert!(t.local_intersect(r).is_none());
        let r = Ray::new(Point::new(0., -1., -2.), Vector::new(0., 0., 1.));
        assert!(t.local_intersect(r).is_none());

        let r = Ray::new(Point::new(0., 0.5, -2.), Vector::new(0., 0., 1.));
        assert!(t.local_intersect(r).unwrap().t() == 2.);
    }
}
