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
pub struct SmoothTriangle {
    pub(super) transform: Matrix<4, 4>,
    pub(super) material: Material,
    pub(super) p1: Point,
    pub(super) p2: Point,
    pub(super) p3: Point,
    pub(super) e1: Vector,
    pub(super) e2: Vector,
    pub(super) n1: Vector,
    pub(super) n2: Vector,
    pub(super) n3: Vector,
    #[derivative(PartialEq = "ignore")]
    pub(super) parent: Group,
}

impl SmoothTriangle {
    pub const EPSILON: f32 = 1e-4;

    pub fn new(p1: Point, p2: Point, p3: Point, n1: Vector, n2: Vector, n3: Vector) -> Self {
        let e1 = p2 - p1;
        let e2 = p3 - p1;
        Self {
            transform: Matrix::identity(),
            material: Material::default(),
            p1,
            p2,
            p3,
            e1,
            e2,
            n1,
            n2,
            n3,
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

    pub fn local_normal_at(&self, _local_point: Point, hit: Intersection) -> Vector {
        hit.u() * self.n2 + hit.v() * self.n3 + (1. - hit.u() - hit.v()) * self.n1
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
        Some(Intersection::new(t, self.into()).with_uv(u, v))
    }
}

impl From<SmoothTriangle> for Shape {
    fn from(st: SmoothTriangle) -> Self {
        Self::SmoothTriangle(st)
    }
}

impl From<&SmoothTriangle> for Shape {
    fn from(st: &SmoothTriangle) -> Self {
        Self::SmoothTriangle(*st)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::intersections::Computations;
    use crate::tuples::Vector;

    #[test]
    fn smooth_triangle_basics() {
        let tri = SmoothTriangle::new(
            Point::new(0., 1., 0.),
            Point::new(-1., 0., 0.),
            Point::new(1., 0., 0.),
            Vector::new(0., 1., 0.),
            Vector::new(-1., 0., 0.),
            Vector::new(1., 0., 0.),
        );

        let r = Ray::new(Point::new(-0.2, 0.3, -2.), Vector::new(0., 0., 1.));
        let xs = tri.local_intersect(r).unwrap();
        assert!(xs.u() == 0.45);
        assert!(xs.v() == 0.25);

        let xs = Intersection::new(1., tri.into()).with_uv(0.45, 0.25);
        let n = Shape::from(tri).normal_at(Point::new(0., 0., 0.), xs);
        assert!(n.equal_approx(Vector::new(-0.5547, 0.83205, 0.)));

        let xs = Intersection::new(1., tri.into()).with_uv(0.45, 0.25);
        let r = Ray::new(Point::new(-0.2, 0.3, -2.), Vector::new(0., 0., 1.));
        let comps = Computations::prepare(xs, r, &vec![xs].into());
        assert!(comps
            .normalv
            .equal_approx(Vector::new(-0.5547, 0.83205, 0.)));
    }
}
