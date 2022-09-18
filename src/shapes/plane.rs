use derivative::Derivative;
use slotmap::Key;

use crate::bounds::BoundingBox;
use crate::groups::Group;
use crate::intersections::Intersection;
use crate::materials::Material;
use crate::matrices::Matrix;
use crate::rays::Ray;
use crate::shapes::Shape;
use crate::tuples::{Point, Vector};

/// A xz plane (left-handed coordinates)
#[derive(Debug, Clone, Copy, Derivative)]
#[derivative(PartialEq)]
pub struct Plane {
    pub(super) transform: Matrix<4, 4>,
    pub(super) material: Material,
    pub(super) shadow: bool,
    #[derivative(PartialEq = "ignore")]
    pub(super) parent: Group,
}

impl Plane {
    pub const EPSILON: f32 = 1e-4;

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

    pub fn local_intersect(&self, local_ray: Ray) -> Option<Intersection> {
        if local_ray.direction.y().abs() < Self::EPSILON {
            None
        } else {
            let t = -local_ray.origin.y() / local_ray.direction.y();
            Some(Intersection::new(t, self.into()))
        }
    }

    pub fn local_normal_at(&self, _local_point: Point) -> Vector {
        Vector::new(0., 1., 0.)
    }

    pub fn bounds(&self) -> BoundingBox {
        BoundingBox::default()
            .with_min(Point::new(f32::NEG_INFINITY, 0., f32::NEG_INFINITY))
            .with_max(Point::new(f32::INFINITY, 0., f32::INFINITY))
    }
}

impl Default for Plane {
    fn default() -> Self {
        Self {
            transform: Matrix::identity(),
            material: Material::default(),
            shadow: true,
            parent: Group::null(),
        }
    }
}

impl From<Plane> for Shape {
    fn from(p: Plane) -> Self {
        Self::Plane(p)
    }
}

impl From<&Plane> for Shape {
    fn from(p: &Plane) -> Self {
        Self::Plane(*p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tuples::Vector;

    #[test]
    fn plane_basics() {
        let p = Plane::default();
        assert!(p.local_normal_at(Point::new(0., 0., 0.)) == Vector::new(0., 1., 0.));
        assert!(p.local_normal_at(Point::new(10., 0., -10.)) == Vector::new(0., 1., 0.));
        assert!(p.local_normal_at(Point::new(-5., 0., 150.)) == Vector::new(0., 1., 0.));

        let r = Ray::new(Point::new(0., 10., 0.), Vector::new(0., 0., 1.));
        assert!(p.local_intersect(r).is_none());

        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        assert!(p.local_intersect(r).is_none());

        let r = Ray::new(Point::new(0., 1., 0.), Vector::new(0., -1., 0.));
        assert!(p.local_intersect(r).unwrap() == Intersection::new(1., p.into()));

        let r = Ray::new(Point::new(0., -1., 0.), Vector::new(0., 1., 0.));
        assert!(p.local_intersect(r).unwrap() == Intersection::new(1., p.into()));
    }
}
