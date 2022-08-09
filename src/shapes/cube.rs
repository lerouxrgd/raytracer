use derivative::Derivative;
use ordered_float::OrderedFloat;
use slotmap::Key;

use crate::groups::Group;
use crate::intersections::Intersection;
use crate::materials::Material;
use crate::matrices::Matrix;
use crate::rays::Ray;
use crate::shapes::Shape;
use crate::tuples::{Point, Vector};

/// An axis-aligned bounding box (AABB) centered on the origin
#[derive(Debug, Clone, Copy, Derivative)]
#[derivative(PartialEq)]
pub struct Cube {
    pub(super) transform: Matrix<4, 4>,
    pub(super) material: Material,
    #[derivative(PartialEq = "ignore")]
    pub(super) parent: Group,
}

impl Cube {
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

    fn check_axis(origin: f32, direction: f32) -> (f32, f32) {
        let tmin_numerator = -1. - origin;
        let tmax_numerator = 1. - origin;

        let (tmin, tmax) = if direction.abs() >= Self::EPSILON {
            (tmin_numerator / direction, tmax_numerator / direction)
        } else {
            (
                tmin_numerator * f32::INFINITY,
                tmax_numerator * f32::INFINITY,
            )
        };

        if tmin > tmax {
            (tmax, tmin)
        } else {
            (tmin, tmax)
        }
    }

    pub fn local_intersect(&self, local_ray: Ray) -> Option<[Intersection; 2]> {
        let (xtmin, xtmax) = Self::check_axis(local_ray.origin.x(), local_ray.direction.x());
        let (ytmin, ytmax) = Self::check_axis(local_ray.origin.y(), local_ray.direction.y());
        let (ztmin, ztmax) = Self::check_axis(local_ray.origin.z(), local_ray.direction.z());

        let tmin = [xtmin, ytmin, ztmin]
            .into_iter()
            .map(OrderedFloat::from)
            .max()
            .unwrap();
        let tmax = [xtmax, ytmax, ztmax]
            .into_iter()
            .map(OrderedFloat::from)
            .min()
            .unwrap();

        if tmin > tmax {
            None
        } else {
            Some([
                Intersection::new(tmin.into(), self.into()),
                Intersection::new(tmax.into(), self.into()),
            ])
        }
    }

    pub fn local_normal_at(&self, local_point: Point) -> Vector {
        let maxc: f32 = [local_point.x(), local_point.y(), local_point.z()]
            .into_iter()
            .map(|v| OrderedFloat(v.abs()))
            .max()
            .unwrap()
            .into();

        if local_point.x().abs() == maxc {
            Vector::new(local_point.x(), 0., 0.)
        } else if local_point.y().abs() == maxc {
            Vector::new(0., local_point.y(), 0.)
        } else if local_point.z().abs() == maxc {
            Vector::new(0., 0., local_point.z())
        } else {
            unreachable!()
        }
    }
}

impl Default for Cube {
    fn default() -> Self {
        Self {
            transform: Matrix::identity(),
            material: Material::default(),
            parent: Group::null(),
        }
    }
}

impl From<Cube> for Shape {
    fn from(c: Cube) -> Self {
        Self::Cube(c)
    }
}

impl From<&Cube> for Shape {
    fn from(c: &Cube) -> Self {
        Self::Cube(*c)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tuples::Vector;

    #[test]
    fn cube_basics() {
        let cases = vec![
            (Point::new(5., 0.5, 0.), Vector::new(-1., 0., 0.), 4., 6.), // +x
            (Point::new(-5., 0.5, 0.), Vector::new(1., 0., 0.), 4., 6.), // -x
            (Point::new(0.5, 5., 0.), Vector::new(0., -1., 0.), 4., 6.), // +y
            (Point::new(0.5, -5., 0.), Vector::new(0., 1., 0.), 4., 6.), // -y
            (Point::new(0.5, 0., 5.), Vector::new(0., 0., -1.), 4., 6.), // +z
            (Point::new(0.5, 0., -5.), Vector::new(0., 0., 1.), 4., 6.), // -z
            (Point::new(0., 0.5, 0.), Vector::new(0., 0., 1.), -1., 1.), // inside
        ];
        for (origin, direction, t1, t2) in cases.into_iter() {
            let c = Cube::default();
            let r = Ray::new(origin, direction);
            let xs = c.local_intersect(r).unwrap();
            assert!(xs[0].t() == t1);
            assert!(xs[1].t() == t2);
        }

        let cases = vec![
            (Point::new(-2., 0., 0.), Vector::new(0.2673, 0.5345, 0.8018)),
            (Point::new(0., -2., 0.), Vector::new(0.8018, 0.2673, 0.5345)),
            (Point::new(0., 0., -2.), Vector::new(0.5345, 0.8018, 0.2673)),
            (Point::new(2., 0., 2.), Vector::new(0., 0., -1.)),
            (Point::new(0., 2., 2.), Vector::new(0., -1., 0.)),
            (Point::new(2., 2., 0.), Vector::new(-1., 0., 0.)),
        ];
        for (origin, direction) in cases.into_iter() {
            let c = Cube::default();
            let r = Ray::new(origin, direction);
            let xs = c.local_intersect(r);
            assert!(xs.is_none());
        }

        let cases = vec![
            (Point::new(1., 0.5, -0.8), Vector::new(1., 0., 0.)),
            (Point::new(-1., -0.2, -0.9), Vector::new(-1., 0., 0.)),
            (Point::new(-0.4, 1., -0.1), Vector::new(0., 1., 0.)),
            (Point::new(0.3, -1., 0.7), Vector::new(0., -1., 0.)),
            (Point::new(-0.6, 0.3, 1.), Vector::new(0., 0., 1.)),
            (Point::new(0.4, 0.4, -1.), Vector::new(0., 0., -1.)),
            (Point::new(1., 1., 1.), Vector::new(1., 0., 0.)),
            (Point::new(-1., -1., -1.), Vector::new(-1., 0., 0.)),
        ];
        for (point, normal) in cases.into_iter() {
            let c = Cube::default();
            assert!(c.local_normal_at(point) == normal);
        }
    }
}
