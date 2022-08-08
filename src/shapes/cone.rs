use std::mem;

use derivative::Derivative;
use slotmap::Key;

use crate::groups::Group;
use crate::intersections::Intersection;
use crate::materials::Material;
use crate::matrices::Matrix;
use crate::rays::Ray;
use crate::shapes::Shape;
use crate::tuples::{Point, Vector};

/// A double-napped cone along the y axis
#[derive(Debug, Clone, Copy, Derivative)]
#[derivative(PartialEq)]
pub struct Cone {
    pub(super) transform: Matrix<4, 4>,
    pub(super) material: Material,
    pub(super) min: f32,
    pub(super) max: f32,
    pub(super) closed: bool,
    #[derivative(PartialEq = "ignore")]
    pub(super) parent: Group,
}

impl Cone {
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

    pub fn min(mut self, min: f32) -> Self {
        self.min = min;
        self
    }

    pub fn max(mut self, max: f32) -> Self {
        self.max = max;
        self
    }

    pub fn closed(mut self, closed: bool) -> Self {
        self.closed = closed;
        self
    }

    /// A helper to reduce duplication.
    ///
    /// Checks to see if the intersection at `t` is within a radius of 1 from the axis.
    fn check_cap(ray: Ray, t: f32, y: f32) -> bool {
        let x = ray.origin.x() + t * ray.direction.x();
        let z = ray.origin.z() + t * ray.direction.z();
        (x.powi(2) + z.powi(2)) <= y.abs()
    }

    fn intersect_caps(&self, ray: Ray, xs: &mut [Option<Intersection>; 4]) {
        // caps only matters if the cone is closed, and might possibly be
        // intersected by the ray
        if !self.closed || ray.direction.y().abs() < Self::EPSILON {
            return;
        }

        // check for an intersection with the lower end cap by intersecting the ray with
        // the plane at y=self.min
        let t = (self.min - ray.origin.y()) / ray.direction.y();
        if Self::check_cap(ray, t, self.min) {
            if let Some(xs) = xs.iter_mut().find(|xs| xs.is_none()) {
                *xs = Some(Intersection::new(t, self.into()));
            }
        }

        // check for an intersection with the upper end cap by intersecting the ray with
        // the plane at y=self.max
        let t = (self.max - ray.origin.y()) / ray.direction.y();
        if Self::check_cap(ray, t, self.max) {
            if let Some(xs) = xs.iter_mut().find(|xs| xs.is_none()) {
                *xs = Some(Intersection::new(t, self.into()));
            }
        }
    }

    pub fn local_intersect(&self, local_ray: Ray) -> [Option<Intersection>; 4] {
        let a = local_ray.direction.x().powi(2) - local_ray.direction.y().powi(2)
            + local_ray.direction.z().powi(2);
        let b = 2. * local_ray.origin.x() * local_ray.direction.x()
            - 2. * local_ray.origin.y() * local_ray.direction.y()
            + 2. * local_ray.origin.z() * local_ray.direction.z();
        let c = local_ray.origin.x().powi(2) - local_ray.origin.y().powi(2)
            + local_ray.origin.z().powi(2);

        if a.abs() < Self::EPSILON {
            if b.abs() < Self::EPSILON {
                return [None; 4];
            } else {
                let mut xs = [None; 4];
                let t = -c / (2. * b);
                xs[0] = Some(Intersection::new(t, self.into()));
                self.intersect_caps(local_ray, &mut xs);
                return xs;
            }
        }

        let disc = b.powi(2) - 4. * a * c;
        if disc < 0. {
            return [None; 4];
        }

        let mut t0 = (-b - disc.sqrt()) / (2. * a);
        let mut t1 = (-b + disc.sqrt()) / (2. * a);
        if t0 > t1 {
            mem::swap(&mut t0, &mut t1);
        }

        let mut xs = [None; 4];
        let y0 = local_ray.origin.y() + t0 * local_ray.direction.y();
        if self.min < y0 && y0 < self.max {
            xs[0] = Some(Intersection::new(t0, self.into()))
        }
        let y1 = local_ray.origin.y() + t1 * local_ray.direction.y();
        if self.min < y1 && y1 < self.max {
            xs[1] = Some(Intersection::new(t1, self.into()))
        }
        self.intersect_caps(local_ray, &mut xs);
        xs
    }

    pub fn local_normal_at(&self, local_point: Point) -> Vector {
        let dist = local_point.x().powi(2) + local_point.z().powi(2);
        if dist < 1. && local_point.y() >= self.max - Self::EPSILON {
            Vector::new(0., 1., 0.)
        } else if dist < 1. && local_point.y() <= self.min + Self::EPSILON {
            Vector::new(0., -1., 0.)
        } else {
            let mut y = (local_point.x().powi(2) + local_point.z().powi(2)).sqrt();
            if local_point.y() > 0. {
                y = -y;
            }
            Vector::new(local_point.x(), y, local_point.z())
        }
    }
}

impl Default for Cone {
    fn default() -> Self {
        Self {
            transform: Matrix::identity(),
            material: Material::default(),
            min: f32::NEG_INFINITY,
            max: f32::INFINITY,
            closed: false,
            parent: Group::null(),
        }
    }
}

impl From<Cone> for Shape {
    fn from(c: Cone) -> Self {
        Self::Cone(c)
    }
}

impl From<&Cone> for Shape {
    fn from(c: &Cone) -> Self {
        Self::Cone(*c)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tuples::Vector;

    #[test]
    fn cone_basics() {
        let cases = vec![
            (Point::new(0., 0., -5.), Vector::new(0., 0., 1.), 5., 5.),
            // (
            //     Point::new(0., 0., -5.),
            //     Vector::new(1., 1., 1.),
            //     8.66025,
            //     8.66025,
            // ),
            (
                Point::new(1., 1., -5.),
                Vector::new(-0.5, -1., 1.),
                4.55006,
                49.44994,
            ),
        ];
        for (origin, direction, t1, t2) in cases.into_iter() {
            let c = Cone::default();
            let direction = direction.normalize();
            let r = Ray::new(origin, direction);
            let xs = c
                .local_intersect(r)
                .into_iter()
                .filter_map(|xs| xs)
                .collect::<Vec<_>>();
            assert!((xs[0].t() - t1).abs() < Cone::EPSILON);
            assert!((xs[1].t() - t2).abs() < Cone::EPSILON);
        }

        let c = Cone::default();
        let direction = Vector::new(0., 1., 1.).normalize();
        let r = Ray::new(Point::new(0., 0., -1.), direction);
        let xs = c.local_intersect(r);
        assert!((xs[0].unwrap().t() - 0.35355).abs() < Cone::EPSILON);

        let cases = vec![
            (Point::new(0., 0., -5.), Vector::new(0., 1., 0.), 0),
            (Point::new(0., 0., -0.25), Vector::new(0., 1., 1.), 2),
            (Point::new(0., 0., -0.25), Vector::new(0., 1., 0.), 4),
        ];
        for (origin, direction, count) in cases.into_iter() {
            let c = Cone::default().min(-0.5).max(0.5).closed(true);
            let direction = direction.normalize();
            let r = Ray::new(origin, direction);
            assert!(c.local_intersect(r).into_iter().filter_map(|xs| xs).count() == count);
        }

        let cases = vec![
            (Point::new(0., 0., 0.), Vector::new(0., 0., 0.)),
            (Point::new(1., 1., 1.), Vector::new(1., -f32::sqrt(2.), 1.)),
            (Point::new(-1., -1., 0.), Vector::new(-1., 1., 0.)),
        ];
        for (point, normal) in cases.into_iter() {
            let c = Cone::default();
            assert!(c.local_normal_at(point) == normal);
        }
    }
}
