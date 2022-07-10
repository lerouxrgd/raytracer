use ordered_float::OrderedFloat;

use crate::intersections::Intersection;
use crate::materials::Material;
use crate::matrices::Matrix;
use crate::rays::Ray;
use crate::tuples::{Point, Tuple, Vector};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Shape {
    Sphere(Sphere),
    Plane(Plane),
    Cube(Cube),
}

impl Shape {
    pub fn transform(&self) -> Matrix<4, 4> {
        match self {
            &Self::Sphere(Sphere { transform, .. })
            | &Self::Plane(Plane { transform, .. })
            | &Self::Cube(Cube { transform, .. }) => transform,
        }
    }

    pub fn transform_mut(&mut self) -> &mut Matrix<4, 4> {
        match self {
            Self::Sphere(Sphere {
                ref mut transform, ..
            })
            | Self::Plane(Plane {
                ref mut transform, ..
            })
            | Self::Cube(Cube {
                ref mut transform, ..
            }) => transform,
        }
    }

    pub fn material(&self) -> Material {
        match self {
            &Self::Sphere(Sphere { material, .. })
            | &Self::Plane(Plane { material, .. })
            | &Self::Cube(Cube { material, .. }) => material,
        }
    }

    pub fn material_mut(&mut self) -> &mut Material {
        match self {
            Self::Sphere(Sphere {
                ref mut material, ..
            })
            | Self::Plane(Plane {
                ref mut material, ..
            })
            | Self::Cube(Cube {
                ref mut material, ..
            }) => material,
        }
    }

    pub fn normal_at(&self, world_point: Point) -> Vector {
        let inverse = self.transform().inverse().unwrap();
        let local_point = inverse * world_point;
        let local_normal = match self {
            Self::Sphere(s) => s.local_normal_at(local_point),
            Self::Plane(p) => p.local_normal_at(local_point),
            Self::Cube(c) => c.local_normal_at(local_point),
        };
        let world_normal = Tuple::from(inverse.transpose() * local_normal);
        let world_normal = Vector::new(world_normal[0], world_normal[1], world_normal[2]);
        world_normal.normalize()
    }

    pub fn intersect(&self, intersections: &mut Vec<Intersection>, world_ray: Ray) {
        let inverse = self.transform().inverse().unwrap();
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
            Shape::Cube(s) => {
                if let Some(xs) = s.local_intersect(local_ray) {
                    intersections.extend(xs)
                }
            }
        }
    }
}

impl From<Sphere> for Shape {
    fn from(s: Sphere) -> Self {
        Self::Sphere(s)
    }
}

impl From<&Sphere> for Shape {
    fn from(s: &Sphere) -> Self {
        Self::Sphere(*s)
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

/// A unit sphere centered on the origin
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sphere {
    pub transform: Matrix<4, 4>,
    pub material: Material,
}

impl Sphere {
    pub fn new() -> Self {
        Self {
            transform: Matrix::identity(),
            material: Material::default(),
        }
    }

    pub fn local_intersect(&self, local_ray: Ray) -> Option<[Intersection; 2]> {
        let sphere_to_ray = local_ray.origin - Point::new(0., 0., 0.);
        let a = local_ray.direction.dot(local_ray.direction);
        let b = 2. * local_ray.direction.dot(sphere_to_ray);
        let c = sphere_to_ray.dot(sphere_to_ray) - 1.;
        let discriminant = b.powi(2) - 4. * a * c;
        if discriminant < 0. {
            None
        } else {
            let t1 = (-b - discriminant.sqrt()) / (2. * a);
            let t2 = (-b + discriminant.sqrt()) / (2. * a);
            Some([
                Intersection::new(t1, self.into()),
                Intersection::new(t2, self.into()),
            ])
        }
    }

    pub fn local_normal_at(&self, local_point: Point) -> Vector {
        local_point - Point::new(0., 0., 0.)
    }
}

impl Default for Sphere {
    fn default() -> Self {
        Self::new()
    }
}

/// A xz plane (left-handed coordinates)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane {
    pub transform: Matrix<4, 4>,
    pub material: Material,
}

impl Plane {
    pub const EPSILON: f32 = 1e-4;

    pub fn new() -> Self {
        Self {
            transform: Matrix::identity(),
            material: Material::default(),
        }
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
}

impl Default for Plane {
    fn default() -> Self {
        Self::new()
    }
}

/// An axis-aligned bounding box (AABB) centered on the origin
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cube {
    pub transform: Matrix<4, 4>,
    pub material: Material,
}

impl Cube {
    pub const EPSILON: f32 = 1e-4;

    pub fn new() -> Self {
        Self {
            transform: Matrix::identity(),
            material: Material::default(),
        }
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
        let maxc = [local_point.x(), local_point.y(), local_point.z()]
            .into_iter()
            .map(|v| OrderedFloat(v.abs()))
            .max()
            .unwrap();

        if local_point.x().abs() == maxc.into() {
            return Vector::new(local_point.x(), 0., 0.);
        } else if local_point.y().abs() == maxc.into() {
            return Vector::new(0., local_point.y(), 0.);
        } else if local_point.z().abs() == maxc.into() {
            return Vector::new(0., 0., local_point.z());
        } else {
            unreachable!()
        }
    }
}

impl Default for Cube {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transformations::*;
    use crate::tuples::Vector;
    use std::f32::consts::PI;

    #[test]
    fn sphere_ray_intersect() {
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.local_intersect(r);
        assert!(xs.unwrap()[0].object() == Shape::Sphere(s));
        assert!(xs.unwrap()[1].object() == Shape::Sphere(s));

        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.local_intersect(r);
        assert!(xs.unwrap()[0].t() == 4.);
        assert!(xs.unwrap()[1].t() == 6.);

        let r = Ray::new(Point::new(0., 1., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.local_intersect(r);
        assert!(xs.unwrap()[0].t() == 5.);
        assert!(xs.unwrap()[1].t() == 5.);

        let r = Ray::new(Point::new(0., 2., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.local_intersect(r);
        assert!(xs.is_none());

        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.local_intersect(r);
        assert!(xs.unwrap()[0].t() == -1.);
        assert!(xs.unwrap()[1].t() == 1.);

        let r = Ray::new(Point::new(0., 0., 5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.local_intersect(r);
        assert!(xs.unwrap()[0].t() == -6.);
        assert!(xs.unwrap()[1].t() == -4.);

        let mut s = Sphere::new();
        s.transform = scaling(2., 2., 2.);
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let mut xs = vec![];
        Shape::Sphere(s).intersect(&mut xs, r);
        assert!(xs[0].t() == 3.);
        assert!(xs[1].t() == 7.);

        let mut s = Sphere::new();
        s.transform = translation(5., 0., 0.);
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let mut xs = vec![];
        Shape::Sphere(s).intersect(&mut xs, r);
        assert!(xs.is_empty());
    }

    #[test]
    fn sphere_normal() {
        let s = Sphere::new();
        assert!(s.local_normal_at(Point::new(1., 0., 0.)) == Vector::new(1., 0., 0.));
        assert!(s.local_normal_at(Point::new(0., 1., 0.)) == Vector::new(0., 1., 0.));
        assert!(s.local_normal_at(Point::new(0., 0., 1.)) == Vector::new(0., 0., 1.));

        let s = Sphere::new();
        let p = Point::new(f32::sqrt(3.) / 3., f32::sqrt(3.) / 3., f32::sqrt(3.) / 3.);
        let n = Vector::new(f32::sqrt(3.) / 3., f32::sqrt(3.) / 3., f32::sqrt(3.) / 3.);
        assert!(s.local_normal_at(p).equal_approx(n));
        assert!(s
            .local_normal_at(p)
            .equal_approx(s.local_normal_at(p).normalize()));

        let mut s = Sphere::new();
        s.transform = translation(0., 1., 0.);
        let p = Point::new(0., 1.70711, -0.70711);
        let n = Shape::Sphere(s).normal_at(p);
        assert!(n.equal_approx(Vector::new(0., 0.70711, -0.70711)));

        let mut s = Sphere::new();
        s.transform = Transform::new()
            .roation_z(PI / 5.)
            .scaling(1., 0.5, 1.)
            .into();
        let p = Point::new(0., f32::sqrt(2.) / 2., -f32::sqrt(2.) / 2.);
        let n = Shape::Sphere(s).normal_at(p);
        assert!(n.equal_approx(Vector::new(0., 0.97014, -0.24254)));
    }

    #[test]
    fn plane_basics() {
        let p = Plane::new();
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
            let c = Cube::new();
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
            let c = Cube::new();
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
            let c = Cube::new();
            assert!(c.local_normal_at(point) == normal);
        }
    }
}
