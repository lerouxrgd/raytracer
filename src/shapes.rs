use crate::intersections::Intersection;
use crate::materials::Material;
use crate::matrices::Matrix;
use crate::rays::Ray;
use crate::tuples::{Point, Tuple, Vector};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Shape {
    Sphere(Sphere),
    Plane(Plane),
}

impl Shape {
    pub fn transform(&self) -> Matrix<4, 4> {
        match self {
            &Self::Sphere(Sphere { transform, .. }) | &Self::Plane(Plane { transform, .. }) => {
                transform
            }
        }
    }

    pub fn transform_mut(&mut self) -> &mut Matrix<4, 4> {
        match self {
            Self::Sphere(Sphere {
                ref mut transform, ..
            })
            | Self::Plane(Plane {
                ref mut transform, ..
            }) => transform,
        }
    }

    pub fn material(&self) -> Material {
        match self {
            &Self::Sphere(Sphere { material, .. }) | &Self::Plane(Plane { material, .. }) => {
                material
            }
        }
    }

    pub fn material_mut(&mut self) -> &mut Material {
        match self {
            Self::Sphere(Sphere {
                ref mut material, ..
            })
            | Self::Plane(Plane {
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
    fn sphere_material() {
        let s = Sphere::new();
        let m = Material::default();
        assert!(s.material == m);
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
}
