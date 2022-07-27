use std::mem;

use derivative::Derivative;
use ordered_float::OrderedFloat;
use slotmap::Key;

use crate::groups::Group;
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
    Cylinder(Cylinder),
    Cone(Cone),
}

impl Shape {
    pub fn transform(&self) -> Matrix<4, 4> {
        match self {
            &Self::Sphere(Sphere { transform, .. })
            | &Self::Plane(Plane { transform, .. })
            | &Self::Cylinder(Cylinder { transform, .. })
            | &Self::Cone(Cone { transform, .. })
            | &Self::Cube(Cube { transform, .. }) => transform,
        }
    }

    pub fn material(&self) -> Material {
        match self {
            &Self::Sphere(Sphere { material, .. })
            | &Self::Plane(Plane { material, .. })
            | &Self::Cylinder(Cylinder { material, .. })
            | &Self::Cone(Cone { material, .. })
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
            | Self::Cylinder(Cylinder {
                ref mut material, ..
            })
            | Self::Cone(Cone {
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
            Self::Cube(c) => c.local_normal_at(local_point),
        };
        self.normal_to_world(local_normal)
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
            | &Self::Cube(Cube { parent, .. }) => parent,
        }
    }

    pub(crate) fn parent_mut(&mut self) -> &mut Group {
        match self {
            Self::Sphere(Sphere { ref mut parent, .. })
            | Self::Plane(Plane { ref mut parent, .. })
            | Self::Cylinder(Cylinder { ref mut parent, .. })
            | Self::Cone(Cone { ref mut parent, .. })
            | Self::Cube(Cube { ref mut parent, .. }) => parent,
        }
    }

    pub fn world_to_object(&self, point: Point) -> Point {
        let point = if !self.parent().is_null() {
            self.parent().world_to_object(point)
        } else {
            point
        };
        self.transform().inverse().unwrap() * point
    }

    pub fn normal_to_world(&self, normal: Vector) -> Vector {
        let inverse = self.transform().inverse().unwrap();
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

impl From<Cylinder> for Shape {
    fn from(c: Cylinder) -> Self {
        Self::Cylinder(c)
    }
}

impl From<&Cylinder> for Shape {
    fn from(c: &Cylinder) -> Self {
        Self::Cylinder(*c)
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

/// A unit sphere centered on the origin
#[derive(Debug, Clone, Copy, Derivative)]
#[derivative(PartialEq)]
pub struct Sphere {
    pub transform: Matrix<4, 4>,
    pub material: Material,
    #[derivative(PartialEq = "ignore")]
    parent: Group,
}

impl Sphere {
    pub fn new() -> Self {
        Self {
            transform: Matrix::identity(),
            material: Material::default(),
            parent: Group::null(),
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
#[derive(Debug, Clone, Copy, Derivative)]
#[derivative(PartialEq)]
pub struct Plane {
    pub transform: Matrix<4, 4>,
    pub material: Material,
    #[derivative(PartialEq = "ignore")]
    parent: Group,
}

impl Plane {
    pub const EPSILON: f32 = 1e-4;

    pub fn new() -> Self {
        Self {
            transform: Matrix::identity(),
            material: Material::default(),
            parent: Group::null(),
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
#[derive(Debug, Clone, Copy, Derivative)]
#[derivative(PartialEq)]
pub struct Cube {
    pub transform: Matrix<4, 4>,
    pub material: Material,
    #[derivative(PartialEq = "ignore")]
    parent: Group,
}

impl Cube {
    pub const EPSILON: f32 = 1e-4;

    pub fn new() -> Self {
        Self {
            transform: Matrix::identity(),
            material: Material::default(),
            parent: Group::null(),
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
            Vector::new(local_point.x(), 0., 0.)
        } else if local_point.y().abs() == maxc.into() {
            Vector::new(0., local_point.y(), 0.)
        } else if local_point.z().abs() == maxc.into() {
            Vector::new(0., 0., local_point.z())
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

/// Infinite cynlinder of radius 1 along the y axis
#[derive(Debug, Clone, Copy, Derivative)]
#[derivative(PartialEq)]
pub struct Cylinder {
    pub transform: Matrix<4, 4>,
    pub material: Material,
    pub min: f32,
    pub max: f32,
    pub closed: bool,
    #[derivative(PartialEq = "ignore")]
    parent: Group,
}

impl Cylinder {
    pub const EPSILON: f32 = 1e-4;

    pub fn new() -> Self {
        Self {
            transform: Matrix::identity(),
            material: Material::default(),
            min: f32::NEG_INFINITY,
            max: f32::INFINITY,
            closed: false,
            parent: Group::null(),
        }
    }

    /// A helper to reduce duplication.
    ///
    /// Checks to see if the intersection at `t` is within a radius of 1 from the axis.
    fn check_cap(ray: Ray, t: f32) -> bool {
        let x = ray.origin.x() + t * ray.direction.x();
        let z = ray.origin.z() + t * ray.direction.z();
        (x.powi(2) + z.powi(2)) <= 1.
    }

    fn intersect_caps(&self, ray: Ray, xs: &mut [Option<Intersection>; 2]) {
        // caps only matters if the cylinder is closed, and might possibly be
        // intersected by the ray
        if !self.closed || ray.direction.y().abs() < Self::EPSILON {
            return;
        }

        // check for an intersection with the lower end cap by intersecting the ray with
        // the plane at y=self.min
        let t = (self.min - ray.origin.y()) / ray.direction.y();
        if Self::check_cap(ray, t) {
            if let Some(xs) = xs.iter_mut().find(|xs| xs.is_none()) {
                *xs = Some(Intersection::new(t, self.into()));
            }
        }

        // check for an intersection with the upper end cap by intersecting the ray with
        // the plane at y=self.max
        let t = (self.max - ray.origin.y()) / ray.direction.y();
        if Self::check_cap(ray, t) {
            if let Some(xs) = xs.iter_mut().find(|xs| xs.is_none()) {
                *xs = Some(Intersection::new(t, self.into()));
            }
        }
    }

    pub fn local_intersect(&self, local_ray: Ray) -> [Option<Intersection>; 2] {
        let a = local_ray.direction.x().powi(2) + local_ray.direction.z().powi(2);
        if a.abs() < Self::EPSILON {
            // ray is parallel to the y axis
            let mut xs = [None; 2];
            self.intersect_caps(local_ray, &mut xs);
            return xs;
        }

        let b = 2. * local_ray.origin.x() * local_ray.direction.x()
            + 2. * local_ray.origin.z() * local_ray.direction.z();
        let c = local_ray.origin.x().powi(2) + local_ray.origin.z().powi(2) - 1.;
        let disc = b.powi(2) - 4. * a * c;
        if disc < 0. {
            // ray does not intersect the cylinder
            return [None, None];
        }

        let mut t0 = (-b - disc.sqrt()) / (2. * a);
        let mut t1 = (-b + disc.sqrt()) / (2. * a);
        if t0 > t1 {
            mem::swap(&mut t0, &mut t1);
        }

        let mut xs = [None; 2];
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
            Vector::new(local_point.x(), 0., local_point.z())
        }
    }
}

impl Default for Cylinder {
    fn default() -> Self {
        Self::new()
    }
}

/// A double-napped cone along the y axis
#[derive(Debug, Clone, Copy, Derivative)]
#[derivative(PartialEq)]
pub struct Cone {
    pub transform: Matrix<4, 4>,
    pub material: Material,
    pub min: f32,
    pub max: f32,
    pub closed: bool,
    #[derivative(PartialEq = "ignore")]
    parent: Group,
}

impl Cone {
    pub const EPSILON: f32 = 1e-4;

    pub fn new() -> Self {
        Self {
            transform: Matrix::identity(),
            material: Material::default(),
            min: f32::NEG_INFINITY,
            max: f32::INFINITY,
            closed: false,
            parent: Group::null(),
        }
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
            .rotation_z(PI / 5.)
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

    #[test]
    fn cylinder_basics() {
        let cases = vec![
            (Point::new(1., 0., 0.), Vector::new(0., 1., 0.)),
            (Point::new(0., 0., 0.), Vector::new(0., 1., 0.)),
            (Point::new(0., 0., -5.), Vector::new(1., 1., 1.)),
        ];
        for (origin, direction) in cases.into_iter() {
            let c = Cylinder::new();
            let direction = direction.normalize();
            let r = Ray::new(origin, direction);
            assert!(c.local_intersect(r) == [None, None]);
        }

        let cases = vec![
            (Point::new(1., 0., -5.), Vector::new(0., 0., 1.), 5., 5.),
            (Point::new(0., 0., -5.), Vector::new(0., 0., 1.), 4., 6.),
            (
                Point::new(0.5, 0., -5.),
                Vector::new(0.1, 1., 1.),
                6.80798,
                7.08872,
            ),
        ];
        for (origin, direction, t1, t2) in cases.into_iter() {
            let c = Cylinder::new();
            let direction = direction.normalize();
            let r = Ray::new(origin, direction);
            let xs = c
                .local_intersect(r)
                .into_iter()
                .filter_map(|xs| xs)
                .collect::<Vec<_>>();
            assert!((xs[0].t() - t1).abs() < Cylinder::EPSILON);
            assert!((xs[1].t() - t2).abs() < Cylinder::EPSILON);
        }

        let cases = vec![
            (Point::new(1., 0., 0.), Vector::new(1., 0., 0.)),
            (Point::new(0., 5., -1.), Vector::new(0., 0., -1.)),
            (Point::new(0., -2., 1.), Vector::new(0., 0., 1.)),
            (Point::new(-1., 1., 0.), Vector::new(-1., 0., 0.)),
        ];
        for (point, normal) in cases.into_iter() {
            let c = Cylinder::new();
            assert!(c.local_normal_at(point) == normal);
        }
    }

    #[test]
    fn cylinder_advanced() {
        let cases = vec![
            (Point::new(0., 1.5, 0.), Vector::new(0.1, 1., 0.), 0),
            (Point::new(0., 3., -5.), Vector::new(0., 0., 1.), 0),
            (Point::new(0., 0., -5.), Vector::new(0., 0., 1.), 0),
            (Point::new(0., 2., -5.), Vector::new(0., 0., 1.), 0),
            (Point::new(0., 1., -5.), Vector::new(0., 0., 1.), 0),
            (Point::new(0., 1.5, -2.), Vector::new(0., 0., 1.), 2),
        ];
        for (origin, direction, count) in cases.into_iter() {
            let mut c = Cylinder::new();
            c.min = 1.;
            c.max = 2.;
            let direction = direction.normalize();
            let r = Ray::new(origin, direction);
            assert!(c.local_intersect(r).into_iter().filter_map(|xs| xs).count() == count);
        }

        let cases = vec![
            (Point::new(0., 3., 0.), Vector::new(0., -1., 0.), 2),
            (Point::new(0., 3., -2.), Vector::new(0., -1., 2.), 2),
            // (Point::new(0., 4., -2.), Vector::new(0., -1., 1.), 2),
            (Point::new(0., 0., -2.), Vector::new(0., 1., 2.), 2),
            // (Point::new(0., -1., -2.), Vector::new(0., 1., 1.), 2),
        ];
        for (point, direction, count) in cases.into_iter() {
            let mut c = Cylinder::new();
            c.min = 1.;
            c.max = 2.;
            c.closed = true;
            let direction = direction.normalize();
            let r = Ray::new(point, direction);
            assert!(c.local_intersect(r).into_iter().filter_map(|xs| xs).count() == count);
        }

        let cases = vec![
            (Point::new(0., 1., 0.), Vector::new(0., -1., 0.)),
            (Point::new(0.5, 1., 0.), Vector::new(0., -1., 0.)),
            (Point::new(0., 1., 0.5), Vector::new(0., -1., 0.)),
            (Point::new(0., 2., 0.), Vector::new(0., 1., 0.)),
            (Point::new(0.5, 2., 0.), Vector::new(0., 1., 0.)),
            (Point::new(0., 2., 0.5), Vector::new(0., 1., 0.)),
        ];
        for (point, normal) in cases.into_iter() {
            let mut c = Cylinder::new();
            c.min = 1.;
            c.max = 2.;
            c.closed = true;
            assert!(c.local_normal_at(point) == normal);
        }
    }

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
            let c = Cone::new();
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

        let c = Cone::new();
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
            let mut c = Cone::new();
            c.min = -0.5;
            c.max = 0.5;
            c.closed = true;
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
            let c = Cone::new();
            assert!(c.local_normal_at(point) == normal);
        }
    }
}
