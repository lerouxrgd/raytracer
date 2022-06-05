use std::cmp;
use std::collections::BTreeMap;
use std::ops::Index;

use ordered_float::OrderedFloat;

use crate::materials::Material;
use crate::rays::Ray;
use crate::spheres::Sphere;
use crate::tuples::{Point, Vector};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Shape {
    Sphere(Sphere),
}

impl Shape {
    pub fn normal_at(&self, world_point: Point) -> Vector {
        match self {
            Self::Sphere(s) => s.normal_at(world_point),
        }
    }

    pub fn material(&self) -> Material {
        match self {
            Self::Sphere(s) => s.material,
        }
    }

    pub fn material_mut(&mut self) -> &mut Material {
        match self {
            Self::Sphere(s) => &mut s.material,
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

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Intersection {
    t: OrderedFloat<f32>,
    object: Shape,
}

impl Eq for Intersection {}

impl Ord for Intersection {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        self.t.cmp(&other.t)
    }
}

impl PartialOrd for Intersection {
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Intersection {
    pub fn new(t: f32, object: Shape) -> Self {
        Self {
            t: t.into(),
            object,
        }
    }

    pub fn t(&self) -> f32 {
        self.t.into()
    }

    pub fn object(&self) -> Shape {
        self.object
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Intersections(BTreeMap<Intersection, ()>);

impl From<Vec<Intersection>> for Intersections {
    fn from(intersections: Vec<Intersection>) -> Self {
        Self(intersections.into_iter().map(|i| (i, ())).collect())
    }
}

impl Index<usize> for Intersections {
    type Output = Intersection;

    fn index(&self, i: usize) -> &Self::Output {
        self.0.keys().nth(i).unwrap()
    }
}

impl Intersections {
    pub fn count(&self) -> usize {
        self.0.len()
    }

    pub fn hit(&self) -> Option<Intersection> {
        self.0.keys().find(|i| i.t() > 0.).copied()
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Computations {
    pub t: f32,
    pub object: Shape,
    pub point: Point,
    pub eyev: Vector,
    pub normalv: Vector,
    pub inside: bool,
}

impl Computations {
    pub fn prepare(intersection: Intersection, ray: Ray) -> Self {
        let world_point = ray.position(intersection.t());

        let t = intersection.t();
        let object = intersection.object;
        let point = world_point;
        let eyev = -ray.direction;
        let mut normalv = intersection.object.normal_at(world_point);

        let inside = if normalv.dot(eyev) < 0. {
            normalv = -normalv;
            true
        } else {
            false
        };

        Self {
            t,
            object,
            point,
            eyev,
            normalv,
            inside,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intersections_sphere() {
        let s = Sphere::new();
        let i = Intersection::new(3.5, s.into());
        assert!(i.object() == Shape::Sphere(s));

        let s = Sphere::new();
        let i1 = Intersection::new(1., s.into());
        let i2 = Intersection::new(2., s.into());
        let xs = Intersections::from(vec![i1, i2]);
        assert!(xs.count() == 2);
        assert!(xs[0].t == 1.);
        assert!(xs[1].t == 2.);

        let s = Sphere::new();
        let i1 = Intersection::new(1., s.into());
        let i2 = Intersection::new(2., s.into());
        let xs = Intersections::from(vec![i2, i1]);
        assert!(xs.hit().unwrap() == i1);

        let s = Sphere::new();
        let i1 = Intersection::new(-1., s.into());
        let i2 = Intersection::new(2., s.into());
        let xs = Intersections::from(vec![i2, i1]);
        assert!(xs.hit().unwrap() == i2);

        let s = Sphere::new();
        let i1 = Intersection::new(-2., s.into());
        let i2 = Intersection::new(-1., s.into());
        let xs = Intersections::from(vec![i2, i1]);
        assert!(xs.hit().is_none());

        let s = Sphere::new();
        let i1 = Intersection::new(5., s.into());
        let i2 = Intersection::new(7., s.into());
        let i3 = Intersection::new(-3., s.into());
        let i4 = Intersection::new(2., s.into());
        let xs = Intersections::from(vec![i1, i2, i3, i4]);
        assert!(xs.hit().unwrap() == i4);
    }

    #[test]
    fn intersection_precomputation() {
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let shape = Sphere::new();
        let i = Intersection::new(4.0, shape.into());
        let comps = Computations::prepare(i, r);
        assert!(comps.point == Point::new(0., 0., -1.));
        assert!(comps.eyev == Vector::new(0., 0., -1.));
        assert!(comps.normalv == Vector::new(0., 0., -1.));
        assert!(comps.inside == false);

        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let shape = Sphere::new();
        let i = Intersection::new(1.0, shape.into());
        let comps = Computations::prepare(i, r);
        assert!(comps.point == Point::new(0., 0., 1.));
        assert!(comps.eyev == Vector::new(0., 0., -1.));
        assert!(comps.normalv == Vector::new(0., 0., -1.));
        assert!(comps.inside == true);
    }
}
