use std::cmp;
use std::collections::BTreeMap;
use std::ops::Index;

use ordered_float::OrderedFloat;

use crate::rays::Ray;
use crate::shapes::Shape;
use crate::tuples::{Point, Vector};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Intersection {
    t: OrderedFloat<f32>,
    shape: Shape,
    u: f32,
    v: f32,
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
    pub fn new(t: f32, shape: Shape) -> Self {
        Self {
            t: t.into(),
            shape,
            u: f32::MAX,
            v: f32::MAX,
        }
    }

    pub fn with_uv(mut self, u: f32, v: f32) -> Self {
        self.u = u;
        self.v = v;
        self
    }

    pub fn t(&self) -> f32 {
        self.t.into()
    }

    pub fn shape(&self) -> Shape {
        self.shape
    }

    pub fn u(&self) -> f32 {
        self.u
    }

    pub fn v(&self) -> f32 {
        self.v
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Intersections(BTreeMap<Intersection, ()>);

impl From<Vec<Intersection>> for Intersections {
    fn from(intersections: Vec<Intersection>) -> Self {
        Self(intersections.into_iter().map(|i| (i, ())).collect())
    }
}

impl From<Intersections> for Vec<Intersection> {
    fn from(xs: Intersections) -> Self {
        xs.0.into_keys().collect()
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

    pub fn iter(&self) -> impl Iterator<Item = &Intersection> {
        self.0.keys()
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Computations {
    pub t: f32,
    pub shape: Shape,
    pub point: Point,
    pub over_point: Point,
    pub under_point: Point,
    pub eyev: Vector,
    pub normalv: Vector,
    pub inside: bool,
    pub reflectv: Vector,
    pub n1: f32,
    pub n2: f32,
}

impl Computations {
    pub const EPSILON: f32 = 1e-4;

    pub fn prepare(intersection: Intersection, ray: Ray, xs: &Intersections) -> Self {
        let world_point = ray.position(intersection.t());

        let t = intersection.t();
        let shape = intersection.shape;
        let point = world_point;
        let eyev = -ray.direction;
        let mut normalv = intersection.shape.normal_at(world_point, intersection);

        let inside = if normalv.dot(eyev) < 0. {
            normalv = -normalv;
            true
        } else {
            false
        };

        let reflectv = ray.direction.reflect(normalv);
        let over_point = point + Self::EPSILON * normalv;
        let under_point = point - Self::EPSILON * normalv;

        let mut n1 = 1.;
        let mut n2 = 1.;
        let mut containers: Vec<Shape> = vec![];
        for x in xs.iter() {
            if x == &intersection {
                if let Some(shape) = containers.last() {
                    n1 = shape.get_material().refractive_index;
                }
            }
            if let Some(i) = containers.iter().position(|shape| shape == &x.shape) {
                containers.remove(i);
            } else {
                containers.push(x.shape);
            }
            if x == &intersection {
                if let Some(shape) = containers.last() {
                    n2 = shape.get_material().refractive_index;
                }
            }
        }

        Self {
            t,
            shape,
            point,
            over_point,
            under_point,
            eyev,
            normalv,
            inside,
            reflectv,
            n1,
            n2,
        }
    }

    pub fn schlick(&self) -> f32 {
        let mut cos = self.eyev.dot(self.normalv);

        if self.n1 > self.n2 {
            let n = self.n1 / self.n2;
            let sin2_t = n.powi(2) * (1. - cos.powi(2));
            if sin2_t > 1. {
                return 1.;
            }
            let cos_t = (1. - sin2_t).sqrt();
            cos = cos_t;
        }

        let r0 = ((self.n1 - self.n2) / (self.n1 + self.n2)).powi(2);
        r0 + (1. - r0) * (1. - cos).powi(5)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::materials::Material;
    use crate::shapes::*;
    use crate::transformations::*;

    #[test]
    fn intersections_sphere() {
        let s = Sphere::default();
        let i = Intersection::new(3.5, s.into());
        assert!(i.shape() == Shape::Sphere(s));

        let s = Sphere::default();
        let i1 = Intersection::new(1., s.into());
        let i2 = Intersection::new(2., s.into());
        let xs = Intersections::from(vec![i1, i2]);
        assert!(xs.count() == 2);
        assert!(xs[0].t == 1.);
        assert!(xs[1].t == 2.);

        let s = Sphere::default();
        let i1 = Intersection::new(1., s.into());
        let i2 = Intersection::new(2., s.into());
        let xs = Intersections::from(vec![i2, i1]);
        assert!(xs.hit().unwrap() == i1);

        let s = Sphere::default();
        let i1 = Intersection::new(-1., s.into());
        let i2 = Intersection::new(2., s.into());
        let xs = Intersections::from(vec![i2, i1]);
        assert!(xs.hit().unwrap() == i2);

        let s = Sphere::default();
        let i1 = Intersection::new(-2., s.into());
        let i2 = Intersection::new(-1., s.into());
        let xs = Intersections::from(vec![i2, i1]);
        assert!(xs.hit().is_none());

        let s = Sphere::default();
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
        let shape = Sphere::default();
        let i = Intersection::new(4.0, shape.into());
        let comps = Computations::prepare(i, r, &vec![].into());
        assert!(comps.point == Point::new(0., 0., -1.));
        assert!(comps.eyev == Vector::new(0., 0., -1.));
        assert!(comps.normalv == Vector::new(0., 0., -1.));
        assert!(comps.inside == false);

        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let shape = Sphere::default();
        let i = Intersection::new(1.0, shape.into());
        let comps = Computations::prepare(i, r, &vec![].into());
        assert!(comps.point == Point::new(0., 0., 1.));
        assert!(comps.eyev == Vector::new(0., 0., -1.));
        assert!(comps.normalv == Vector::new(0., 0., -1.));
        assert!(comps.inside == true);
    }

    #[test]
    fn intersection_offset() {
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::default().with_transform(translation(0., 0., 1.));
        let i = Intersection::new(5., s.into());
        let comps = Computations::prepare(i, r, &vec![].into());
        assert!(comps.over_point.z() < -Computations::EPSILON / 2.);
        assert!(comps.point.z() > comps.over_point.z());

        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let shape = Sphere::default()
            .with_transform(translation(0., 0., -1.))
            .with_material(Material::default().transparency(1.).refractive_index(1.5));
        let i = Intersection::new(5., shape.into());
        let xs = vec![Intersection::new(5., shape.into())].into();
        let comps = Computations::prepare(i, r, &xs);
        assert!(comps.under_point.z() > -Computations::EPSILON / 2.);
        assert!(comps.point.z() < comps.under_point.z());
    }

    #[test]
    fn reflection_vector() {
        let p = Plane::default();
        let r = Ray::new(
            Point::new(0., 1., -1.),
            Vector::new(0., -f32::sqrt(2.) / 2., f32::sqrt(2.) / 2.),
        );
        let i = Intersection::new(f32::sqrt(2.) / 2., p.into());
        let comps = Computations::prepare(i, r, &vec![].into());
        assert!(comps.reflectv == Vector::new(0., f32::sqrt(2.) / 2., f32::sqrt(2.) / 2.));
    }

    #[test]
    fn refraction_at_intersections() {
        let a = Sphere::default()
            .with_transform(scaling(2., 2., 2.))
            .with_material(Material::default().transparency(1.).refractive_index(1.5));

        let b = Sphere::default()
            .with_transform(translation(0., 0., -0.25))
            .with_material(Material::default().transparency(1.).refractive_index(2.));

        let c = Sphere::default()
            .with_transform(translation(0., 0., 0.25))
            .with_material(Material::default().transparency(1.).refractive_index(2.5));

        let r = Ray::new(Point::new(0., 0., -4.), Vector::new(0., 0., 1.));
        let xs: Intersections = vec![
            Intersection::new(2., a.into()),
            Intersection::new(2.75, b.into()),
            Intersection::new(3.25, c.into()),
            Intersection::new(4.75, b.into()),
            Intersection::new(5.25, c.into()),
            Intersection::new(6., a.into()),
        ]
        .into();

        let cases = vec![
            (0, 1.0, 1.5),
            (1, 1.5, 2.0),
            (2, 2.0, 2.5),
            (3, 2.5, 2.5),
            (4, 2.5, 1.5),
            (5, 1.5, 1.0),
        ];
        for (i, n1, n2) in cases.into_iter() {
            let comps = Computations::prepare(xs[i], r, &xs);
            assert!(comps.n1 == n1);
            assert!(comps.n2 == n2);
        }
    }

    #[test]
    fn fresnel_effect() {
        let shape = Sphere::default()
            .with_material(Material::default().transparency(1.).refractive_index(1.5));
        let r = Ray::new(
            Point::new(0., 0., f32::sqrt(2.) / 2.),
            Vector::new(0., 1., 0.),
        );
        let xs: Intersections = vec![
            Intersection::new(-f32::sqrt(2.) / 2., shape.into()),
            Intersection::new(f32::sqrt(2.) / 2., shape.into()),
        ]
        .into();
        let comps = Computations::prepare(xs[1], r, &xs);
        assert!(comps.schlick() == 1.);

        let shape = Sphere::default()
            .with_material(Material::default().transparency(1.).refractive_index(1.5));
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 1., 0.));
        let xs: Intersections = vec![
            Intersection::new(-1., shape.into()),
            Intersection::new(1., shape.into()),
        ]
        .into();
        let comps = Computations::prepare(xs[1], r, &xs);
        assert!((comps.schlick() - 0.04).abs() < Computations::EPSILON);

        let shape = Sphere::default()
            .with_material(Material::default().transparency(1.).refractive_index(1.5));
        let r = Ray::new(Point::new(0., 0.99, -2.), Vector::new(0., 0., 1.));
        let xs: Intersections = vec![Intersection::new(1.8589, shape.into())].into();
        let comps = Computations::prepare(xs[0], r, &xs);
        assert!((comps.schlick() - 0.48873).abs() < Computations::EPSILON);
    }
}
