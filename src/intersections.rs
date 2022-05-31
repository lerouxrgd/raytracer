use std::cmp;
use std::collections::BTreeMap;
use std::ops::Index;

use ordered_float::OrderedFloat;

use crate::spheres::Sphere;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Intersectable {
    Sphere(Sphere),
}

impl From<Sphere> for Intersectable {
    fn from(s: Sphere) -> Self {
        Self::Sphere(s)
    }
}

impl From<&Sphere> for Intersectable {
    fn from(s: &Sphere) -> Self {
        Self::Sphere(*s)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Intersection {
    t: OrderedFloat<f32>,
    object: Intersectable,
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
    pub fn new(t: f32, object: Intersectable) -> Self {
        Self {
            t: t.into(),
            object,
        }
    }

    pub fn t(&self) -> f32 {
        self.t.into()
    }

    pub fn object(&self) -> Intersectable {
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
        self.0.keys().find(|i| i.t() > 0.).map(|&i| i)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intersections_sphere() {
        let s = Sphere::new();
        let i = Intersection::new(3.5, s.into());
        assert!(i.object() == Intersectable::Sphere(s));

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
}
