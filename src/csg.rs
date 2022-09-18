use crate::bounds::BoundingBox;
use crate::intersections::{Intersection, Intersections};
use crate::rays::Ray;
use crate::shapes::{Cone, Cube, Cylinder, Plane, Shape, SmoothTriangle, Triangle};

#[derive(Debug, Clone, PartialEq)]
pub struct Csg {
    op: CsgOp,
    left: CsgChild,
    right: CsgChild,
}

impl Csg {
    pub fn new<L, R>(op: CsgOp, left: L, right: R) -> Self
    where
        L: Into<CsgChild>,
        R: Into<CsgChild>,
    {
        Self {
            op,
            left: left.into(),
            right: right.into(),
        }
    }

    pub fn local_intersect(&self, ray: Ray) -> Vec<Intersection> {
        if !self.bounds().intersects(ray) {
            return vec![];
        }

        let mut xs = vec![];
        match (&self.left, &self.right) {
            (CsgChild::Shape(shape_l), CsgChild::Shape(shape_r)) => {
                shape_l.intersect(&mut xs, ray);
                shape_r.intersect(&mut xs, ray);
            }
            (CsgChild::Shape(shape), CsgChild::Csg(csg))
            | (CsgChild::Csg(csg), CsgChild::Shape(shape)) => {
                shape.intersect(&mut xs, ray);
                xs.extend(csg.local_intersect(ray));
            }
            (CsgChild::Csg(csg_l), CsgChild::Csg(csg_r)) => {
                xs.extend(csg_l.local_intersect(ray));
                xs.extend(csg_r.local_intersect(ray));
            }
        }
        let xs = Intersections::from(xs);
        self.filter_intersections(xs).into()
    }

    fn filter_intersections(&self, xs: Intersections) -> Intersections {
        // Begin outside of both children
        let mut in_l = false;
        let mut in_r = false;
        let mut result = vec![];
        for x in Vec::from(xs) {
            let l_hit = match &self.left {
                CsgChild::Shape(shape) => shape.as_ref() == &x.shape(),
                CsgChild::Csg(csg) => csg.includes(&x.shape()),
            };
            if self.op.intersection_allowed(l_hit, in_l, in_r) {
                result.push(x)
            }
            // Depending on which object was hit, toggle either inl or inr
            if l_hit {
                in_l = !in_l;
            } else {
                in_r = !in_r;
            }
        }
        result.into()
    }

    fn includes(&self, s: &Shape) -> bool {
        match (&self.left, &self.right) {
            (CsgChild::Shape(shape_l), CsgChild::Shape(shape_r)) => {
                shape_l.as_ref() == s || shape_r.as_ref() == s
            }
            (CsgChild::Shape(shape), CsgChild::Csg(csg))
            | (CsgChild::Csg(csg), CsgChild::Shape(shape)) => {
                shape.as_ref() == s || csg.includes(s)
            }
            (CsgChild::Csg(csg_l), CsgChild::Csg(csg_r)) => csg_l.includes(s) || csg_r.includes(s),
        }
    }

    pub fn bounds(&self) -> BoundingBox {
        let mut bb = BoundingBox::default();
        match (&self.left, &self.right) {
            (CsgChild::Shape(shape_l), CsgChild::Shape(shape_r)) => {
                bb.add_box(shape_l.parent_space_bounds());
                bb.add_box(shape_r.parent_space_bounds());
            }
            (CsgChild::Shape(shape), CsgChild::Csg(csg))
            | (CsgChild::Csg(csg), CsgChild::Shape(shape)) => {
                bb.add_box(shape.parent_space_bounds());
                bb.add_box(csg.bounds());
            }
            (CsgChild::Csg(csg_l), CsgChild::Csg(csg_r)) => {
                bb.add_box(csg_l.bounds());
                bb.add_box(csg_r.bounds());
            }
        }
        bb
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum CsgOp {
    Union,
    Intersect,
    Difference,
}

impl CsgOp {
    pub fn intersection_allowed(&self, l_hit: bool, in_l: bool, in_r: bool) -> bool {
        match self {
            Self::Union => (l_hit && !in_r) || (!l_hit && !in_l),
            Self::Intersect => (l_hit && in_r) || (!l_hit && in_l),
            Self::Difference => (l_hit && !in_r) || (!l_hit && in_l),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum CsgChild {
    Shape(Box<Shape>),
    Csg(Box<Csg>),
}

impl From<Shape> for CsgChild {
    fn from(shape: Shape) -> Self {
        Self::Shape(Box::new(shape))
    }
}

impl From<Plane> for CsgChild {
    fn from(shape: Plane) -> Self {
        Self::Shape(Box::new(shape.into()))
    }
}

impl From<Cube> for CsgChild {
    fn from(shape: Cube) -> Self {
        Self::Shape(Box::new(shape.into()))
    }
}

impl From<Cylinder> for CsgChild {
    fn from(shape: Cylinder) -> Self {
        Self::Shape(Box::new(shape.into()))
    }
}

impl From<Cone> for CsgChild {
    fn from(shape: Cone) -> Self {
        Self::Shape(Box::new(shape.into()))
    }
}

impl From<Triangle> for CsgChild {
    fn from(shape: Triangle) -> Self {
        Self::Shape(Box::new(shape.into()))
    }
}

impl From<SmoothTriangle> for CsgChild {
    fn from(shape: SmoothTriangle) -> Self {
        Self::Shape(Box::new(shape.into()))
    }
}

impl From<Csg> for CsgChild {
    fn from(csg: Csg) -> Self {
        Self::Csg(Box::new(csg))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transformations::translation;
    use crate::tuples::{Point, Vector};

    #[test]
    fn csg_rules() {
        let cases = vec![
            (CsgOp::Union, true, true, true, false),
            (CsgOp::Union, true, true, false, true),
            (CsgOp::Union, true, false, true, false),
            (CsgOp::Union, true, false, false, true),
            (CsgOp::Union, false, true, true, false),
            (CsgOp::Union, false, true, false, false),
            (CsgOp::Union, false, false, true, true),
            (CsgOp::Union, false, false, false, true),
            (CsgOp::Intersect, true, true, true, true),
            (CsgOp::Intersect, true, true, false, false),
            (CsgOp::Intersect, true, false, true, true),
            (CsgOp::Intersect, true, false, false, false),
            (CsgOp::Intersect, false, true, true, true),
            (CsgOp::Intersect, false, true, false, true),
            (CsgOp::Intersect, false, false, true, false),
            (CsgOp::Intersect, false, false, false, false),
            (CsgOp::Difference, true, true, true, false),
            (CsgOp::Difference, true, true, false, true),
            (CsgOp::Difference, true, false, true, false),
            (CsgOp::Difference, true, false, false, true),
            (CsgOp::Difference, false, true, true, true),
            (CsgOp::Difference, false, true, false, true),
            (CsgOp::Difference, false, false, true, false),
            (CsgOp::Difference, false, false, false, false),
        ];
        for (op, l_hit, in_l, in_r, res) in cases {
            assert!(op.intersection_allowed(l_hit, in_l, in_r) == res)
        }
    }

    #[test]
    fn csg_intersections() {
        let s1 = Shape::sphere();
        let s2 = Shape::cube();
        let xs = vec![
            Intersection::new(1., s1),
            Intersection::new(2., s2),
            Intersection::new(3., s1),
            Intersection::new(4., s2),
        ];
        let cases = vec![
            (CsgOp::Union, xs[0].t(), xs[3].t()),
            (CsgOp::Intersect, xs[1].t(), xs[2].t()),
            (CsgOp::Difference, xs[0].t(), xs[1].t()),
        ];
        for (op, t1, t2) in cases {
            let c = Csg::new(op, s1, s2);
            let res = Vec::from(c.filter_intersections(xs.clone().into()));
            assert!(res.len() == 2);
            assert!(res[0].t() == t1);
            assert!(res[1].t() == t2);
        }

        let c = Csg::new(CsgOp::Union, Shape::sphere(), Shape::cube());
        let r = Ray::new(Point::new(0., 2., -5.), Vector::new(0., 1., 0.));
        let xs = c.local_intersect(r);
        assert!(xs.is_empty());

        let s1 = Shape::sphere();
        let s2 = Shape::sphere().with_transform(translation(0., 0., 0.5));
        let c = Csg::new(CsgOp::Union, s1, s2);
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let xs = c.local_intersect(r);
        assert!(xs.len() == 2);
        assert!(xs[0].t() == 4.);
        assert!(xs[0].shape() == s1);
        assert!(xs[1].t() == 6.5);
        assert!(xs[1].shape() == s2);
    }

    #[test]
    fn csg_bounding_box() {
        let left = Shape::sphere();
        let right = Shape::sphere().with_transform(translation(2., 3., 4.));
        let csg = Csg::new(CsgOp::Difference, left, right);
        let bb = csg.bounds();
        assert!(bb.min == Point::new(-1., -1., -1.));
        assert!(bb.max == Point::new(3., 4., 5.));
    }
}
