use std::cmp;

use ordered_float::OrderedFloat;

use crate::matrices::Matrix;
use crate::rays::Ray;
use crate::shapes::Cube;
use crate::tuples::Point;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundingBox {
    pub min: Point,
    pub max: Point,
}

impl Default for BoundingBox {
    fn default() -> Self {
        Self {
            min: Point::new(f32::INFINITY, f32::INFINITY, f32::INFINITY),
            max: Point::new(f32::NEG_INFINITY, f32::NEG_INFINITY, f32::NEG_INFINITY),
        }
    }
}

impl BoundingBox {
    pub fn with_min(mut self, min: Point) -> Self {
        self.min = min;
        self
    }

    pub fn with_max(mut self, max: Point) -> Self {
        self.max = max;
        self
    }

    pub fn transform<T: Into<Matrix<4, 4>>>(&self, transform: T) -> Self {
        let t = transform.into();

        let p1 = self.min;
        let p2 = Point::new(self.min.x(), self.min.y(), self.max.z());
        let p3 = Point::new(self.min.x(), self.max.y(), self.min.z());
        let p4 = Point::new(self.min.x(), self.max.y(), self.max.z());
        let p5 = Point::new(self.max.x(), self.min.y(), self.min.z());
        let p6 = Point::new(self.max.x(), self.min.y(), self.max.z());
        let p7 = Point::new(self.max.x(), self.max.y(), self.min.z());
        let p8 = self.max;

        let mut bb = Self::default();
        for p in [p1, p2, p3, p4, p5, p6, p7, p8] {
            bb.add_point(t * p);
        }
        bb
    }

    pub fn add_point(&mut self, p: Point) {
        if p.x() < self.min.x() {
            self.min = self.min.with_x(p.x());
        }
        if p.y() < self.min.y() {
            self.min = self.min.with_y(p.y());
        }
        if p.z() < self.min.z() {
            self.min = self.min.with_z(p.z());
        }
        if p.x() > self.max.x() {
            self.max = self.max.with_x(p.x());
        }
        if p.y() > self.max.y() {
            self.max = self.max.with_y(p.y());
        }
        if p.z() > self.max.z() {
            self.max = self.max.with_z(p.z());
        }
    }

    pub fn add_box(&mut self, b: BoundingBox) {
        self.add_point(b.min);
        self.add_point(b.max);
    }

    pub fn contains_point(&self, p: Point) -> bool {
        (self.min.x() <= p.x() && p.x() <= self.max.x())
            && (self.min.y() <= p.y() && p.y() <= self.max.y())
            && (self.min.z() <= p.z() && p.z() <= self.max.z())
    }

    pub fn contains_box(&self, b: BoundingBox) -> bool {
        self.contains_point(b.min) && self.contains_point(b.max)
    }

    pub fn intersects(&self, ray: Ray) -> bool {
        let (xtmin, xtmax) = Cube::check_axis(
            ray.origin.x(),
            ray.direction.x(),
            self.min.x(),
            self.max.x(),
        );
        let (ytmin, ytmax) = Cube::check_axis(
            ray.origin.y(),
            ray.direction.y(),
            self.min.y(),
            self.max.y(),
        );
        let (ztmin, ztmax) = Cube::check_axis(
            ray.origin.z(),
            ray.direction.z(),
            self.min.z(),
            self.max.z(),
        );

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

        tmin <= tmax
    }

    pub fn split(&self) -> (BoundingBox, BoundingBox) {
        let dx: OrderedFloat<f32> = (self.max.x() - self.min.x()).abs().into();
        let dy: OrderedFloat<f32> = (self.max.y() - self.min.y()).abs().into();
        let dz: OrderedFloat<f32> = (self.max.z() - self.min.z()).abs().into();

        let greatest = cmp::max(cmp::max(dx, dy), dz);

        let [mut x0, mut y0, mut z0]: [f32; 3] = self.min.into();
        let [mut x1, mut y1, mut z1]: [f32; 3] = self.max.into();

        if greatest == dx {
            x1 = x0 + f32::from(dx) / 2.;
            x0 += f32::from(dx) / 2.;
        } else if greatest == dy {
            y1 = y0 + f32::from(dy) / 2.;
            y0 += f32::from(dy) / 2.;
        } else {
            z1 = z0 + f32::from(dz) / 2.;
            z0 += f32::from(dz) / 2.;
        }

        let mid_min = Point::new(x0, y0, z0);
        let mid_max = Point::new(x1, y1, z1);

        let left = BoundingBox::default().with_min(self.min).with_max(mid_max);
        let right = BoundingBox::default().with_min(mid_min).with_max(self.max);

        (left, right)
    }
}

#[cfg(test)]
mod tests {
    use std::f32::consts::PI;

    use crate::{transformations::Transform, tuples::Vector};

    use super::*;

    #[test]
    fn bounding_box_basics() {
        let mut bb = BoundingBox::default();
        bb.add_point(Point::new(-5., 2., 0.));
        bb.add_point(Point::new(7., 0., -3.));
        assert!(bb.min == Point::new(-5., 0., -3.));
        assert!(bb.max == Point::new(7., 2., 0.));

        let mut bb1 = BoundingBox::default()
            .with_min(Point::new(-5., -2., 0.))
            .with_max(Point::new(7., 4., 4.));
        let bb2 = BoundingBox::default()
            .with_min(Point::new(8., -7., -2.))
            .with_max(Point::new(14., 2., 8.));
        bb1.add_box(bb2);
        assert!(bb1.min == Point::new(-5., -7., -2.));
        assert!(bb1.max == Point::new(14., 4., 8.));

        let bb = BoundingBox::default()
            .with_min(Point::new(5., -2., 0.))
            .with_max(Point::new(11., 4., 7.));
        let cases = vec![
            (Point::new(5., -2., 0.), true),
            (Point::new(11., 4., 7.), true),
            (Point::new(8., 1., 3.), true),
            (Point::new(3., 0., 3.), false),
            (Point::new(8., -4., 3.), false),
            (Point::new(8., 1., -1.), false),
            (Point::new(13., 1., 3.), false),
            (Point::new(8., 5., 3.), false),
            (Point::new(8., 1., 8.), false),
        ];
        for (p, res) in cases {
            assert!(bb.contains_point(p) == res);
        }

        let bb1 = BoundingBox::default()
            .with_min(Point::new(5., -2., 0.))
            .with_max(Point::new(11., 4., 7.));
        let cases = vec![
            (Point::new(5., -2., 0.), Point::new(11., 4., 7.), true),
            (Point::new(6., -1., 1.), Point::new(10., 3., 6.), true),
            (Point::new(4., -3., -1.), Point::new(10., 3., 6.), false),
            (Point::new(6., -1., 1.), Point::new(12., 5., 8.), false),
        ];
        for (min, max, res) in cases {
            let bb2 = BoundingBox::default().with_min(min).with_max(max);
            assert!(bb1.contains_box(bb2) == res);
        }

        let bb = BoundingBox::default()
            .with_min(Point::new(-1., -1., -1.))
            .with_max(Point::new(1., 1., 1.))
            .transform(Transform::default().rotation_y(PI / 4.).rotation_x(PI / 4.));
        assert!(bb.min.equal_approx(Point::new(-1.4142, -1.7071, -1.7071)));
        assert!(bb.max.equal_approx(Point::new(1.4142, 1.7071, 1.7071)));
    }

    #[test]
    fn bounding_box_intersect() {
        let bb = BoundingBox::default()
            .with_min(Point::new(-1., -1., -1.))
            .with_max(Point::new(1., 1., 1.));
        let cases = vec![
            (Point::new(5., 0.5, 0.), Vector::new(-1., 0., 0.), true),
            (Point::new(-5., 0.5, 0.), Vector::new(1., 0., 0.), true),
            (Point::new(0.5, 5., 0.), Vector::new(0., -1., 0.), true),
            (Point::new(0.5, -5., 0.), Vector::new(0., 1., 0.), true),
            (Point::new(0.5, 0., 5.), Vector::new(0., 0., -1.), true),
            (Point::new(0.5, 0., -5.), Vector::new(0., 0., 1.), true),
            (Point::new(0., 0.5, 0.), Vector::new(0., 0., 1.), true),
            (Point::new(-2., 0., 0.), Vector::new(2., 4., 6.), false),
            (Point::new(0., -2., 0.), Vector::new(6., 2., 4.), false),
            (Point::new(0., 0., -2.), Vector::new(4., 6., 2.), false),
            (Point::new(2., 0., 2.), Vector::new(0., 0., -1.), false),
            (Point::new(0., 2., 2.), Vector::new(0., -1., 0.), false),
            (Point::new(2., 2., 0.), Vector::new(-1., 0., 0.), false),
        ];
        for (origin, direction, res) in cases {
            let direction = direction.normalize();
            let ray = Ray::new(origin, direction);
            assert!(bb.intersects(ray) == res);
        }

        let bb = BoundingBox::default()
            .with_min(Point::new(5., -2., 0.))
            .with_max(Point::new(11., 4., 7.));
        let cases = vec![
            (Point::new(15., 1., 2.), Vector::new(-1., 0., 0.), true),
            (Point::new(-5., -1., 4.), Vector::new(1., 0., 0.), true),
            (Point::new(7., 6., 5.), Vector::new(0., -1., 0.), true),
            (Point::new(9., -5., 6.), Vector::new(0., 1., 0.), true),
            (Point::new(8., 2., 1.2), Vector::new(0., 0., -1.), true),
            (Point::new(6., 0., -5.), Vector::new(0., 0., 1.), true),
            (Point::new(8., 1., 3.5), Vector::new(0., 0., 1.), true),
            (Point::new(9., -1., -8.), Vector::new(2., 4., 6.), false),
            (Point::new(8., 3., -4.), Vector::new(6., 2., 4.), false),
            (Point::new(9., -1., -2.), Vector::new(4., 6., 2.), false),
            (Point::new(4., 0., 9.), Vector::new(0., 0., -1.), false),
            (Point::new(8., 6., -1.), Vector::new(0., -1., 0.), false),
            (Point::new(12., 5., 4.), Vector::new(-1., 0., 0.), false),
        ];
        for (origin, direction, res) in cases {
            let direction = direction.normalize();
            let ray = Ray::new(origin, direction);
            assert!(bb.intersects(ray) == res);
        }
    }

    #[test]
    fn bounding_box_splits() {
        let bb = BoundingBox::default()
            .with_min(Point::new(-1., -4., -5.))
            .with_max(Point::new(9., 6., 5.));
        let (left, right) = bb.split();
        assert!(left.min == Point::new(-1., -4., -5.));
        assert!(left.max == Point::new(4., 6., 5.));
        assert!(right.min == Point::new(4., -4., -5.));
        assert!(right.max == Point::new(9., 6., 5.));

        let bb = BoundingBox::default()
            .with_min(Point::new(-1., -2., -3.))
            .with_max(Point::new(9., 5.5, 3.));
        let (left, right) = bb.split();
        assert!(left.min == Point::new(-1., -2., -3.));
        assert!(left.max == Point::new(4., 5.5, 3.));
        assert!(right.min == Point::new(4., -2., -3.));
        assert!(right.max == Point::new(9., 5.5, 3.));

        let bb = BoundingBox::default()
            .with_min(Point::new(-1., -2., -3.))
            .with_max(Point::new(5., 8., 3.));
        let (left, right) = bb.split();
        assert!(left.min == Point::new(-1., -2., -3.));
        assert!(left.max == Point::new(5., 3., 3.));
        assert!(right.min == Point::new(-1., 3., -3.));
        assert!(right.max == Point::new(5., 8., 3.));

        let bb = BoundingBox::default()
            .with_min(Point::new(-1., -2., -3.))
            .with_max(Point::new(5., 3., 7.));
        let (left, right) = bb.split();
        assert!(left.min == Point::new(-1., -2., -3.));
        assert!(left.max == Point::new(5., 3., 2.));
        assert!(right.min == Point::new(-1., -2., 2.));
        assert!(right.max == Point::new(5., 3., 7.));
    }
}
