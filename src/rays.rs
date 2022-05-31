use crate::matrices::Matrix;
use crate::tuples::{Point, Vector};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ray {
    pub origin: Point,
    pub direction: Vector,
}

impl Ray {
    pub fn new(origin: Point, direction: Vector) -> Self {
        Self { origin, direction }
    }

    pub fn position(&self, t: f32) -> Point {
        self.origin + t * self.direction
    }

    pub fn transform(&self, m: Matrix<4, 4>) -> Ray {
        Self {
            origin: m * self.origin,
            direction: m * self.direction,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transformations::*;

    #[test]
    fn ray_basics() {
        let r = Ray::new(Point::new(2., 3., 4.), Vector::new(1., 0., 0.));
        assert!(r.position(0.) == Point::new(2., 3., 4.));
        assert!(r.position(1.) == Point::new(3., 3., 4.));
        assert!(r.position(-1.) == Point::new(1., 3., 4.));
        assert!(r.position(2.5) == Point::new(4.5, 3., 4.));
    }

    #[test]
    fn ray_transform() {
        let r1 = Ray::new(Point::new(1., 2., 3.), Vector::new(0., 1., 0.));
        let m = translation(3., 4., 5.);
        let r2 = r1.transform(m);
        assert!(r2.origin == Point::new(4., 6., 8.));
        assert!(r2.direction == Vector::new(0., 1., 0.));

        let r1 = Ray::new(Point::new(1., 2., 3.), Vector::new(0., 1., 0.));
        let m = scaling(2., 3., 4.);
        let r2 = r1.transform(m);
        assert!(r2.origin == Point::new(2., 6., 12.));
        assert!(r2.direction == Vector::new(0., 3., 0.));
    }
}
