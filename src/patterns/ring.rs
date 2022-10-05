use crate::matrices::Matrix;
use crate::tuples::{Color, Point};

use super::Pattern;

/// An xz ring
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ring {
    pub a: Color,
    pub b: Color,
    pub transform: Matrix<4, 4>,
}

impl Ring {
    pub fn new(a: Color, b: Color) -> Self {
        Self {
            a,
            b,
            transform: Matrix::identity(),
        }
    }

    pub fn with_transform<T: Into<Matrix<4, 4>>>(mut self, transform: T) -> Self {
        self.transform = transform.into();
        self
    }

    pub fn pattern_at(&self, pattern_point: Point) -> Color {
        if (pattern_point.x().powi(2) + pattern_point.z().powi(2))
            .sqrt()
            .floor()
            % 2.
            == 0.
        {
            self.a
        } else {
            self.b
        }
    }
}

impl From<Ring> for Pattern {
    fn from(s: Ring) -> Self {
        Pattern::Ring(s)
    }
}

impl From<Ring> for Option<Pattern> {
    fn from(s: Ring) -> Self {
        Some(Pattern::Ring(s))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ring_basics() {
        let p = Ring::new(Color::white(), Color::black());
        assert!(p.pattern_at(Point::new(0., 0., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(1., 0., 0.)) == Color::black());
        assert!(p.pattern_at(Point::new(0., 0., 1.)) == Color::black());
        assert!(p.pattern_at(Point::new(0.708, 0., 0.708)) == Color::black());
    }
}
