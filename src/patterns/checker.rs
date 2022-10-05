use crate::matrices::Matrix;
use crate::tuples::{Color, Point};

use super::Pattern;

/// Repeating pattern of squares (adjacent squares have different colors)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Checker {
    pub a: Color,
    pub b: Color,
    pub transform: Matrix<4, 4>,
}

impl Checker {
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
        if (pattern_point.x().floor() + pattern_point.y().floor() + pattern_point.z().floor()) % 2.
            == 0.
        {
            self.a
        } else {
            self.b
        }
    }
}

impl From<Checker> for Pattern {
    fn from(s: Checker) -> Self {
        Pattern::Checker(s)
    }
}

impl From<Checker> for Option<Pattern> {
    fn from(s: Checker) -> Self {
        Some(Pattern::Checker(s))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checker_basics() {
        let p = Checker::new(Color::white(), Color::black());
        assert!(p.pattern_at(Point::new(0., 0., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(0.99, 0., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(1.01, 0., 0.)) == Color::black());
        assert!(p.pattern_at(Point::new(0., 0., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(0., 0.99, 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(0., 1.01, 0.)) == Color::black());
        assert!(p.pattern_at(Point::new(0., 0., 0.99)) == Color::white());
        assert!(p.pattern_at(Point::new(0., 0., 1.01)) == Color::black());
    }
}
