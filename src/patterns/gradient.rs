use crate::matrices::Matrix;
use crate::tuples::{Color, Point};

use super::Pattern;

/// Linearly blend two colors
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Gradient {
    pub a: Color,
    pub b: Color,
    pub transform: Matrix<4, 4>,
}

impl Gradient {
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
        let distance = self.b - self.a;
        let fraction = pattern_point.x() - pattern_point.x().floor();
        self.a + distance * fraction
    }
}

impl From<Gradient> for Pattern {
    fn from(s: Gradient) -> Self {
        Pattern::Gradient(s)
    }
}

impl From<Gradient> for Option<Pattern> {
    fn from(s: Gradient) -> Self {
        Some(Pattern::Gradient(s))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gradient_basics() {
        let p = Gradient::new(Color::white(), Color::black());
        assert!(p.pattern_at(Point::new(0., 0., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(0.25, 0., 0.)) == Color::new(0.75, 0.75, 0.75));
        assert!(p.pattern_at(Point::new(0.5, 0., 0.)) == Color::new(0.5, 0.5, 0.5));
        assert!(p.pattern_at(Point::new(0.75, 0., 0.)) == Color::new(0.25, 0.25, 0.25));
    }
}
