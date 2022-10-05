use crate::matrices::Matrix;
use crate::tuples::{Color, Point};

use super::Pattern;

/// Stripes based on x value
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Striped {
    pub a: Color,
    pub b: Color,
    pub transform: Matrix<4, 4>,
}

impl Striped {
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
        if pattern_point.x().floor() % 2. == 0. {
            self.a
        } else {
            self.b
        }
    }
}

impl From<Striped> for Pattern {
    fn from(s: Striped) -> Self {
        Pattern::Striped(s)
    }
}

impl From<Striped> for Option<Pattern> {
    fn from(s: Striped) -> Self {
        Some(Pattern::Striped(s))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shapes::Sphere;
    use crate::transformations::{scaling, translation};

    #[test]
    fn striped_basics() {
        let p = Striped::new(Color::white(), Color::black());
        assert!(p.pattern_at(Point::new(0., 0., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(0., 1., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(0., 2., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(0., 0., 1.)) == Color::white());
        assert!(p.pattern_at(Point::new(0., 0., 2.)) == Color::white());
        assert!(p.pattern_at(Point::new(0.9, 0., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(1., 0., 0.)) == Color::black());
        assert!(p.pattern_at(Point::new(-0.1, 0., 0.)) == Color::black());
        assert!(p.pattern_at(Point::new(-1., 0., 0.)) == Color::black());
        assert!(p.pattern_at(Point::new(-1.1, 0., 0.)) == Color::white());
    }

    #[test]
    fn striped_transform() {
        let shape = Sphere::default().with_transform(scaling(2., 2., 2.));
        let p: Pattern = Striped::new(Color::white(), Color::black()).into();
        let c = p.pattern_at_shape(shape.into(), Point::new(1.5, 0., 0.));
        assert!(c == Color::white());

        let shape = Sphere::default();
        let p: Pattern = Striped::new(Color::white(), Color::black())
            .with_transform(scaling(2., 2., 2.))
            .into();
        let c = p.pattern_at_shape(shape.into(), Point::new(1.5, 0., 0.));
        assert!(c == Color::white());

        let shape = Sphere::default().with_transform(scaling(2., 2., 2.));
        let p: Pattern = Striped::new(Color::white(), Color::black())
            .with_transform(translation(0.5, 0., 0.))
            .into();
        let c = p.pattern_at_shape(shape.into(), Point::new(2.5, 0., 0.));
        assert!(c == Color::white());
    }
}
