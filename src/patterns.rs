use crate::matrices::Matrix;
use crate::shapes::Shape;
use crate::tuples::{Color, Point};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Pattern {
    Striped(Striped),
    Gradient(Gradient),
    Ring(Ring),
    Checker(Checker),
    XyzRgb(XyzRgb),
}

impl Pattern {
    pub fn transform(&self) -> Matrix<4, 4> {
        match self {
            &Self::Striped(Striped { transform, .. })
            | &Self::Gradient(Gradient { transform, .. })
            | &Self::Ring(Ring { transform, .. })
            | &Self::Checker(Checker { transform, .. })
            | &Self::XyzRgb(XyzRgb { transform, .. }) => transform,
        }
    }

    pub fn transform_mut(&mut self) -> &mut Matrix<4, 4> {
        match self {
            Self::Striped(Striped {
                ref mut transform, ..
            })
            | Self::Gradient(Gradient {
                ref mut transform, ..
            })
            | Self::Ring(Ring {
                ref mut transform, ..
            })
            | Self::Checker(Checker {
                ref mut transform, ..
            })
            | Self::XyzRgb(XyzRgb {
                ref mut transform, ..
            }) => transform,
        }
    }

    pub fn pattern_at_shape(&self, shape: Shape, world_point: Point) -> Color {
        let object_point = shape.transform().inverse().unwrap() * world_point;
        let pattern_point = self.transform().inverse().unwrap() * object_point;
        match self {
            Self::Striped(s) => s.pattern_at(pattern_point),
            Self::Gradient(g) => g.pattern_at(pattern_point),
            Self::Ring(r) => r.pattern_at(pattern_point),
            Self::Checker(c) => c.pattern_at(pattern_point),
            Self::XyzRgb(xyz) => xyz.pattern_at(pattern_point),
        }
    }
}

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

/// Coordinates are the color components
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct XyzRgb {
    pub transform: Matrix<4, 4>,
}

impl XyzRgb {
    pub fn new() -> Self {
        Self {
            transform: Matrix::identity(),
        }
    }

    pub fn pattern_at(&self, pattern_point: Point) -> Color {
        Color::new(pattern_point.x(), pattern_point.y(), pattern_point.z())
    }
}

impl From<XyzRgb> for Pattern {
    fn from(xyz: XyzRgb) -> Self {
        Pattern::XyzRgb(xyz)
    }
}

impl From<XyzRgb> for Option<Pattern> {
    fn from(xyz: XyzRgb) -> Self {
        Some(Pattern::XyzRgb(xyz))
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
        let mut shape = Sphere::new();
        shape.transform = scaling(2., 2., 2.);
        let p: Pattern = Striped::new(Color::white(), Color::black()).into();
        let c = p.pattern_at_shape(shape.into(), Point::new(1.5, 0., 0.));
        assert!(c == Color::white());

        let shape = Sphere::new();
        let mut p: Pattern = Striped::new(Color::white(), Color::black()).into();
        *p.transform_mut() = scaling(2., 2., 2.);
        let c = p.pattern_at_shape(shape.into(), Point::new(1.5, 0., 0.));
        assert!(c == Color::white());

        let mut shape = Sphere::new();
        shape.transform = scaling(2., 2., 2.);
        let mut p: Pattern = Striped::new(Color::white(), Color::black()).into();
        *p.transform_mut() = translation(0.5, 0., 0.);
        let c = p.pattern_at_shape(shape.into(), Point::new(2.5, 0., 0.));
        assert!(c == Color::white());
    }

    #[test]
    fn gradient_basics() {
        let p = Gradient::new(Color::white(), Color::black());
        assert!(p.pattern_at(Point::new(0., 0., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(0.25, 0., 0.)) == Color::new(0.75, 0.75, 0.75));
        assert!(p.pattern_at(Point::new(0.5, 0., 0.)) == Color::new(0.5, 0.5, 0.5));
        assert!(p.pattern_at(Point::new(0.75, 0., 0.)) == Color::new(0.25, 0.25, 0.25));
    }

    #[test]
    fn ring_basics() {
        let p = Ring::new(Color::white(), Color::black());
        assert!(p.pattern_at(Point::new(0., 0., 0.)) == Color::white());
        assert!(p.pattern_at(Point::new(1., 0., 0.)) == Color::black());
        assert!(p.pattern_at(Point::new(0., 0., 1.)) == Color::black());
        assert!(p.pattern_at(Point::new(0.708, 0., 0.708)) == Color::black());
    }

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
