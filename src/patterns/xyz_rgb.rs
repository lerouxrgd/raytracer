use crate::matrices::Matrix;
use crate::tuples::{Color, Point};

use super::Pattern;

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

    pub fn with_transform<T: Into<Matrix<4, 4>>>(mut self, transform: T) -> Self {
        self.transform = transform.into();
        self
    }

    pub fn pattern_at(&self, pattern_point: Point) -> Color {
        Color::new(pattern_point.x(), pattern_point.y(), pattern_point.z())
    }
}

impl Default for XyzRgb {
    fn default() -> Self {
        Self::new()
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
