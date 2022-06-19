use crate::tuples::{Color, Point};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointLight {
    pub intensity: Color,
    pub position: Point,
}

impl PointLight {
    pub fn new(position: Point, intensity: Color) -> Self {
        Self {
            intensity,
            position,
        }
    }
}
