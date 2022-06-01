use crate::tuples::{Color, Point};

#[derive(Debug, Clone, Copy)]
pub struct LightPoint {
    pub intensity: Color,
    pub position: Point,
}

impl LightPoint {
    pub fn new(position: Point, intensity: Color) -> Self {
        Self {
            intensity,
            position,
        }
    }
}
