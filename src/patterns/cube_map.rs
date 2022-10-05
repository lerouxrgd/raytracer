use std::cmp;

use ordered_float::OrderedFloat;

use crate::matrices::Matrix;
use crate::tuples::{Color, Point};

use super::{Pattern, UvPattern};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CubeMap {
    pub transform: Matrix<4, 4>,
    pub left: UvPattern,
    pub front: UvPattern,
    pub right: UvPattern,
    pub back: UvPattern,
    pub up: UvPattern,
    pub down: UvPattern,
}

impl CubeMap {
    pub fn new(
        left: UvPattern,
        front: UvPattern,
        right: UvPattern,
        back: UvPattern,
        up: UvPattern,
        down: UvPattern,
    ) -> Self {
        Self {
            transform: Matrix::identity(),
            left,
            front,
            right,
            back,
            up,
            down,
        }
    }

    pub fn with_transform<T: Into<Matrix<4, 4>>>(mut self, transform: T) -> Self {
        self.transform = transform.into();
        self
    }

    pub fn pattern_at(&self, pattern_point: Point) -> Color {
        let face = CubeFace::from(pattern_point);
        let (u, v) = face.uv(pattern_point);
        match face {
            CubeFace::Left => self.left.uv_pattern_at(u, v),
            CubeFace::Front => self.front.uv_pattern_at(u, v),
            CubeFace::Right => self.right.uv_pattern_at(u, v),
            CubeFace::Back => self.back.uv_pattern_at(u, v),
            CubeFace::Up => self.up.uv_pattern_at(u, v),
            CubeFace::Down => self.down.uv_pattern_at(u, v),
        }
    }
}

impl From<CubeMap> for Pattern {
    fn from(cm: CubeMap) -> Self {
        Pattern::CubeMap(cm)
    }
}

impl From<CubeMap> for Option<Pattern> {
    fn from(cm: CubeMap) -> Self {
        Some(Pattern::CubeMap(cm))
    }
}

enum CubeFace {
    Left,
    Right,
    Front,
    Back,
    Up,
    Down,
}

impl From<Point> for CubeFace {
    fn from(p: Point) -> Self {
        let abs_x = OrderedFloat(p.x().abs());
        let abs_y = OrderedFloat(p.y().abs());
        let abs_z = OrderedFloat(p.z().abs());
        let coord: f32 = cmp::max(cmp::max(abs_x, abs_y), abs_z).into();
        if coord == p.x() {
            Self::Right
        } else if coord == -p.x() {
            Self::Left
        } else if coord == p.y() {
            Self::Up
        } else if coord == -p.y() {
            Self::Down
        } else if coord == p.z() {
            Self::Front
        } else {
            Self::Back
        }
    }
}

impl CubeFace {
    pub fn uv(&self, p: Point) -> (f32, f32) {
        match self {
            Self::Front => {
                let u = ((p.x() + 1.) % 2.) / 2.;
                let v = ((p.y() + 1.) % 2.) / 2.;
                (u, v)
            }
            Self::Back => {
                let u = ((1. - p.x()) % 2.) / 2.;
                let v = ((p.y() + 1.) % 2.) / 2.;
                (u, v)
            }
            Self::Up => {
                let u = ((p.x() + 1.) % 2.) / 2.;
                let v = ((1. - p.z()) % 2.) / 2.;
                (u, v)
            }
            Self::Down => {
                let u = ((p.x() + 1.) % 2.) / 2.;
                let v = ((p.z() + 1.) % 2.) / 2.;
                (u, v)
            }
            Self::Left => {
                let u = ((p.z() + 1.) % 2.) / 2.;
                let v = ((p.y() + 1.) % 2.) / 2.;
                (u, v)
            }
            Self::Right => {
                let u = ((1. - p.z()) % 2.) / 2.;
                let v = ((p.y() + 1.) % 2.) / 2.;
                (u, v)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uv_cube_faces() {
        let cases = vec![
            (Point::new(-0.5, 0.5, 1.), 0.25, 0.75),
            (Point::new(0.5, -0.5, 1.), 0.75, 0.25),
        ];
        for (p, u, v) in cases {
            assert!(CubeFace::Front.uv(p) == (u, v));
        }

        let cases = vec![
            (Point::new(0.5, 0.5, -1.), 0.25, 0.75),
            (Point::new(-0.5, -0.5, -1.), 0.75, 0.25),
        ];
        for (p, u, v) in cases {
            assert!(CubeFace::Back.uv(p) == (u, v));
        }

        let cases = vec![
            (Point::new(-1., 0.5, -0.5), 0.25, 0.75),
            (Point::new(-1., -0.5, 0.5), 0.75, 0.25),
        ];
        for (p, u, v) in cases {
            assert!(CubeFace::Left.uv(p) == (u, v));
        }

        let cases = vec![
            (Point::new(1., 0.5, 0.5), 0.25, 0.75),
            (Point::new(1., -0.5, -0.5), 0.75, 0.25),
        ];
        for (p, u, v) in cases {
            assert!(CubeFace::Right.uv(p) == (u, v));
        }

        let cases = vec![
            (Point::new(-0.5, 1., -0.5), 0.25, 0.75),
            (Point::new(0.5, 1., 0.5), 0.75, 0.25),
        ];
        for (p, u, v) in cases {
            assert!(CubeFace::Up.uv(p) == (u, v));
        }

        let cases = vec![
            (Point::new(-0.5, -1., 0.5), 0.25, 0.75),
            (Point::new(0.5, -1., -0.5), 0.75, 0.25),
        ];
        for (p, u, v) in cases {
            assert!(CubeFace::Down.uv(p) == (u, v));
        }
    }
}
