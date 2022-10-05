mod checker;
mod cube_map;
mod gradient;
mod ring;
mod striped;
mod texture_map;
mod uv_pattern;
mod xyz_rgb;

use crate::matrices::Matrix;
use crate::shapes::Shape;
use crate::tuples::{Color, Point};

pub use self::checker::Checker;
pub use self::cube_map::CubeMap;
pub use self::gradient::Gradient;
pub use self::ring::Ring;
pub use self::striped::Striped;
pub use self::texture_map::{TextureMap, UvMapping};
pub use self::uv_pattern::{UvAlignCheck, UvChecker, UvImage, UvPattern};
pub use self::xyz_rgb::XyzRgb;

#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Pattern {
    Striped(Striped),
    Gradient(Gradient),
    Ring(Ring),
    Checker(Checker),
    XyzRgb(XyzRgb),
    TextureMap(TextureMap),
    CubeMap(CubeMap),
}

impl Pattern {
    pub fn with_transform<T: Into<Matrix<4, 4>>>(mut self, transform: T) -> Self {
        self.set_transform(transform.into());
        self
    }

    pub fn get_transform(&self) -> Matrix<4, 4> {
        match self {
            &Self::Striped(Striped { transform, .. })
            | &Self::Gradient(Gradient { transform, .. })
            | &Self::Ring(Ring { transform, .. })
            | &Self::Checker(Checker { transform, .. })
            | &Self::XyzRgb(XyzRgb { transform, .. })
            | &Self::TextureMap(TextureMap { transform, .. })
            | &Self::CubeMap(CubeMap { transform, .. }) => transform,
        }
    }

    pub fn set_transform<T: Into<Matrix<4, 4>>>(&mut self, t: T) {
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
            })
            | Self::TextureMap(TextureMap {
                ref mut transform, ..
            })
            | Self::CubeMap(CubeMap {
                ref mut transform, ..
            }) => *transform = t.into(),
        }
    }

    pub fn pattern_at_shape(&self, shape: Shape, world_point: Point) -> Color {
        let object_point = shape.world_to_object(world_point);
        let pattern_point = self.get_transform().inverse().unwrap() * object_point;
        match self {
            Self::Striped(s) => s.pattern_at(pattern_point),
            Self::Gradient(g) => g.pattern_at(pattern_point),
            Self::Ring(r) => r.pattern_at(pattern_point),
            Self::Checker(c) => c.pattern_at(pattern_point),
            Self::XyzRgb(xyz) => xyz.pattern_at(pattern_point),
            Self::TextureMap(tm) => tm.pattern_at(pattern_point),
            Self::CubeMap(cm) => cm.pattern_at(pattern_point),
        }
    }
}
