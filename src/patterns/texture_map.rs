use std::f32::consts::PI;

use serde::{Deserialize, Serialize};

use crate::matrices::Matrix;
use crate::tuples::{Color, Point, Vector};

use super::{Pattern, UvPattern};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TextureMap {
    pub transform: Matrix<4, 4>,
    pub uv_pattern: UvPattern,
    pub uv_mapping: UvMapping,
}

impl TextureMap {
    pub fn new(uv_pattern: UvPattern, uv_mapping: UvMapping) -> Self {
        Self {
            transform: Matrix::identity(),
            uv_pattern,
            uv_mapping,
        }
    }

    pub fn with_transform<T: Into<Matrix<4, 4>>>(mut self, transform: T) -> Self {
        self.transform = transform.into();
        self
    }

    pub fn pattern_at(&self, pattern_point: Point) -> Color {
        let (u, v) = self.uv_mapping.map_point(pattern_point);
        self.uv_pattern.uv_pattern_at(u, v)
    }
}

impl From<TextureMap> for Pattern {
    fn from(tm: TextureMap) -> Self {
        Pattern::TextureMap(tm)
    }
}

impl From<TextureMap> for Option<Pattern> {
    fn from(tm: TextureMap) -> Self {
        Some(Pattern::TextureMap(tm))
    }
}

////////////////////////////////////////////////////////////////////////////////////////

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum UvMapping {
    Spherical,
    Planar,
    Cylindrical,
}

impl UvMapping {
    pub fn map_point(&self, p: Point) -> (f32, f32) {
        match self {
            Self::Spherical => spherical_map(p),
            Self::Planar => planar_map(p),
            Self::Cylindrical => cylindrical_map(p),
        }
    }
}

fn spherical_map(p: Point) -> (f32, f32) {
    // compute the azimuthal angle
    // -π < theta <= π
    // angle increases clockwise as viewed from above,
    // which is opposite of what we want, but we'll fix it later.
    let theta = p.x().atan2(p.z());

    // v is the vector pointing from the sphere's origin (the world origin)
    // to the point, which will also happen to be exactly equal to the sphere's
    // radius.
    let v = Vector::new(p.x(), p.y(), p.z());
    let radius = v.magnitude();

    // compute the polar angle
    // 0 <= phi <= π
    let phi = (p.y() / radius).acos();

    // -0.5 < raw_u <= 0.5
    let raw_u = theta / (2. * PI);

    // 0 <= u < 1
    // here's also where we fix the direction of u. Subtract it from 1,
    // so that it increases counterclockwise as viewed from above.
    let u = 1. - (raw_u + 0.5);

    // we want v to be 0 at the south pole of the sphere,
    // and 1 at the north pole, so we have to "flip it over"
    // by subtracting it from 1.
    let v = 1. - phi / PI;

    (u, v)
}

fn planar_map(p: Point) -> (f32, f32) {
    (p.x().rem_euclid(1.), p.z().rem_euclid(1.))
}

fn cylindrical_map(p: Point) -> (f32, f32) {
    // compute the azimuthal angle, same as with spherical_map()
    let theta = p.x().atan2(p.z());
    let raw_u = theta / (2. * PI);
    let u = 1. - (raw_u + 0.5);
    //let v go from 0 to 1 between whole units of y
    let v = p.y().rem_euclid(1.);
    (u, v)
}

#[cfg(test)]
mod tests {
    use crate::patterns::UvChecker;

    use super::*;

    #[test]
    fn spherical_mapping() {
        let cases = vec![
            (Point::new(0., 0., -1.), 0.0, 0.5),
            (Point::new(1., 0., 0.), 0.25, 0.5),
            (Point::new(0., 0., 1.), 0.5, 0.5),
            (Point::new(-1., 0., 0.), 0.75, 0.5),
            (Point::new(0., 1., 0.), 0.5, 1.0),
            (Point::new(0., -1., 0.), 0.5, 0.0),
            (
                Point::new(f32::sqrt(2.) / 2., f32::sqrt(2.) / 2., 0.),
                0.25,
                0.75,
            ),
        ];
        for (p, u, v) in cases {
            let (res_u, res_v) = spherical_map(p);
            assert!((res_u - u).abs() < 1e-4);
            assert!((res_v - v).abs() < 1e-4);
        }
    }

    #[test]
    fn texture_map_spherical() {
        let cases = vec![
            (Point::new(0.4315, 0.4670, 0.7719), Color::white()),
            (Point::new(-0.9654, 0.2552, -0.0534), Color::black()),
            (Point::new(0.1039, 0.7090, 0.6975), Color::white()),
            (Point::new(-0.4986, -0.7856, -0.3663), Color::black()),
            (Point::new(-0.0317, -0.9395, 0.3411), Color::black()),
            (Point::new(0.4809, -0.7721, 0.4154), Color::black()),
            (Point::new(0.0285, -0.9612, -0.2745), Color::black()),
            (Point::new(-0.5734, -0.2162, -0.7903), Color::white()),
            (Point::new(0.7688, -0.1470, 0.6223), Color::black()),
            (Point::new(-0.7652, 0.2175, 0.6060), Color::black()),
        ];
        let tm = TextureMap::new(
            UvChecker::new(16., 8., Color::black(), Color::white()).into(),
            UvMapping::Spherical,
        );
        for (p, color) in cases {
            assert!(tm.pattern_at(p) == color);
        }
    }

    #[test]
    fn planar_mapping() {
        let cases = vec![
            (Point::new(0.25, 0., 0.5), 0.25, 0.5),
            (Point::new(0.25, 0., -0.25), 0.25, 0.75),
            (Point::new(0.25, 0.5, -0.25), 0.25, 0.75),
            (Point::new(1.25, 0., 0.5), 0.25, 0.5),
            (Point::new(0.25, 0., -1.75), 0.25, 0.25),
            (Point::new(1., 0., -1.), 0.0, 0.0),
            (Point::new(0., 0., 0.), 0.0, 0.0),
        ];
        for (p, u, v) in cases {
            assert!(planar_map(p) == (u, v));
        }
    }

    #[test]
    fn cylindrical_mapping() {
        let cases = vec![
            (Point::new(0., 0., -1.), 0.0, 0.0),
            (Point::new(0., 0.5, -1.), 0.0, 0.5),
            (Point::new(0., 1., -1.), 0.0, 0.0),
            (Point::new(0.70711, 0.5, -0.70711), 0.125, 0.5),
            (Point::new(1., 0.5, 0.), 0.25, 0.5),
            (Point::new(0.70711, 0.5, 0.70711), 0.375, 0.5),
            (Point::new(0., -0.25, 1.), 0.5, 0.75),
            (Point::new(-0.70711, 0.5, 0.70711), 0.625, 0.5),
            (Point::new(-1., 1.25, 0.), 0.75, 0.25),
            (Point::new(-0.70711, 0.5, -0.70711), 0.875, 0.5),
        ];
        for (p, u, v) in cases {
            assert!(cylindrical_map(p) == (u, v));
        }
    }
}
