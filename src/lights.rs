#[cfg(test)]
use std::{cell::RefCell, iter::Cycle, vec::IntoIter};

use crate::tuples::{Color, Point, Vector};
use crate::world::World;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Light {
    PointLight(PointLight),
    AreaLight(AreaLight),
}

impl Light {
    pub fn intensity(&self) -> Color {
        match self {
            Self::PointLight(l) => l.intensity,
            Self::AreaLight(l) => l.intensity,
        }
    }

    pub fn position(&self) -> Point {
        match self {
            Self::PointLight(l) => l.position,
            Self::AreaLight(l) => l.position,
        }
    }

    pub fn intensity_at(&self, p: Point, w: &World) -> f32 {
        match self {
            Self::PointLight(l) => l.intensity_at(p, w),
            Self::AreaLight(l) => l.intensity_at(p, w),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointLight {
    intensity: Color,
    position: Point,
}

impl PointLight {
    pub fn new(position: Point, intensity: Color) -> Self {
        Self {
            intensity,
            position,
        }
    }

    pub fn intensity_at(&self, p: Point, w: &World) -> f32 {
        if w.is_shadowed(self.position, p) {
            0.
        } else {
            1.
        }
    }
}

impl From<PointLight> for Light {
    fn from(l: PointLight) -> Self {
        Self::PointLight(l)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AreaLight {
    position: Point,
    corner: Point,
    uvec: Vector,
    usteps: usize,
    vvec: Vector,
    vsteps: usize,
    intensity: Color,
}

impl AreaLight {
    #[cfg(test)]
    thread_local! {
        static JITTER_BY: RefCell<Cycle<IntoIter<f32>>> =
            RefCell::new(vec![].into_iter().cycle());
    }

    pub fn new(
        corner: Point,
        full_uvec: Vector,
        usteps: usize,
        full_vvec: Vector,
        vsteps: usize,
        intensity: Color,
        #[cfg(test)] jitter_by: Vec<f32>,
    ) -> Self {
        #[cfg(test)]
        Self::JITTER_BY.with(|jitter| *jitter.borrow_mut() = jitter_by.into_iter().cycle());
        Self {
            position: corner + (full_uvec + full_vvec) / 2.,
            corner,
            uvec: full_uvec / usteps as f32,
            usteps,
            vvec: full_vvec / vsteps as f32,
            vsteps,
            intensity,
        }
    }

    pub fn point_on_light(&self, u: usize, v: usize) -> Point {
        #[cfg(test)]
        {
            let (jitter_u, jitter_v) = Self::JITTER_BY.with(|jitter| {
                let mut jitter = jitter.borrow_mut();
                (jitter.next().unwrap(), jitter.next().unwrap())
            });
            self.corner + (u as f32 + jitter_u) * self.uvec + (v as f32 + jitter_v) * self.vvec
        }
        #[cfg(not(test))]
        {
            use rand::prelude::*;
            let mut rng = thread_rng();
            let (jitter_u, jitter_v) = (rng.gen::<f32>(), rng.gen::<f32>());
            self.corner + (u as f32 + jitter_u) * self.uvec + (v as f32 + jitter_v) * self.vvec
        }
    }

    pub fn intensity_at(&self, p: Point, w: &World) -> f32 {
        let mut total = 0.;
        for u in 0..self.usteps {
            for v in 0..self.vsteps {
                let light_pos = self.point_on_light(u, v);
                if !w.is_shadowed(light_pos, p) {
                    total += 1.;
                }
            }
        }
        total / self.nb_samples()
    }

    #[inline(always)]
    pub fn usteps(&self) -> usize {
        self.usteps
    }

    #[inline(always)]
    pub fn vsteps(&self) -> usize {
        self.vsteps
    }

    #[inline(always)]
    pub fn nb_samples(&self) -> f32 {
        (self.usteps * self.vsteps) as f32
    }

    #[inline(always)]
    pub fn intensity(&self) -> Color {
        self.intensity
    }
}

impl From<AreaLight> for Light {
    fn from(l: AreaLight) -> Self {
        Self::AreaLight(l)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_light_intensity() {
        let cases = vec![
            (Point::new(0., 1.0001, 0.), 1.),
            (Point::new(-1.0001, 0., 0.), 1.),
            (Point::new(0., 0., -1.0001), 1.),
            (Point::new(0., 0., 1.0001), 0.),
            (Point::new(1.0001, 0., 0.), 0.),
            (Point::new(0., -1.0001, 0.), 0.),
            (Point::new(0., 0., 0.), 0.),
        ];
        let w = World::default();
        for (point, res) in cases {
            assert!(w.lights[0].intensity_at(point, &w) == res);
        }
    }

    #[test]
    fn area_light_points() {
        let cases = vec![
            (0, 0, Point::new(0.15, 0., 0.35)),
            (1, 0, Point::new(0.65, 0., 0.35)),
            (0, 1, Point::new(0.15, 0., 0.85)),
            (2, 0, Point::new(1.15, 0., 0.35)),
            (3, 1, Point::new(1.65, 0., 0.85)),
        ];
        for (u, v, res) in cases {
            let light = AreaLight::new(
                Point::new(0., 0., 0.),
                Vector::new(2., 0., 0.),
                4,
                Vector::new(0., 0., 1.),
                2,
                Color::white(),
                vec![0.3, 0.7],
            );
            assert!(light.point_on_light(u, v) == res);
        }
    }

    #[test]
    fn area_light_intensity() {
        let cases = vec![
            (Point::new(0., 0., 2.), 0.),
            (Point::new(1., -1., 2.), 0.5),
            // (Point::new(1.5, 0., 2.), 0.75),
            (Point::new(1.25, 1.25, 3.), 0.75),
            (Point::new(0., 0., -2.), 1.),
        ];
        let w = World::default();
        for (point, res) in cases {
            let light = AreaLight::new(
                Point::new(-0.5, -0.5, -5.),
                Vector::new(1., 0., 0.),
                2,
                Vector::new(0., 1., 0.),
                2,
                Color::white(),
                vec![0.7, 0.3, 0.9, 0.1, 0.5],
            );
            assert!(light.intensity_at(point, &w) == res);
        }
    }
}
