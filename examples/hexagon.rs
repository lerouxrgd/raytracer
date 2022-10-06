use std::f32::consts::PI;

use raytracer::camera::Camera;
use raytracer::groups::hexagon;
use raytracer::lights::PointLight;
use raytracer::transformations::*;
use raytracer::tuples::{Color, Point, Vector};
use raytracer::world::World;

fn main() {
    let mut h = hexagon(
        Transform::default()
            .rotation_x(PI / 3.0)
            .translation(0.0, 0.75, 0.0),
    );
    h.cache_bounds();

    let world = World {
        shapes: vec![],
        csgs: vec![],
        groups: vec![h],
        light: PointLight::new(Point::new(-10., 10., -10.), Color::white()).into(),
    };

    let mut camera = Camera::new(1280, 720, PI / 3.);
    camera.transform = view_transform(
        Point::new(0., 1.5, -5.),
        Point::new(0., 1., 0.),
        Vector::new(0., 1., 0.),
    );

    let canvas = camera.render(&world);
    canvas.to_ppm(std::io::stdout()); // pipe this to a .ppm file
}
