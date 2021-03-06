use std::f32::consts::PI;

use raytracer::camera::Camera;
use raytracer::lights::PointLight;
use raytracer::patterns::*;
use raytracer::shapes::{Plane, Sphere};
use raytracer::transformations::*;
use raytracer::tuples::{Color, Point, Vector};
use raytracer::world::World;

// TODO: consider "caching" inverse by making it a regular field

fn main() {
    let mut floor = Plane::new();
    let mut pattern = Checker::new(Color::white(), Color::black());
    pattern.transform = Transform::new()
        .rotation_y(PI / 4.)
        .scaling(0.4, 0.4, 0.4)
        .into();
    floor.material.pattern = pattern.into();
    floor.material.color = Color::new(1., 0.9, 0.9);
    floor.material.specular = 0.;
    floor.material.reflective = 0.3;

    let mut backdrop = Plane::new();
    backdrop.transform = Transform::new()
        .rotation_x(PI / 2.)
        .translation(0., 0., 5.)
        .into();
    backdrop.material.color = Color::new(1., 0.9, 0.9);
    backdrop.material.specular = 0.;

    let mut middle = Sphere::new();
    middle.transform = translation(-0.5, 1., 0.5);
    middle.material.color = Color::new(0.1, 0.4, 0.9);
    middle.material.diffuse = 0.7;
    middle.material.specular = 0.3;
    middle.material.reflective = 0.8;

    let mut right = Sphere::new();
    right.transform = Transform::new()
        .scaling(0.5, 0.5, 0.5)
        .translation(1.5, 0.5, -0.5)
        .into();
    right.material.color = Color::new(0.5, 1., 0.1);
    right.material.diffuse = 0.7;
    right.material.specular = 0.3;

    let mut left = Sphere::new();
    left.transform = Transform::new()
        .scaling(0.33, 0.33, 0.33)
        .translation(-1.5, 0.33, -0.75)
        .into();
    left.material.color = Color::new(1., 0.8, 0.1);
    left.material.diffuse = 0.7;
    left.material.specular = 0.3;

    let objects = vec![
        floor.into(),
        backdrop.into(),
        middle.into(),
        right.into(),
        left.into(),
    ];

    let groups = vec![
        // //
        // raytracer::groups::hexagon(
        //     Transform::new()
        //         .rotation_x(PI / 3.0)
        //         .translation(0.0, 0.75, 0.0)
        //         .into(),
        // ),
    ];

    let world = World {
        objects,
        groups,
        light: PointLight::new(Point::new(-10., 10., -10.), Color::white()),
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
