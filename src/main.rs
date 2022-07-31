use std::f32::consts::PI;

use raytracer::camera::Camera;
use raytracer::lights::PointLight;
use raytracer::materials::Material;
use raytracer::patterns::*;
use raytracer::shapes::Shape;
use raytracer::transformations::*;
use raytracer::tuples::{Color, Point, Vector};
use raytracer::world::World;

fn main() {
    let floor = Shape::plane()
        .with_transform(
            Transform::default()
                .rotation_y(PI / 4.)
                .scaling(0.4, 0.4, 0.4),
        )
        .with_material(
            Material::default()
                .pattern(Checker::new(Color::white(), Color::black()))
                .color(Color::new(1., 0.9, 0.9))
                .specular(0.)
                .reflective(0.3),
        );

    let backdrop = Shape::plane()
        .with_transform(
            Transform::default()
                .rotation_x(PI / 2.)
                .translation(0., 0., 5.),
        )
        .with_material(
            Material::default()
                .color(Color::new(1., 0.9, 0.9))
                .specular(0.),
        );

    let middle = Shape::sphere()
        .with_transform(translation(-0.5, 1., 0.5))
        .with_material(
            Material::default()
                .color(Color::new(0.1, 0.4, 0.9))
                .diffuse(0.7)
                .specular(0.3)
                .reflective(0.8),
        );

    let right = Shape::sphere()
        .with_transform(
            Transform::default()
                .scaling(0.5, 0.5, 0.5)
                .translation(1.5, 0.5, -0.5),
        )
        .with_material(
            Material::default()
                .color(Color::new(0.5, 1., 0.1))
                .diffuse(0.7)
                .specular(0.3),
        );

    let left = Shape::sphere()
        .with_transform(
            Transform::default()
                .scaling(0.33, 0.33, 0.33)
                .translation(-1.5, 0.33, -0.75),
        )
        .with_material(
            Material::default()
                .color(Color::new(1., 0.8, 0.1))
                .diffuse(0.7)
                .specular(0.3),
        );

    let shapes = vec![floor, backdrop, middle, right, left];

    let groups = vec![
        // //
        // raytracer::groups::hexagon(
        //     Transform::default()
        //         .rotation_x(PI / 3.0)
        //         .translation(0.0, 0.75, 0.0),
        // ),
    ];

    let world = World {
        shapes,
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
