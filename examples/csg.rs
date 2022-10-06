use std::f32::consts::PI;

use raytracer::camera::Camera;
use raytracer::csg::{Csg, CsgOp};
use raytracer::lights::PointLight;
use raytracer::materials::Material;
use raytracer::shapes::{Cylinder, Shape};
use raytracer::transformations::*;
use raytracer::tuples::{Color, Point, Vector};
use raytracer::world::World;

fn main() {
    let outer = Csg::new(
        CsgOp::Intersect,
        Shape::sphere().with_transform(scaling(1.3, 1.3, 1.3)),
        Shape::cube().with_transform(rotation_y(PI / 3.)),
    );

    let inner = Csg::new(
        CsgOp::Union,
        Csg::new(
            CsgOp::Union,
            Cylinder::default()
                .with_transform(
                    Transform::default()
                        .scaling(0.6, 0.6, 0.6)
                        .rotation_y(PI / 3.),
                )
                .with_material(Material::default().color(Color::new(1., 0., 0.))),
            Cylinder::default()
                .with_transform(
                    Transform::default()
                        .scaling(0.6, 0.6, 0.6)
                        .rotation_x(PI / 2.)
                        .rotation_y(PI / 3.),
                )
                .with_material(Material::default().color(Color::new(0., 0., 1.))),
        ),
        Cylinder::default()
            .with_transform(
                Transform::default()
                    .scaling(0.6, 0.6, 0.6)
                    .rotation_z(PI / 2.)
                    .rotation_y(PI / 3.),
            )
            .with_material(Material::default().color(Color::new(0., 1., 0.))),
    );

    let world = World {
        shapes: vec![],
        groups: vec![],
        csgs: vec![Csg::new(CsgOp::Difference, outer, inner)],
        light: PointLight::new(Point::new(-10., 10., -10.), Color::white()).into(),
        ..Default::default()
    };

    let mut camera = Camera::new(1280, 720, PI / 3.);
    camera.transform = view_transform(
        Point::new(0., 1.5, -5.),
        Point::new(0., 0., 0.),
        Vector::new(0., 1., 0.),
    );

    let canvas = camera.render(&world);
    canvas.to_ppm(std::io::stdout()); // pipe this to a .ppm file
}
