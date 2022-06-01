use raytracer::canvas::Canvas;
use raytracer::lights::LightPoint;
use raytracer::materials::lighting;
use raytracer::rays::Ray;
use raytracer::spheres::Sphere;
use raytracer::tuples::{Color, Point};

fn main() {
    let ray_origin = Point::new(0., 0., -5.); // start the ray at z = -5
    let wall_z = 10.;
    let wall_size = 7.; // in order to see it completely
    let canvas_pixels = 100;
    let pixel_size = wall_size / canvas_pixels as f32;
    let half = wall_size / 2.;

    let mut canvas = Canvas::new(canvas_pixels, canvas_pixels);

    let mut shape = Sphere::new();
    shape.material.color = Color::new(1., 0.2, 1.);

    let light_position = Point::new(-10., 10., -10.);
    let light_color = Color::white();
    let light = LightPoint::new(light_position, light_color);

    for y in 0..canvas_pixels {
        let world_y = half - pixel_size * y as f32;
        for x in 0..canvas_pixels {
            let world_x = -half + pixel_size * x as f32;
            let position = Point::new(world_x, world_y, wall_z); // on the wall
            let ray_direction = (position - ray_origin).normalize();
            let ray = Ray::new(position, ray_direction);
            let xs = shape.intersect(ray);
            if let Some(hit) = xs {
                let hit = hit[0];
                let point = ray.position(hit.t());
                let normal = hit.object().normal_at(point);
                let eye = -ray.direction;
                let color = lighting(shape.material, light, point, eye, normal);
                canvas.write_pixel(x, y, color);
            }
        }
    }

    canvas.to_ppm(std::io::stdout()); // pipe this to a .ppm file
}
