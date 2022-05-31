use raytracer::canvas::Canvas;
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
    let color = Color::new(1., 0., 0.);
    let shape = Sphere::new();

    for y in 0..canvas_pixels {
        let world_y = half - pixel_size * y as f32;
        for x in 0..canvas_pixels {
            let world_x = -half + pixel_size * x as f32;
            let position = Point::new(world_x, world_y, wall_z); // on the wall
            let ray_direction = (position - ray_origin).normalize();
            let ray = Ray::new(position, ray_direction);
            let xs = shape.intersect(ray);
            if xs.is_some() {
                canvas.write_pixel(x, y, color);
            }
        }
    }

    canvas.to_ppm(std::io::stdout()); // pipe this to a .ppm file
}
