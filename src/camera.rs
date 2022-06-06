use crate::canvas::Canvas;
use crate::matrices::Matrix;
use crate::rays::Ray;
use crate::tuples::Point;
use crate::world::World;

#[derive(Debug, Clone)]
pub struct Camera {
    pub hsize: u32,
    pub vsize: u32,
    pub field_of_view: f32,
    pub transform: Matrix<4, 4>,
    pub pixel_size: f32,
    pub half_width: f32,
    pub half_height: f32,
}

impl Camera {
    pub fn new(hsize: u32, vsize: u32, field_of_view: f32) -> Self {
        let half_view = (field_of_view / 2.).tan();
        let aspect = hsize as f32 / vsize as f32;
        let (half_width, half_height) = if aspect >= 1. {
            (half_view, half_view / aspect)
        } else {
            (half_view * aspect, half_view)
        };
        let pixel_size = (half_width * 2.) / hsize as f32;

        Self {
            hsize,
            vsize,
            field_of_view,
            transform: Matrix::identity(),
            pixel_size,
            half_width,
            half_height,
        }
    }

    pub fn ray_for_pixel(&self, px: u32, py: u32) -> Ray {
        // The offset from the edge of the canvas to the pixel's center
        let xoffset = (px as f32 + 0.5) * self.pixel_size;
        let yoffset = (py as f32 + 0.5) * self.pixel_size;

        // The untransformed coordinates of the pixel in world space
        // (remember that the camera looks toward -z, so +x is to the *left*)
        let world_x = self.half_width - xoffset;
        let world_y = self.half_height - yoffset;

        // Using the camera matrix, transform the canvas point and the origin,
        // and then compute the ray's direction vector.
        // (remember that the canvas is at z=-1)
        let inverse = self.transform.inverse().unwrap();
        let pixel = inverse * Point::new(world_x, world_y, -1.);
        let origin = inverse * Point::new(0., 0., 0.);
        let direction = (pixel - origin).normalize();

        Ray::new(origin, direction)
    }

    pub fn render(&self, world: &World) -> Canvas {
        use rayon::prelude::*;
        let mut canvas = Canvas::new(self.hsize, self.vsize);
        (0..self.hsize)
            .flat_map(move |px| (0..self.vsize).map(move |py| (px, py)))
            .collect::<Vec<_>>()
            .into_par_iter()
            .map(|(px, py)| {
                let ray = self.ray_for_pixel(px, py);
                let color = world.color_at(ray);
                (px, py, color)
            })
            .collect::<Vec<_>>()
            .into_iter()
            .for_each(|(px, py, color)| {
                canvas.write_pixel(px, py, color);
            });
        canvas
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transformations::Transform;
    use crate::tuples::Vector;
    use std::f32::consts::PI;

    #[test]
    fn camera_basics() {
        let c = Camera::new(200, 125, PI / 2.);
        assert!(c.pixel_size == 0.01);
        let c = Camera::new(125, 200, PI / 2.);
        assert!(c.pixel_size == 0.01);
    }

    #[test]
    fn camera_rays() {
        let c = Camera::new(201, 101, PI / 2.);
        let r = c.ray_for_pixel(100, 50);
        assert!(r.origin == Point::new(0., 0., 0.));
        assert!(r.direction.equal_approx(Vector::new(0., 0., -1.)));

        let c = Camera::new(201, 101, PI / 2.);
        let r = c.ray_for_pixel(0, 0);
        assert!(r.origin == Point::new(0., 0., 0.));
        assert!(r
            .direction
            .equal_approx(Vector::new(0.66519, 0.33259, -0.66851)));

        let mut c = Camera::new(201, 101, PI / 2.);
        c.transform = Transform::new()
            .translation(0., -2., 5.)
            .roation_y(PI / 4.)
            .into();
        let r = c.ray_for_pixel(100, 50);
        assert!(r.origin.equal_approx(Point::new(0., 2., -5.)));
        assert!(r
            .direction
            .equal_approx(Vector::new(f32::sqrt(2.) / 2., 0., -f32::sqrt(2.) / 2.)));
    }
}
