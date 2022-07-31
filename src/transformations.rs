use std::ops::Mul;

use crate::matrices::Matrix;
use crate::tuples::{Point, Vector};

pub fn translation(x: f32, y: f32, z: f32) -> Matrix<4, 4> {
    let mut transform = Matrix::identity();
    transform[0][3] = x;
    transform[1][3] = y;
    transform[2][3] = z;
    transform
}

pub fn scaling(x: f32, y: f32, z: f32) -> Matrix<4, 4> {
    let mut transform = Matrix::identity();
    transform[0][0] = x;
    transform[1][1] = y;
    transform[2][2] = z;
    transform
}

pub fn rotation_x(angle: f32) -> Matrix<4, 4> {
    let mut transform = Matrix::identity();
    transform[1][1] = angle.cos();
    transform[1][2] = -angle.sin();
    transform[2][1] = angle.sin();
    transform[2][2] = angle.cos();
    transform
}

pub fn rotation_y(angle: f32) -> Matrix<4, 4> {
    let mut transform = Matrix::identity();
    transform[0][0] = angle.cos();
    transform[0][2] = angle.sin();
    transform[2][0] = -angle.sin();
    transform[2][2] = angle.cos();
    transform
}

pub fn rotation_z(angle: f32) -> Matrix<4, 4> {
    let mut transform = Matrix::identity();
    transform[0][0] = angle.cos();
    transform[0][1] = -angle.sin();
    transform[1][0] = angle.sin();
    transform[1][1] = angle.cos();
    transform
}

pub fn shearing(x_y: f32, x_z: f32, y_x: f32, y_z: f32, z_x: f32, z_y: f32) -> Matrix<4, 4> {
    let mut transform = Matrix::identity();
    transform[0][1] = x_y;
    transform[0][2] = x_z;
    transform[1][0] = y_x;
    transform[1][2] = y_z;
    transform[2][0] = z_x;
    transform[2][1] = z_y;
    transform
}

pub struct Transform(Matrix<4, 4>);

impl Transform {
    pub fn translation(mut self, x: f32, y: f32, z: f32) -> Self {
        self.0 = translation(x, y, z) * self.0;
        self
    }

    pub fn scaling(mut self, x: f32, y: f32, z: f32) -> Self {
        self.0 = scaling(x, y, z) * self.0;
        self
    }

    pub fn rotation_x(mut self, angle: f32) -> Self {
        self.0 = rotation_x(angle) * self.0;
        self
    }

    pub fn rotation_y(mut self, angle: f32) -> Self {
        self.0 = rotation_y(angle) * self.0;
        self
    }

    pub fn rotation_z(mut self, angle: f32) -> Self {
        self.0 = rotation_z(angle) * self.0;
        self
    }

    pub fn shearing(mut self, x_y: f32, x_z: f32, y_x: f32, y_z: f32, z_x: f32, z_y: f32) -> Self {
        self.0 = shearing(x_y, x_z, y_x, y_z, z_x, z_y) * self.0;
        self
    }
}

impl Default for Transform {
    fn default() -> Self {
        Self(Matrix::identity())
    }
}

impl From<Transform> for Matrix<4, 4> {
    fn from(t: Transform) -> Self {
        t.0
    }
}

impl Mul<Point> for Transform {
    type Output = Point;

    fn mul(self, rhs: Point) -> Self::Output {
        self.0 * rhs
    }
}

impl Mul<Vector> for Transform {
    type Output = Vector;

    fn mul(self, rhs: Vector) -> Self::Output {
        self.0 * rhs
    }
}

pub fn view_transform(from: Point, to: Point, up: Vector) -> Matrix<4, 4> {
    let forward = (to - from).normalize();
    let upn = up.normalize();
    let left = forward.cross(upn);
    let true_up = left.cross(forward);
    let orientation = Matrix::<4, 4>::new([
        [left.x(), left.y(), left.z(), 0.],
        [true_up.x(), true_up.y(), true_up.z(), 0.],
        [-forward.x(), -forward.y(), -forward.z(), 0.],
        [0., 0., 0., 1.],
    ]);
    orientation * translation(-from.x(), -from.y(), -from.z())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    #[test]
    fn translation_basics() {
        let t = translation(5., -3., 2.);
        let p = Point::new(-3., 4., 5.);
        let v = Vector::new(-3., 4., 5.);
        assert!(t * p == Point::new(2., 1., 7.));
        assert!(t.inverse().unwrap() * p == Point::new(-8., 7., 3.));
        assert!(t * v == v);
    }

    #[test]
    fn scaling_basics() {
        let s = scaling(2., 3., 4.);
        let p = Point::new(-4., 6., 8.);
        let v = Vector::new(-4., 6., 8.);
        assert!(s * p == Point::new(-8., 18., 32.));
        assert!(s * v == Vector::new(-8., 18., 32.));
        assert!(s.inverse().unwrap() * v == Vector::new(-2., 2., 2.));

        let s = scaling(-1., 1., 1.);
        let p = Point::new(2., 3., 4.);
        assert!(s * p == Point::new(-2., 3., 4.));
    }

    #[test]
    fn rotations_basics() {
        let p = Point::new(0., 1., 0.);
        let r1 = rotation_x(PI / 4.);
        let r2 = rotation_x(PI / 2.);
        assert!((r1 * p).equal_approx(Point::new(0., f32::sqrt(2.) / 2., f32::sqrt(2.) / 2.)));
        assert!((r2 * p).equal_approx(Point::new(0., 0., 1.)));
        assert!((r1.inverse().unwrap() * p).equal_approx(Point::new(
            0.,
            f32::sqrt(2.) / 2.,
            -f32::sqrt(2.) / 2.
        )));

        let p = Point::new(0., 0., 1.);
        let r1 = rotation_y(PI / 4.);
        let r2 = rotation_y(PI / 2.);
        assert!((r1 * p).equal_approx(Point::new(f32::sqrt(2.) / 2., 0., f32::sqrt(2.) / 2.)));
        assert!((r2 * p).equal_approx(Point::new(1., 0., 0.)));

        let p = Point::new(0., 1., 0.);
        let r1 = rotation_z(PI / 4.);
        let r2 = rotation_z(PI / 2.);
        assert!((r1 * p).equal_approx(Point::new(-f32::sqrt(2.) / 2., f32::sqrt(2.) / 2., 0.)));
        assert!((r2 * p).equal_approx(Point::new(-1., 0., 0.)));
    }

    #[test]
    fn shearing_basics() {
        let p = Point::new(2., 3., 4.);

        let s = shearing(0., 1., 0., 0., 0., 0.);
        assert!(s * p == Point::new(6., 3., 4.));

        let s = shearing(0., 0., 1., 0., 0., 0.);
        assert!(s * p == Point::new(2., 5., 4.));

        let s = shearing(0., 0., 0., 1., 0., 0.);
        assert!(s * p == Point::new(2., 7., 4.));

        let s = shearing(0., 0., 0., 0., 1., 0.);
        assert!(s * p == Point::new(2., 3., 6.));

        let s = shearing(0., 0., 0., 0., 0., 1.);
        assert!(s * p == Point::new(2., 3., 7.));
    }

    #[test]
    fn transformations_chaining() {
        let t1 = rotation_x(PI / 2.);
        let t2 = scaling(5., 5., 5.);
        let t3 = translation(10., 5., 7.);
        let p1 = Point::new(1., 0., 1.);
        let p2 = t1 * p1;
        let p3 = t2 * p2;
        let p4 = t3 * p3;
        assert!(p2.equal_approx(Point::new(1., -1., 0.)));
        assert!(p3.equal_approx(Point::new(5., -5., 0.)));
        assert!(p4.equal_approx(Point::new(15., 0., 7.)));

        let p = Point::new(1., 0., 1.);
        let t = Transform::default()
            .rotation_x(PI / 2.)
            .scaling(5., 5., 5.)
            .translation(10., 5., 7.);
        assert!((t * p).equal_approx(Point::new(15., 0., 7.)));
    }

    #[test]
    fn view_transformations() {
        let from = Point::new(0., 0., 0.);
        let to = Point::new(0., 0., -1.);
        let up = Vector::new(0., 1., 0.);
        let t = view_transform(from, to, up);
        assert!(t == Matrix::identity());

        let from = Point::new(0., 0., 0.);
        let to = Point::new(0., 0., 1.);
        let up = Vector::new(0., 1., 0.);
        let t = view_transform(from, to, up);
        assert!(t == scaling(-1., 1., -1.));

        let from = Point::new(0., 0., 8.);
        let to = Point::new(0., 0., 0.);
        let up = Vector::new(0., 1., 0.);
        let t = view_transform(from, to, up);
        assert!(t == translation(0., 0., -8.));

        let from = Point::new(1., 3., 2.);
        let to = Point::new(4., -2., 8.);
        let up = Vector::new(1., 1., 0.);
        let t = view_transform(from, to, up);
        let exepected = Matrix::<4, 4>::new([
            [-0.50709, 0.50709, 0.67612, -2.36643],
            [0.76772, 0.60609, 0.12122, -2.82843],
            [-0.35857, 0.59761, -0.71714, 0.],
            [0., 0., 0., 1.],
        ]);
        assert!(t.equal_approx(exepected));
    }
}
