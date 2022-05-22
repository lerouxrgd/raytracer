use std::ops::{Add, Deref, DerefMut, Div, Mul, Neg, Sub};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tuple {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub w: f32,
}

impl Tuple {
    pub fn new(x: f32, y: f32, z: f32, w: f32) -> Self {
        Self { x, y, z, w }
    }

    pub fn equal_approx(&self, other: Self) -> bool {
        const EPSILON: f32 = 1e-5;
        (self.x - other.x).abs() < EPSILON
            && (self.y - other.y).abs() < EPSILON
            && (self.z - other.z).abs() < EPSILON
            && (self.w - other.w).abs() < EPSILON
    }
}

impl Add for Tuple {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
            w: self.w + rhs.w,
        }
    }
}

impl Sub for Tuple {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
            w: self.w - rhs.w,
        }
    }
}

impl Neg for Tuple {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: -self.w,
        }
    }
}

impl Mul<Tuple> for f32 {
    type Output = Tuple;

    fn mul(self, rhs: Tuple) -> Self::Output {
        Tuple {
            x: self * rhs.x,
            y: self * rhs.y,
            z: self * rhs.z,
            w: self * rhs.w,
        }
    }
}

impl Div<f32> for Tuple {
    type Output = Self;

    fn div(self, rhs: f32) -> Self::Output {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
            w: self.w / rhs,
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point(Tuple);

impl Point {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self(Tuple::new(x, y, z, 1.))
    }

    pub fn equal_approx(&self, other: Self) -> bool {
        self.0.equal_approx(other.0)
    }
}

impl Deref for Point {
    type Target = Tuple;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Point {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Sub for Point {
    type Output = Vector;

    fn sub(self, rhs: Self) -> Self::Output {
        Vector(self.0 - rhs.0)
    }
}

impl Sub<Vector> for Point {
    type Output = Self;

    fn sub(self, rhs: Vector) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

////////////////////////////////////////////////////////////////////////////////////////

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vector(Tuple);

impl Vector {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self(Tuple::new(x, y, z, 0.))
    }

    pub fn equal_approx(&self, other: Self) -> bool {
        self.0.equal_approx(other.0)
    }

    pub fn magnitude(&self) -> f32 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2) + self.w.powi(2)).sqrt()
    }

    pub fn normalize(&self) -> Self {
        Self(self.0 / self.magnitude())
    }

    pub fn dot(&self, other: Self) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z + self.w * other.w
    }

    pub fn cross(&self, other: Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }
}

impl Deref for Vector {
    type Target = Tuple;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Vector {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Add for Vector {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl Sub for Vector {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl Neg for Vector {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

impl Mul<Vector> for f32 {
    type Output = Vector;

    fn mul(self, rhs: Vector) -> Self::Output {
        Vector(self * rhs.0)
    }
}

impl Div<f32> for Vector {
    type Output = Self;

    fn div(self, rhs: f32) -> Self::Output {
        Self(self.0 / rhs)
    }
}

////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_ops() {
        let t1 = Tuple::new(3., -2., 5., 1.);
        let t2 = Tuple::new(-2., 3., 1., 0.);
        assert!(t1 + t2 == Tuple::new(1., 1., 6., 1.));

        let p1 = Point::new(3., 2., 1.);
        let p2 = Point::new(5., 6., 7.);
        assert!(p1 - p2 == Vector::new(-2., -4., -6.));

        let p = Point::new(3., 2., 1.);
        let v = Vector::new(5., 6., 7.);
        assert!(p - v == Point::new(-2., -4., -6.));

        let v1 = Vector::new(3., 2., 1.);
        let v2 = Vector::new(5., 6., 7.);
        assert!(v1 - v2 == Vector::new(-2., -4., -6.));

        let t = Tuple::new(1., -2., 3., -4.);
        assert!(-t == Tuple::new(-1., 2., -3., 4.));
        let v = Vector::new(1., -2., 3.);
        assert!(-v == Vector::new(-1., 2., -3.));

        let t = Tuple::new(1., -2., 3., -4.);
        assert!(3.5 * t == Tuple::new(3.5, -7., 10.5, -14.));
        let v = Vector::new(1., -2., 3.);
        assert!(3.5 * v == Vector::new(3.5, -7., 10.5));

        let t = Tuple::new(1., -2., 3., -4.);
        assert!(t / 2. == Tuple::new(0.5, -1., 1.5, -2.));
        let v = Vector::new(1., -2., 3.);
        assert!(v / 2. == Vector::new(0.5, -1., 1.5));

        assert!(Vector::new(1., 0., 0.).magnitude() == 1.);
        assert!(Vector::new(0., 1., 0.).magnitude() == 1.);
        assert!(Vector::new(0., 0., 1.).magnitude() == 1.);
        assert!(Vector::new(1., 2., 3.).magnitude() == 14_f32.sqrt());
        assert!(Vector::new(-1., -2., -3.).magnitude() == 14_f32.sqrt());

        assert!(Vector::new(4., 0., 0.).normalize() == Vector::new(1., 0., 0.));
        let v1 = Vector::new(1., 2., 3.);
        let v2 = Vector::new(1. / 14_f32.sqrt(), 2. / 14_f32.sqrt(), 3. / 14_f32.sqrt());
        assert!(v1.normalize().equal_approx(v2));

        let v1 = Vector::new(1., 2., 3.);
        let v2 = Vector::new(2., 3., 4.);
        assert!(v1.dot(v2) == 20.);

        let v1 = Vector::new(1., 2., 3.);
        let v2 = Vector::new(2., 3., 4.);
        assert!(v1.cross(v2) == Vector::new(-1., 2., -1.));
        assert!(v2.cross(v1) == Vector::new(1., -2., 1.));
    }
}
