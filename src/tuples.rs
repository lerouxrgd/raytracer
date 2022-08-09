use std::ops::{Add, Deref, Div, Mul, Neg, Sub};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tuple([f32; 4]);

impl Tuple {
    pub fn equal_approx(&self, other: Self) -> bool {
        const EPSILON: f32 = 1e-4;
        (self.0[0] - other.0[0]).abs() < EPSILON
            && (self.0[1] - other.0[1]).abs() < EPSILON
            && (self.0[2] - other.0[2]).abs() < EPSILON
            && (self.0[3] - other.0[3]).abs() < EPSILON
    }
}

impl From<[f32; 4]> for Tuple {
    fn from(t: [f32; 4]) -> Self {
        Self(t)
    }
}

impl Deref for Tuple {
    type Target = [f32; 4];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl Add for Tuple {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self([
            self.0[0] + rhs.0[0],
            self.0[1] + rhs.0[1],
            self.0[2] + rhs.0[2],
            self.0[3] + rhs.0[3],
        ])
    }
}

impl Sub for Tuple {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self([
            self.0[0] - rhs.0[0],
            self.0[1] - rhs.0[1],
            self.0[2] - rhs.0[2],
            self.0[3] - rhs.0[3],
        ])
    }
}

impl Neg for Tuple {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self([-self.0[0], -self.0[1], -self.0[2], -self.0[3]])
    }
}

impl Mul for Tuple {
    type Output = Tuple;

    fn mul(self, rhs: Tuple) -> Self::Output {
        Tuple([
            self.0[0] * rhs.0[0],
            self.0[1] * rhs.0[1],
            self.0[2] * rhs.0[2],
            self.0[3] * rhs.0[3],
        ])
    }
}

impl Mul<Tuple> for f32 {
    type Output = Tuple;

    fn mul(self, rhs: Tuple) -> Self::Output {
        Tuple([
            self * rhs.0[0],
            self * rhs.0[1],
            self * rhs.0[2],
            self * rhs.0[3],
        ])
    }
}

impl Div<f32> for Tuple {
    type Output = Self;

    fn div(self, rhs: f32) -> Self::Output {
        Self([
            self.0[0] / rhs,
            self.0[1] / rhs,
            self.0[2] / rhs,
            self.0[3] / rhs,
        ])
    }
}

////////////////////////////////////////////////////////////////////////////////////////

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point(Tuple);

impl Point {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self(Tuple([x, y, z, 1.]))
    }

    pub fn x(&self) -> f32 {
        self.0[0]
    }

    pub fn y(&self) -> f32 {
        self.0[1]
    }

    pub fn z(&self) -> f32 {
        self.0[2]
    }

    pub fn w(&self) -> f32 {
        self.0[3]
    }

    pub fn equal_approx(&self, other: Self) -> bool {
        self.0.equal_approx(other.0)
    }
}

impl From<[f32; 3]> for Point {
    fn from(t: [f32; 3]) -> Self {
        Self::new(t[0], t[1], t[2])
    }
}

impl Sub for Point {
    type Output = Vector;

    fn sub(self, rhs: Self) -> Self::Output {
        Vector(self.0 - rhs.0)
    }
}

impl Add<Vector> for Point {
    type Output = Self;

    fn add(self, rhs: Vector) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl Sub<Vector> for Point {
    type Output = Self;

    fn sub(self, rhs: Vector) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl From<Point> for Tuple {
    fn from(p: Point) -> Self {
        Tuple(p.0 .0)
    }
}

////////////////////////////////////////////////////////////////////////////////////////

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vector(Tuple);

impl Vector {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self(Tuple([x, y, z, 0.]))
    }

    pub fn x(&self) -> f32 {
        self.0[0]
    }

    pub fn y(&self) -> f32 {
        self.0[1]
    }

    pub fn z(&self) -> f32 {
        self.0[2]
    }

    pub fn w(&self) -> f32 {
        self.0[3]
    }

    pub fn equal_approx(&self, other: Vector) -> bool {
        self.0.equal_approx(other.0)
    }

    pub fn magnitude(&self) -> f32 {
        (self.x().powi(2) + self.y().powi(2) + self.z().powi(2) + self.w().powi(2)).sqrt()
    }

    pub fn normalize(&self) -> Vector {
        Vector(self.0 / self.magnitude())
    }

    pub fn dot(&self, other: Vector) -> f32 {
        self.x() * other.x() + self.y() * other.y() + self.z() * other.z() + self.w() * other.w()
    }

    pub fn cross(&self, other: Vector) -> Vector {
        Vector::new(
            self.y() * other.z() - self.z() * other.y(),
            self.z() * other.x() - self.x() * other.z(),
            self.x() * other.y() - self.y() * other.x(),
        )
    }

    pub fn reflect(&self, normal: Vector) -> Vector {
        *self - 2. * self.dot(normal) * normal
    }
}

impl From<[f32; 3]> for Vector {
    fn from(t: [f32; 3]) -> Self {
        Self::new(t[0], t[1], t[2])
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

impl From<Vector> for Tuple {
    fn from(v: Vector) -> Self {
        Tuple(v.0 .0)
    }
}

////////////////////////////////////////////////////////////////////////////////////////

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Color(Tuple);

impl Color {
    pub fn new(r: f32, g: f32, b: f32) -> Self {
        Self(Tuple([r, g, b, 0.]))
    }

    pub fn white() -> Self {
        Self::new(1., 1., 1.)
    }

    pub fn black() -> Self {
        Self::new(0., 0., 0.)
    }

    pub fn r(&self) -> f32 {
        self.0[0]
    }

    pub fn g(&self) -> f32 {
        self.0[1]
    }

    pub fn b(&self) -> f32 {
        self.0[2]
    }

    pub fn r_u8(&self) -> u8 {
        (0_f32.max(self.0[0]).min(1.) * 255.).round() as u8
    }

    pub fn g_u8(&self) -> u8 {
        (0_f32.max(self.0[1]).min(1.) * 255.).round() as u8
    }

    pub fn b_u8(&self) -> u8 {
        (0_f32.max(self.0[2]).min(1.) * 255.).round() as u8
    }

    pub fn equal_approx(&self, other: Self) -> bool {
        self.0.equal_approx(other.0)
    }
}

impl From<[f32; 3]> for Color {
    fn from(t: [f32; 3]) -> Self {
        Self::new(t[0], t[1], t[2])
    }
}

impl Add for Color {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl Sub for Color {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl Mul for Color {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}

impl Mul<f32> for Color {
    type Output = Self;

    fn mul(self, rhs: f32) -> Self::Output {
        Self(rhs * self.0)
    }
}

////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tuples_basics() {
        let t1 = Tuple([3., -2., 5., 1.]);
        let t2 = Tuple([-2., 3., 1., 0.]);
        assert!(t1 + t2 == Tuple([1., 1., 6., 1.]));

        let p1 = Point::new(3., 2., 1.);
        let p2 = Point::new(5., 6., 7.);
        assert!(p1 - p2 == Vector::new(-2., -4., -6.));

        let p = Point::new(3., 2., 1.);
        let v = Vector::new(5., 6., 7.);
        assert!(p - v == Point::new(-2., -4., -6.));

        let v1 = Vector::new(3., 2., 1.);
        let v2 = Vector::new(5., 6., 7.);
        assert!(v1 - v2 == Vector::new(-2., -4., -6.));

        let t = Tuple([1., -2., 3., -4.]);
        assert!(-t == Tuple([-1., 2., -3., 4.]));
        let v = Vector::new(1., -2., 3.);
        assert!(-v == Vector::new(-1., 2., -3.));

        let t = Tuple([1., -2., 3., -4.]);
        assert!(3.5 * t == Tuple([3.5, -7., 10.5, -14.]));
        let v = Vector::new(1., -2., 3.);
        assert!(3.5 * v == Vector::new(3.5, -7., 10.5));

        let t = Tuple([1., -2., 3., -4.]);
        assert!(t / 2. == Tuple([0.5, -1., 1.5, -2.]));
        let v = Vector::new(1., -2., 3.);
        assert!(v / 2. == Vector::new(0.5, -1., 1.5));

        assert!(Vector::new(1., 0., 0.).magnitude() == 1.);
        assert!(Vector::new(0., 1., 0.).magnitude() == 1.);
        assert!(Vector::new(0., 0., 1.).magnitude() == 1.);
        assert!(Vector::new(1., 2., 3.).magnitude() == f32::sqrt(14.));
        assert!(Vector::new(-1., -2., -3.).magnitude() == f32::sqrt(14.));

        assert!(Vector::new(4., 0., 0.).normalize() == Vector::new(1., 0., 0.));
        let v1 = Vector::new(1., 2., 3.);
        let v2 = Vector::new(
            1. / f32::sqrt(14.),
            2. / f32::sqrt(14.),
            3. / f32::sqrt(14.),
        );
        assert!(v1.normalize().equal_approx(v2));

        let v1 = Vector::new(1., 2., 3.);
        let v2 = Vector::new(2., 3., 4.);
        assert!(v1.dot(v2) == 20.);

        let v1 = Vector::new(1., 2., 3.);
        let v2 = Vector::new(2., 3., 4.);
        assert!(v1.cross(v2) == Vector::new(-1., 2., -1.));
        assert!(v2.cross(v1) == Vector::new(1., -2., 1.));

        let c = Color::new(0.2, 0.3, 0.4);
        assert!(c * 2. == Color::new(0.4, 0.6, 0.8));

        let c1 = Color::new(0.9, 0.6, 0.75);
        let c2 = Color::new(0.7, 0.1, 0.25);
        assert!((c1 + c2).equal_approx(Color::new(1.6, 0.7, 1.)));
        assert!((c1 - c2).equal_approx(Color::new(0.2, 0.5, 0.5)));

        let c1 = Color::new(1., 0.2, 0.4);
        let c2 = Color::new(0.9, 1., 0.1);
        assert!((c1 * c2).equal_approx(Color::new(0.9, 0.2, 0.04)));
    }

    #[test]
    fn vector_reflect() {
        let v = Vector::new(1., -1., 0.);
        let n = Vector::new(0., 1., 0.);
        assert!(v.reflect(n) == Vector::new(1., 1., 0.));

        let v = Vector::new(0., -1., 0.);
        let n = Vector::new(f32::sqrt(2.) / 2., f32::sqrt(2.) / 2., 0.);
        assert!(v.reflect(n).equal_approx(Vector::new(1., 0., 0.)));
    }
}
