#![allow(clippy::flat_map_identity)]
#![allow(clippy::needless_range_loop)]

use std::ops::{Deref, DerefMut, Mul};

use crate::tuples::{Point, Tuple, Vector};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix<const N: usize, const M: usize>([[f32; M]; N]);

impl<const N: usize, const M: usize> Matrix<N, M> {
    const N: usize = N;
    const M: usize = M;

    pub fn equal_approx(&self, other: Self) -> bool {
        const EPSILON: f32 = 1e-4;
        self.0
            .iter()
            .flat_map(|x| x)
            .zip(other.0.iter().flat_map(|y| y))
            .all(|(x, y)| (x - y).abs() < EPSILON)
    }

    pub fn transpose(&self) -> Matrix<M, N> {
        let mut res = [[0.; N]; M];
        for i in 0..N {
            for j in 0..M {
                res[j][i] = self[i][j]
            }
        }
        Matrix(res)
    }
}

impl<const N: usize, const M: usize> Deref for Matrix<N, M> {
    type Target = [[f32; M]; N];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<const N: usize, const M: usize> DerefMut for Matrix<N, M> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<const N1: usize, const M: usize, const N2: usize> Mul<Matrix<M, N2>> for Matrix<N1, M> {
    type Output = Matrix<N1, N2>;

    fn mul(self, rhs: Matrix<M, N2>) -> Self::Output {
        let mut res = [[0.; N2]; N1];
        for row in 0..N1 {
            for col in 0..N2 {
                res[row][col] = (0..M).map(|m| self[row][m] * rhs[m][col]).sum();
            }
        }
        Matrix(res)
    }
}

impl<const N: usize> Matrix<N, N> {
    pub fn identity() -> Self {
        let mut res = [[0.; N]; N];
        for i in 0..N {
            for j in 0..N {
                if i == j {
                    res[i][j] = 1.;
                }
            }
        }
        Self(res)
    }
}

impl Matrix<2, 2> {
    pub fn new(rows: [[f32; Self::M]; Self::N]) -> Self {
        Self(rows)
    }

    pub fn determinant(&self) -> f32 {
        self[0][0] * self[1][1] - self[1][0] * self[0][1]
    }
}

impl Matrix<3, 3> {
    pub fn new(rows: [[f32; Self::M]; Self::N]) -> Self {
        Self(rows)
    }

    pub fn submatrix(&self, row: usize, col: usize) -> Matrix<{ Self::N - 1 }, { Self::M - 1 }> {
        let mut res = [[0.; Self::M - 1]; Self::N - 1];
        for i in 0..Self::N {
            if i == row {
                continue;
            };
            let ii = if i > row { i - 1 } else { i };
            for j in 0..Self::M {
                if j == col {
                    continue;
                }
                let jj = if j > col { j - 1 } else { j };
                res[ii][jj] = self[i][j]
            }
        }
        Matrix(res)
    }

    pub fn minor(&self, row: usize, col: usize) -> f32 {
        self.submatrix(row, col).determinant()
    }

    pub fn cofactor(&self, row: usize, col: usize) -> f32 {
        if (row + col) % 2 == 0 {
            self.minor(row, col)
        } else {
            -self.minor(row, col)
        }
    }

    pub fn determinant(&self) -> f32 {
        let mut det = 0.;
        for j in 0..Self::M {
            det += self[0][j] * self.cofactor(0, j);
        }
        det
    }
}

impl Matrix<4, 4> {
    pub fn new(rows: [[f32; Self::M]; Self::N]) -> Self {
        Self(rows)
    }

    pub fn submatrix(&self, row: usize, col: usize) -> Matrix<{ Self::N - 1 }, { Self::M - 1 }> {
        let mut res = [[0.; Self::M - 1]; Self::N - 1];
        for i in 0..Self::N {
            if i == row {
                continue;
            };
            let ii = if i > row { i - 1 } else { i };
            for j in 0..Self::M {
                if j == col {
                    continue;
                }
                let jj = if j > col { j - 1 } else { j };
                res[ii][jj] = self[i][j]
            }
        }
        Matrix(res)
    }

    pub fn minor(&self, row: usize, col: usize) -> f32 {
        self.submatrix(row, col).determinant()
    }

    pub fn cofactor(&self, row: usize, col: usize) -> f32 {
        if (row + col) % 2 == 0 {
            self.minor(row, col)
        } else {
            -self.minor(row, col)
        }
    }

    pub fn determinant(&self) -> f32 {
        let mut det = 0.;
        for j in 0..Self::M {
            det += self[0][j] * self.cofactor(0, j);
        }
        det
    }

    pub fn inverse(&self) -> Option<Self> {
        let det = self.determinant();
        if det == 0. {
            return None;
        }
        let mut res = [[0.; Self::M]; Self::N];
        for i in 0..Self::N {
            for j in 0..Self::M {
                res[j][i] = self.cofactor(i, j) / det;
            }
        }
        Some(Self(res))
    }
}

impl Mul<Tuple> for Matrix<4, 4> {
    type Output = Tuple;

    fn mul(self, rhs: Tuple) -> Self::Output {
        Tuple::from(self * Matrix::<4, 1>::from(rhs))
    }
}

impl From<Tuple> for Matrix<4, 1> {
    fn from(t: Tuple) -> Self {
        Self([[t[0]], [t[1]], [t[2]], [t[3]]])
    }
}

impl From<Matrix<4, 1>> for Tuple {
    fn from(m: Matrix<4, 1>) -> Self {
        Tuple::from([m[0][0], m[1][0], m[2][0], m[3][0]])
    }
}

impl Mul<Point> for Matrix<4, 4> {
    type Output = Point;

    fn mul(self, rhs: Point) -> Self::Output {
        let t = self * Tuple::from(rhs);
        Point::new(t[0], t[1], t[2])
    }
}

impl Mul<Vector> for Matrix<4, 4> {
    type Output = Vector;

    fn mul(self, rhs: Vector) -> Self::Output {
        let t = self * Tuple::from(rhs);
        Vector::new(t[0], t[1], t[2])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matrix_basics() {
        let m1 = Matrix::<3, 3>::new([[0., 1., 2.], [3., 4., 5.], [6., 7., 8.]]);
        let m2 = Matrix::<3, 3>::new([[0., 1., 2.], [3., 4., 5.], [6., 7., 8.]]);
        let m3 = Matrix::<3, 3>::new([[6., 7., 8.], [3., 4., 5.], [0., 1., 2.]]);
        assert!(m1[0][0] == 0.);
        assert!(m1[1][0] == 3.);
        assert!(m1[2][1] == 7.);
        assert!(m1.equal_approx(m2));
        assert!(!m1.equal_approx(m3));

        let m1 = Matrix::<4, 4>::new([
            [1., 2., 3., 4.],
            [5., 6., 7., 8.],
            [9., 8., 7., 6.],
            [5., 4., 3., 2.],
        ]);
        let m2 = Matrix::<4, 4>::new([
            [-2., 1., 2., 3.],
            [3., 2., 1., -1.],
            [4., 3., 6., 5.],
            [1., 2., 7., 8.],
        ]);
        let m3 = Matrix::<4, 4>::new([
            [20., 22., 50., 48.],
            [44., 54., 114., 108.],
            [40., 58., 110., 102.],
            [16., 26., 46., 42.],
        ]);
        assert!(m1 * m2 == m3);

        let m1 = Matrix::<4, 4>::new([
            [1., 2., 3., 4.],
            [2., 4., 4., 2.],
            [8., 6., 4., 1.],
            [0., 0., 0., 1.],
        ]);
        let t = Tuple::from([1., 2., 3., 1.]);
        assert!(m1 * t == Tuple::from([18., 24., 33., 1.]));

        let m = Matrix::<4, 4>::new([
            [1., 2., 3., 4.],
            [5., 6., 7., 8.],
            [9., 8., 7., 6.],
            [5., 4., 3., 2.],
        ]);
        let t = Tuple::from([1., 2., 3., 1.]);
        assert!(m * Matrix::identity() == m);
        assert!(Matrix::identity() * m == m);
        assert!(Matrix::identity() * t == t);

        let m1 = Matrix::<4, 4>::new([
            [1., 2., 3., 4.],
            [5., 6., 7., 8.],
            [9., 8., 7., 6.],
            [5., 4., 3., 2.],
        ]);
        let m2 = Matrix::<4, 4>::new([
            [1., 5., 9., 5.],
            [2., 6., 8., 4.],
            [3., 7., 7., 3.],
            [4., 8., 6., 2.],
        ]);
        assert!(m1.transpose() == m2);
        assert!(Matrix::<4, 4>::identity().transpose() == Matrix::identity());
    }

    #[test]
    fn matrix_advanced() {
        let m = Matrix::<2, 2>::new([[1., 5.], [-3., 2.]]);
        assert!(m.determinant() == 17.);

        let m1 = Matrix::<4, 4>::new([
            [1., 2., 3., 4.],
            [5., 6., 7., 8.],
            [9., 8., 7., 6.],
            [5., 4., 3., 2.],
        ]);
        let m2 = Matrix::<3, 3>::new([[1., 2., 4.], [9., 8., 6.], [5., 4., 2.]]);
        let m3 = Matrix::<2, 2>::new([[2., 4.], [8., 6.]]);
        assert!(m1.submatrix(1, 2) == m2);
        assert!(m2.submatrix(2, 0) == m3);

        let m1 = Matrix::<3, 3>::new([[3., 5., 0.], [2., -1., -7.], [6., -1., 5.]]);
        let m2 = m1.submatrix(1, 0);
        assert!(m2.determinant() == 25.);
        assert!(m1.minor(1, 0) == 25.);

        let m = Matrix::<3, 3>::new([[3., 5., 0.], [2., -1., -7.], [6., -1., 5.]]);
        assert!(m.minor(0, 0) == -12.);
        assert!(m.cofactor(0, 0) == -12.);
        assert!(m.minor(1, 0) == 25.);
        assert!(m.cofactor(1, 0) == -25.);

        let m = Matrix::<3, 3>::new([[1., 2., 6.], [-5., 8., -4.], [2., 6., 4.]]);
        assert!(m.cofactor(0, 0) == 56.);
        assert!(m.cofactor(0, 1) == 12.);
        assert!(m.cofactor(0, 2) == -46.);
        assert!(m.determinant() == -196.);

        let m = Matrix::<4, 4>::new([
            [-2., -8., 3., 5.],
            [-3., 1., 7., 3.],
            [1., 2., -9., 6.],
            [-6., 7., 7., -9.],
        ]);
        assert!(m.cofactor(0, 0) == 690.);
        assert!(m.cofactor(0, 1) == 447.);
        assert!(m.cofactor(0, 2) == 210.);
        assert!(m.cofactor(0, 3) == 51.);
        assert!(m.determinant() == -4071.);

        let m1 = Matrix::<4, 4>::new([
            [-5., 2., 6., -8.],
            [1., -5., 1., 8.],
            [7., 7., -6., -7.],
            [1., -3., 7., 4.],
        ]);
        let m2 = m1.inverse().unwrap();
        let m3 = Matrix::<4, 4>::new([
            [0.21805, 0.45113, 0.24060, -0.04511],
            [-0.80827, -1.45677, -0.44361, 0.52068],
            [-0.07895, -0.22368, -0.05263, 0.19737],
            [-0.52256, -0.81391, -0.30075, 0.30639],
        ]);
        assert!(m1.determinant() == 532.);
        assert!(m1.cofactor(2, 3) == -160.);
        assert!(m2[3][2] == -160. / 532.);
        assert!(m1.cofactor(3, 2) == 105.);
        assert!(m2[2][3] == 105. / 532.);
        assert!(m2.equal_approx(m3));

        let m1 = Matrix::<4, 4>::new([
            [1., 2., 3., 4.],
            [5., 6., 7., 8.],
            [9., 8., 7., 6.],
            [5., 4., 3., 2.],
        ]);
        let m2 = Matrix::<4, 4>::new([
            [-2., 1., 2., 3.],
            [3., 2., 1., -1.],
            [4., 3., 6., 5.],
            [1., 2., 7., 8.],
        ]);
        let m3 = m1 * m2;
        assert!(m1.equal_approx(m3 * m2.inverse().unwrap()));
    }
}
