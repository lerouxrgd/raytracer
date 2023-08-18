use crate::matrices::Matrix;

#[derive(Debug, Clone, clap::ValueEnum)]
pub enum Dither {
    Bayer2,
    Bayer4,
    Bayer8,
    Bayer16,
    BayerColor,
}

pub fn bayer_matrix<const N: usize>() -> Matrix<N, N> {
    assert!(is_pow2(N), "N must be a power of 2");

    let mut size = 2;
    let mut res = Matrix::<N, N>::identity();
    res[0][0] = 0.;
    res[0][1] = 2.;
    res[1][0] = 3.;
    res[1][1] = 1.;
    for _ in 1..res.size().0.ilog2() as usize {
        for j in 0..size {
            for i in 0..size {
                res[i][j] *= 2_f32.powf(2.);
            }
        }
        for j in 0..size {
            for i in 0..size {
                res[size + i][j] = res[i][j] + 3.;
                res[i][size + j] = res[i][j] + 2.;
                res[size + i][size + j] = res[i][j] + 1.;
            }
        }
        size *= 2;
    }
    (1. / res.size().0.pow(2) as f32) * res
}

pub fn is_pow2(n: usize) -> bool {
    n >= 2 && (n as f32).log2().ceil() == (n as f32).log2().floor()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn check_pow2() {
        assert!(is_pow2(2));
        assert!(is_pow2(4));
        assert!(is_pow2(8));
        assert!(is_pow2(1024));
        assert!(!is_pow2(0));
        assert!(!is_pow2(1));
        assert!(!is_pow2(3));
        assert!(!is_pow2(77));
    }

    #[test]
    fn valid_bayer_matrix() {
        let bayer_2 = bayer_matrix::<2>();
        assert_eq!(bayer_2, (1. / 4.) * Matrix::from([[0.0, 2.0], [3.0, 1.0],]));

        let bayer_4 = bayer_matrix::<4>();
        assert_eq!(
            bayer_4,
            (1. / 16.)
                * Matrix::from([
                    [0.0, 8.0, 2.0, 10.0],
                    [12.0, 4.0, 14.0, 6.0],
                    [3.0, 11.0, 1.0, 9.0],
                    [15.0, 7.0, 13.0, 5.0],
                ])
        );

        let bayer_8 = bayer_matrix::<8>();
        assert_eq!(
            bayer_8,
            (1. / 64.)
                * Matrix::from([
                    [0.0, 32.0, 8.0, 40.0, 2.0, 34.0, 10.0, 42.0],
                    [48.0, 16.0, 56.0, 24.0, 50.0, 18.0, 58.0, 26.0],
                    [12.0, 44.0, 4.0, 36.0, 14.0, 46.0, 6.0, 38.0],
                    [60.0, 28.0, 52.0, 20.0, 62.0, 30.0, 54.0, 22.0],
                    [3.0, 35.0, 11.0, 43.0, 1.0, 33.0, 9.0, 41.0],
                    [51.0, 19.0, 59.0, 27.0, 49.0, 17.0, 57.0, 25.0],
                    [15.0, 47.0, 7.0, 39.0, 13.0, 45.0, 5.0, 37.0],
                    [63.0, 31.0, 55.0, 23.0, 61.0, 29.0, 53.0, 21.0],
                ])
        );
    }
}
