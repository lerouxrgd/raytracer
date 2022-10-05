use lazy_static::lazy_static;
use parking_lot::RwLock;
use slotmap::{new_key_type, SlotMap};

use crate::{canvas::Canvas, tuples::Color};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum UvPattern {
    Checker(UvChecker),
    AlignCheck(UvAlignCheck),
    Image(UvImage),
}

impl UvPattern {
    pub fn uv_pattern_at(&self, u: f32, v: f32) -> Color {
        match self {
            Self::Checker(p) => p.uv_pattern_at(u, v),
            Self::AlignCheck(p) => p.uv_pattern_at(u, v),
            Self::Image(p) => p.uv_pattern_at(u, v),
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UvChecker {
    pub width: f32,
    pub height: f32,
    pub a: Color,
    pub b: Color,
}

impl UvChecker {
    pub fn new(width: f32, height: f32, a: Color, b: Color) -> Self {
        Self {
            width,
            height,
            a,
            b,
        }
    }

    pub fn uv_pattern_at(&self, u: f32, v: f32) -> Color {
        let u2 = (u * self.width).floor();
        let v2 = (v * self.height).floor();
        if (u2 + v2) % 2. == 0. {
            self.a
        } else {
            self.b
        }
    }
}

impl From<UvChecker> for UvPattern {
    fn from(uv_checker: UvChecker) -> Self {
        UvPattern::Checker(uv_checker)
    }
}

////////////////////////////////////////////////////////////////////////////////////////

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UvAlignCheck {
    main: Color,
    ul: Color,
    ur: Color,
    bl: Color,
    br: Color,
}

impl From<UvAlignCheck> for UvPattern {
    fn from(uv_align_check: UvAlignCheck) -> Self {
        UvPattern::AlignCheck(uv_align_check)
    }
}

impl UvAlignCheck {
    pub fn new(main: Color, ul: Color, ur: Color, bl: Color, br: Color) -> Self {
        Self {
            main,
            ul,
            ur,
            bl,
            br,
        }
    }

    pub fn uv_pattern_at(&self, u: f32, v: f32) -> Color {
        if v > 0.8 {
            if u < 0.2 {
                return self.ul;
            } else if u > 0.8 {
                return self.ur;
            }
        } else if v < 0.2 {
            if u < 0.2 {
                return self.bl;
            } else if u > 0.8 {
                return self.br;
            }
        }
        self.main
    }
}

////////////////////////////////////////////////////////////////////////////////////////

lazy_static! {
    static ref UV_IMAGES: RwLock<SlotMap<UvImage, Canvas>> = RwLock::new(SlotMap::with_key());
}

new_key_type! {
    pub struct UvImage;
}

impl UvImage {
    pub fn new(canvas: Canvas) -> Self {
        UV_IMAGES.write().insert(canvas)
    }

    pub fn uv_pattern_at(&self, u: f32, v: f32) -> Color {
        // flip v over so it matches the image layout, with y at the top
        let v = 1. - v;

        let images = UV_IMAGES.read();
        let canvas = images.get(*self).unwrap();
        let x = u * (canvas.width() as f32 - 1.);
        let y = v * (canvas.height() as f32 - 1.);

        canvas.pixel_at(x.round() as u32, y.round() as u32)
    }
}

impl From<UvImage> for UvPattern {
    fn from(uv_image: UvImage) -> Self {
        UvPattern::Image(uv_image)
    }
}

////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    #[test]
    fn uv_checkers_basics() {
        let pattern = UvChecker::new(2., 2., Color::black(), Color::white());
        let cases = vec![
            (0.0, 0.0, Color::black()),
            (0.5, 0.0, Color::white()),
            (0.0, 0.5, Color::white()),
            (0.5, 0.5, Color::black()),
            (1.0, 1.0, Color::black()),
        ];
        for (u, v, res) in cases {
            assert!(pattern.uv_pattern_at(u, v) == res);
        }
    }

    #[test]
    fn uv_align_check_basics() {
        let main = Color::white();
        let ul = Color::new(1., 0., 0.);
        let ur = Color::new(1., 1., 0.);
        let bl = Color::new(0., 1., 0.);
        let br = Color::new(0., 1., 1.);
        let pattern = UvAlignCheck::new(main, ul, ur, bl, br);
        let cases = vec![
            (0.5, 0.5, main),
            (0.1, 0.9, ul),
            (0.9, 0.9, ur),
            (0.1, 0.1, bl),
            (0.9, 0.1, br),
        ];
        for (u, v, res) in cases {
            assert!(pattern.uv_pattern_at(u, v) == res);
        }
    }

    #[test]
    fn uv_image_basics() {
        let ppm_data = b"P3
10 10
10
0 0 0  1 1 1  2 2 2  3 3 3  4 4 4  5 5 5  6 6 6  7 7 7  8 8 8  9 9 9
1 1 1  2 2 2  3 3 3  4 4 4  5 5 5  6 6 6  7 7 7  8 8 8  9 9 9  0 0 0
2 2 2  3 3 3  4 4 4  5 5 5  6 6 6  7 7 7  8 8 8  9 9 9  0 0 0  1 1 1
3 3 3  4 4 4  5 5 5  6 6 6  7 7 7  8 8 8  9 9 9  0 0 0  1 1 1  2 2 2
4 4 4  5 5 5  6 6 6  7 7 7  8 8 8  9 9 9  0 0 0  1 1 1  2 2 2  3 3 3
5 5 5  6 6 6  7 7 7  8 8 8  9 9 9  0 0 0  1 1 1  2 2 2  3 3 3  4 4 4
6 6 6  7 7 7  8 8 8  9 9 9  0 0 0  1 1 1  2 2 2  3 3 3  4 4 4  5 5 5
7 7 7  8 8 8  9 9 9  0 0 0  1 1 1  2 2 2  3 3 3  4 4 4  5 5 5  6 6 6
8 8 8  9 9 9  0 0 0  1 1 1  2 2 2  3 3 3  4 4 4  5 5 5  6 6 6  7 7 7
9 9 9  0 0 0  1 1 1  2 2 2  3 3 3  4 4 4  5 5 5  6 6 6  7 7 7  8 8 8";

        let reader = Cursor::new(ppm_data);
        let canvas = Canvas::from_ppm(reader).unwrap();
        let pattern = UvImage::new(canvas);
        let cases = vec![
            (0., 0., Color::new(0.9, 0.9, 0.9)),
            (0.3, 0., Color::new(0.2, 0.2, 0.2)),
            (0.6, 0.3, Color::new(0.1, 0.1, 0.1)),
            (1., 1., Color::new(0.9, 0.9, 0.9)),
        ];
        for (u, v, res) in cases {
            assert!(pattern.uv_pattern_at(u, v) == res);
        }
    }
}
