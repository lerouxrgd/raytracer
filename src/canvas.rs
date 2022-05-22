use std::io::Write;

use image::{
    codecs::pnm::{PixmapHeader, PnmEncoder, SampleEncoding},
    ColorType, Rgb32FImage,
};

use crate::tuples::Color;

pub struct Canvas {
    buffer: Rgb32FImage,
}

impl Canvas {
    pub fn new(width: u32, height: u32) -> Self {
        Self {
            buffer: Rgb32FImage::new(width, height),
        }
    }

    pub fn width(&self) -> u32 {
        self.buffer.width()
    }

    pub fn height(&self) -> u32 {
        self.buffer.height()
    }

    pub fn write_pixel(&mut self, x: u32, y: u32, color: Color) {
        self.buffer.put_pixel(x, y, color.into())
    }

    pub fn pixel_at(&self, x: u32, y: u32) -> Color {
        self.buffer.get_pixel(x, y).into()
    }

    pub fn to_ppm<W: Write>(&self, mut writer: W) {
        let header = PixmapHeader {
            encoding: SampleEncoding::Ascii,
            height: self.height(),
            width: self.width(),
            maxval: 255,
        };

        let mut encoder = PnmEncoder::new(&mut writer).with_header(header.into());
        encoder
            .encode(
                &self
                    .buffer
                    .pixels()
                    .flat_map(|p| {
                        [
                            (0_f32.max(p.0[0]).min(1.) * 255.).round() as u8,
                            (0_f32.max(p.0[1]).min(1.) * 255.).round() as u8,
                            (0_f32.max(p.0[2]).min(1.) * 255.).round() as u8,
                        ]
                    })
                    .collect::<Vec<_>>()[..],
                self.width(),
                self.height(),
                ColorType::Rgb8,
            )
            .ok();
    }
}

impl From<Color> for image::Rgb<f32> {
    fn from(color: Color) -> Self {
        Self([color.r(), color.g(), color.b()])
    }
}

impl From<&image::Rgb<f32>> for Color {
    fn from(rgb: &image::Rgb<f32>) -> Self {
        Self::new(rgb.0[0], rgb.0[1], rgb.0[2])
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    #[test]
    fn canvas_ops() {
        let mut canvas = Canvas::new(10, 20);

        for i in 0..10 {
            for j in 0..20 {
                assert!(canvas.pixel_at(i, j) == Color::new(0., 0., 0.));
            }
        }

        let red = Color::new(1., 0., 0.);
        canvas.write_pixel(2, 3, red);
        assert!(canvas.pixel_at(2, 3) == red);
    }

    #[test]
    fn canvas_ppm() {
        let mut canvas = Canvas::new(5, 3);
        canvas.write_pixel(0, 0, Color::new(1.5, 0., 0.));
        canvas.write_pixel(2, 1, Color::new(0., 0.5, 0.));
        canvas.write_pixel(4, 2, Color::new(-0.5, 0., 1.));

        let mut out = Cursor::new(Vec::new());
        canvas.to_ppm(&mut out);
        let content = String::from_utf8(out.into_inner());
        assert!(content.is_ok());
    }
}
