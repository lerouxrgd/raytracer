use std::error::Error;
use std::io::{BufRead, Seek, Write};

use image::codecs::pnm::{PixmapHeader, PnmDecoder, PnmEncoder, SampleEncoding};
use image::{ColorType, ImageFormat, Rgb32FImage};

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

    pub fn from_ppm<R>(reader: R) -> Result<Self, Box<dyn Error>>
    where
        R: BufRead + Seek,
    {
        let (mut reader, header) = PnmDecoder::new(reader)?.into_inner();
        reader.rewind()?;

        let max_val = header.maximal_sample() as f32;
        let mut buffer = Rgb32FImage::new(header.width(), header.height());

        let pixels = image::io::Reader::with_format(reader, ImageFormat::Pnm)
            .decode()?
            .into_rgb8();

        for (x, y, rgb) in pixels.enumerate_pixels() {
            let (r, g, b) = (rgb.0[0], rgb.0[1], rgb.0[2]);
            let (r, g, b) = (r as f32 / max_val, g as f32 / max_val, b as f32 / max_val);
            buffer.put_pixel(x, y, image::Rgb([r, g, b]));
        }

        Ok(Self { buffer })
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
    fn canvas_basics() {
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
    fn canvas_to_ppm() {
        let mut canvas = Canvas::new(5, 3);
        canvas.write_pixel(0, 0, Color::new(1.5, 0., 0.));
        canvas.write_pixel(2, 1, Color::new(0., 0.5, 0.));
        canvas.write_pixel(4, 2, Color::new(-0.5, 0., 1.));

        let mut out = Cursor::new(Vec::new());
        canvas.to_ppm(&mut out);
        let content = String::from_utf8(out.into_inner());
        assert!(content.is_ok());
    }

    #[test]
    fn canvas_from_ppm() {
        let ppm_data = b"P3
2 2
100
100 100 100  50 50 50
75 50 25  0 0 0
";
        let reader = Cursor::new(ppm_data);
        let canvas = Canvas::from_ppm(reader).unwrap();
        assert!(dbg!(canvas.pixel_at(0, 1)) == Color::new(0.75, 0.5, 0.25));
    }
}
