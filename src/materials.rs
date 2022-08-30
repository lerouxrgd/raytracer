use crate::lights::Light;
use crate::patterns::Pattern;
use crate::shapes::Shape;
use crate::tuples::{Color, Point, Vector};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Material {
    pub pattern: Option<Pattern>,
    pub color: Color,
    pub ambient: f32,
    pub diffuse: f32,
    pub specular: f32,
    pub shininess: f32,
    pub reflective: f32,
    pub transparency: f32,
    pub refractive_index: f32,
}

impl Default for Material {
    fn default() -> Self {
        Material {
            pattern: None,
            color: Color::white(),
            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.9,
            shininess: 200.,
            reflective: 0.,
            transparency: 0.,
            refractive_index: 1.,
        }
    }
}

impl Material {
    pub fn pattern<P: Into<Pattern>>(mut self, pattern: P) -> Self {
        self.pattern = Some(pattern.into());
        self
    }

    pub fn color(mut self, color: Color) -> Self {
        self.color = color;
        self
    }

    pub fn ambient(mut self, ambient: f32) -> Self {
        self.ambient = ambient;
        self
    }

    pub fn diffuse(mut self, diffuse: f32) -> Self {
        self.diffuse = diffuse;
        self
    }

    pub fn specular(mut self, specular: f32) -> Self {
        self.specular = specular;
        self
    }

    pub fn shininess(mut self, shininess: f32) -> Self {
        self.shininess = shininess;
        self
    }

    pub fn reflective(mut self, reflective: f32) -> Self {
        self.reflective = reflective;
        self
    }

    pub fn transparency(mut self, transparency: f32) -> Self {
        self.transparency = transparency;
        self
    }

    pub fn refractive_index(mut self, refractive_index: f32) -> Self {
        self.refractive_index = refractive_index;
        self
    }
}

pub fn lighting(
    material: Material,
    shape: Shape,
    light: Light,
    point: Point,
    eyev: Vector,
    normalv: Vector,
    light_intensity: f32,
) -> Color {
    let color = match material.pattern {
        Some(pattern) => pattern.pattern_at_shape(shape, point),
        None => material.color,
    };

    // Combine the surface color with the light's color/intensity
    let effective_color = color * light.intensity();

    // Compute the ambient contribution
    let ambient = effective_color * material.ambient;

    match light {
        Light::PointLight(_) => {
            // Find the direction to the light source
            let lightv = (light.position() - point).normalize();

            // `light_dot_normal` represents the cosine of the angle between the light
            // vector and the normal vector. A negative number means the light is on the
            // other side of the surface.
            let light_dot_normal = lightv.dot(normalv);

            let (diffuse, specular) = if light_dot_normal >= 0. {
                // Compute the diffuse contribution
                let diffuse = effective_color * material.diffuse * light_dot_normal;

                // `reflect_dot_eye` represents the cosine of the angle between the
                // reflection vector and the eye vector. A negative value means that the
                // lights reflects away from the eye.
                let reflectv = (-lightv).reflect(normalv);
                let reflect_dot_eye = reflectv.dot(eyev);

                // Compute the specular contribution
                let specular = if reflect_dot_eye > 0. {
                    let factor = reflect_dot_eye.powf(material.shininess);
                    light.intensity() * material.specular * factor
                } else {
                    Color::black()
                };

                (diffuse, specular)
            } else {
                (Color::black(), Color::black())
            };

            ambient + (diffuse + specular) * light_intensity
        }
        Light::AreaLight(light) => {
            let mut sum = Color::black();

            for u in 0..light.usteps() {
                for v in 0..light.vsteps() {
                    // Find the direction to the light source
                    let light_pos = light.point_on_light(u, v);
                    let lightv = (light_pos - point).normalize();

                    // `light_dot_normal` represents the cosine of the angle between the
                    // light vector and the normal vector. A negative number means the
                    // light is on the other side of the surface.
                    let light_dot_normal = lightv.dot(normalv);

                    if light_dot_normal >= 0. {
                        // Compute the diffuse contribution
                        let diffuse = effective_color * material.diffuse * light_dot_normal;

                        // `reflect_dot_eye` represents the cosine of the angle between
                        // the reflection vector and the eye vector. A negative value
                        // means that the lights reflects away from the eye.
                        let reflectv = (-lightv).reflect(normalv);
                        let reflect_dot_eye = reflectv.dot(eyev);

                        // Compute the specular contribution
                        let specular = if reflect_dot_eye > 0. {
                            let factor = reflect_dot_eye.powf(material.shininess);
                            light.intensity() * material.specular * factor
                        } else {
                            Color::black()
                        };

                        sum = sum + diffuse + specular;
                    };
                }
            }

            sum = sum * (1. / light.nb_samples());
            ambient + sum * light_intensity
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::lights::{AreaLight, PointLight};
    use crate::patterns::Striped;
    use crate::shapes::{Shape, Sphere};
    use crate::world::World;

    use super::*;

    #[test]
    fn lighting_basics() {
        let m = Material::default();
        let s = Shape::Sphere(Sphere::default());
        let pos = Point::new(0., 0., 0.);

        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::white());
        let res = lighting(m, s, light.into(), pos, eyev, normalv, 1.);
        assert!(res == Color::new(1.9, 1.9, 1.9));

        let eyev = Vector::new(0., f32::sqrt(2.) / 2., -f32::sqrt(2.) / 2.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::white());
        let res = lighting(m, s, light.into(), pos, eyev, normalv, 1.);
        assert!(res == Color::white());

        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 10., -10.), Color::white());
        let res = lighting(m, s, light.into(), pos, eyev, normalv, 1.);
        assert!(res.equal_approx(Color::new(0.7364, 0.7364, 0.7364)));

        let eyev = Vector::new(0., -f32::sqrt(2.) / 2., -f32::sqrt(2.) / 2.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 10., -10.), Color::white());
        let res = lighting(m, s, light.into(), pos, eyev, normalv, 1.);
        assert!(res.equal_approx(Color::new(1.6364, 1.6364, 1.6364)));

        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., 10.), Color::white());
        let res = lighting(m, s, light.into(), pos, eyev, normalv, 1.);
        assert!(res == Color::new(0.1, 0.1, 0.1));

        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::white());
        let res = lighting(m, s, light.into(), pos, eyev, normalv, 0.);
        assert!(res == Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn lighting_patterns() {
        let s = Shape::Sphere(Sphere::default());
        let mut m = Material::default();
        m.pattern = Striped::new(Color::white(), Color::black()).into();
        m.ambient = 1.;
        m.diffuse = 0.;
        m.specular = 0.;
        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::white()).into();
        let c1 = lighting(m, s, light, Point::new(0.9, 0., 0.), eyev, normalv, 1.);
        let c2 = lighting(m, s, light, Point::new(1.1, 0., 0.), eyev, normalv, 1.);
        assert!(c1 == Color::white());
        assert!(c2 == Color::black());
    }

    #[test]
    fn lighting_intensity() {
        let mut w = World::default();
        w.light = PointLight::new(Point::new(0., 0., -10.), Color::white()).into();
        w.shapes[0].material_mut().ambient = 0.1;
        w.shapes[0].material_mut().diffuse = 0.9;
        w.shapes[0].material_mut().specular = 0.;
        w.shapes[0].material_mut().color = Color::white();

        let pt = Point::new(0., 0., -1.);
        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);

        let cases = vec![
            (1., Color::new(1., 1., 1.)),
            (0.5, Color::new(0.55, 0.55, 0.55)),
            (0., Color::new(0.1, 0.1, 0.1)),
        ];
        for (intensity, res) in cases {
            assert!(
                lighting(
                    w.shapes[0].get_material(),
                    w.shapes[0],
                    w.light,
                    pt,
                    eyev,
                    normalv,
                    intensity
                ) == res
            )
        }
    }

    #[test]
    fn lighting_sample_area_light() {
        let s = Shape::Sphere(Sphere::default());
        let mut m = Material::default();
        m.pattern = Striped::new(Color::white(), Color::black()).into();
        m.ambient = 0.1;
        m.diffuse = 0.9;
        m.specular = 0.;
        m.color = Color::white();

        let light = AreaLight::new(
            Point::new(-0.5, -0.5, -5.),
            Vector::new(1., 0., 0.),
            2,
            Vector::new(0., 1., 0.),
            2,
            Color::white(),
            vec![0.5],
        );

        let eye = Point::new(0., 0., -5.);
        let cases = vec![
            (Point::new(0., 0., -1.), Color::new(0.9965, 0.9965, 0.9965)),
            (
                Point::new(0., 0.7071, -0.7071),
                Color::new(0.6232, 0.6232, 0.6232),
            ),
        ];
        for (pt, res) in cases {
            let eyev = -(eye - pt).normalize();
            let normalv = Vector::new(pt.x(), pt.y(), pt.z());
            let lighting = lighting(s.get_material(), s, light.into(), pt, eyev, normalv, 1.);
            assert!(lighting.equal_approx(res));
        }
    }
}
