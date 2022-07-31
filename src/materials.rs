use crate::lights::PointLight;
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
    light: PointLight,
    point: Point,
    eyev: Vector,
    normalv: Vector,
    in_shadow: bool,
) -> Color {
    let color = match material.pattern {
        Some(pattern) => pattern.pattern_at_shape(shape, point),
        None => material.color,
    };

    // Combine the surface color with the light's color/intensity
    let effective_color = color * light.intensity;

    // Compute the ambient contribution
    let ambient = effective_color * material.ambient;

    // Find the direction to the light source
    let lightv = (light.position - point).normalize();

    // light_dot_normal represents the cosine of the angle between the light vector and
    // the normal vector. A negative number means the light is on the other side of the surface.
    let light_dot_normal = lightv.dot(normalv);
    let (diffuse, specular) = if light_dot_normal < 0. {
        let diffuse = Color::black();
        let specular = Color::black();
        (diffuse, specular)
    } else {
        // Compute the diffuse contribution
        let diffuse = effective_color * material.diffuse * light_dot_normal;

        // reflect_dot_eye represents the cosine of the angle between the reflection
        // vector and the eye vector. A negative value means that the lights reflects
        // away from the eye.
        let reflectv = (-lightv).reflect(normalv);
        let reflect_dot_eye = reflectv.dot(eyev);
        let specular = if reflect_dot_eye < 0. {
            Color::black()
        } else {
            // Compute the specular contribution
            let factor = reflect_dot_eye.powf(material.shininess);
            light.intensity * material.specular * factor
        };
        (diffuse, specular)
    };

    if in_shadow {
        ambient
    } else {
        ambient + diffuse + specular
    }
}

#[cfg(test)]
mod tests {
    use crate::patterns::Striped;
    use crate::shapes::{Shape, Sphere};

    use super::*;

    #[test]
    fn lighting_basics() {
        let m = Material::default();
        let s = Shape::Sphere(Sphere::default());
        let pos = Point::new(0., 0., 0.);

        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::white());
        let res = lighting(m, s, light, pos, eyev, normalv, false);
        assert!(res == Color::new(1.9, 1.9, 1.9));

        let eyev = Vector::new(0., f32::sqrt(2.) / 2., -f32::sqrt(2.) / 2.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::white());
        let res = lighting(m, s, light, pos, eyev, normalv, false);
        assert!(res == Color::white());

        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 10., -10.), Color::white());
        let res = lighting(m, s, light, pos, eyev, normalv, false);
        assert!(res.equal_approx(Color::new(0.7364, 0.7364, 0.7364)));

        let eyev = Vector::new(0., -f32::sqrt(2.) / 2., -f32::sqrt(2.) / 2.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 10., -10.), Color::white());
        let res = lighting(m, s, light, pos, eyev, normalv, false);
        assert!(res.equal_approx(Color::new(1.6364, 1.6364, 1.6364)));

        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., 10.), Color::white());
        let res = lighting(m, s, light, pos, eyev, normalv, false);
        assert!(res == Color::new(0.1, 0.1, 0.1));

        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::white());
        let res = lighting(m, s, light, pos, eyev, normalv, true);
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
        let light = PointLight::new(Point::new(0., 0., -10.), Color::white());
        let c1 = lighting(m, s, light, Point::new(0.9, 0., 0.), eyev, normalv, false);
        let c2 = lighting(m, s, light, Point::new(1.1, 0., 0.), eyev, normalv, false);
        assert!(c1 == Color::white());
        assert!(c2 == Color::black());
    }
}
