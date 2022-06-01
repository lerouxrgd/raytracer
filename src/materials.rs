use crate::lights::LightPoint;
use crate::tuples::{Color, Point, Vector};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Material {
    pub color: Color,
    pub ambient: f32,
    pub diffuse: f32,
    pub specular: f32,
    pub shininess: f32,
}

impl Default for Material {
    fn default() -> Self {
        Material {
            color: Color::white(),
            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.9,
            shininess: 200.,
        }
    }
}

pub fn lighting(
    material: Material,
    light: LightPoint,
    point: Point,
    eyev: Vector,
    normalv: Vector,
) -> Color {
    // Combine the surface color with the light's color/intensity
    let effective_color = material.color * light.intensity;

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

    ambient + diffuse + specular
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lighting_basics() {
        let m = Material::default();
        let pos = Point::new(0., 0., 0.);

        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = LightPoint::new(Point::new(0., 0., -10.), Color::white());
        let res = lighting(m, light, pos, eyev, normalv);
        assert!(res == Color::new(1.9, 1.9, 1.9));

        let eyev = Vector::new(0., f32::sqrt(2.) / 2., -f32::sqrt(2.) / 2.);
        let normalv = Vector::new(0., 0., -1.);
        let light = LightPoint::new(Point::new(0., 0., -10.), Color::white());
        let res = lighting(m, light, pos, eyev, normalv);
        assert!(res == Color::white());

        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = LightPoint::new(Point::new(0., 10., -10.), Color::white());
        let res = lighting(m, light, pos, eyev, normalv);
        assert!(res.equal_approx(Color::new(0.7364, 0.7364, 0.7364)));

        let eyev = Vector::new(0., -f32::sqrt(2.) / 2., -f32::sqrt(2.) / 2.);
        let normalv = Vector::new(0., 0., -1.);
        let light = LightPoint::new(Point::new(0., 10., -10.), Color::white());
        let res = lighting(m, light, pos, eyev, normalv);
        assert!(res.equal_approx(Color::new(1.6364, 1.6364, 1.6364)));

        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = LightPoint::new(Point::new(0., 0., 10.), Color::white());
        let res = lighting(m, light, pos, eyev, normalv);
        assert!(res == Color::new(0.1, 0.1, 0.1));
    }
}
