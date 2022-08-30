use crate::csg::Csg;
use crate::groups::Group;
use crate::intersections::{Computations, Intersections};
use crate::lights::{Light, PointLight};
use crate::materials::{lighting, Material};
use crate::rays::Ray;
use crate::shapes::{Shape, Sphere};
use crate::transformations::scaling;
use crate::tuples::{Color, Point};

#[derive(Debug, Clone, PartialEq)]
pub struct World {
    pub light: Light,
    pub shapes: Vec<Shape>,
    pub groups: Vec<Group>,
    pub csgs: Vec<Csg>,
}

impl Default for World {
    fn default() -> Self {
        let light = PointLight::new(Point::new(-10., 10., -10.), Color::white());

        let s1 = Sphere::default().with_material(
            Material::default()
                .color(Color::new(0.8, 1., 0.6))
                .diffuse(0.7)
                .specular(0.2),
        );

        let s2 = Sphere::default().with_transform(scaling(0.5, 0.5, 0.5));

        Self {
            light: light.into(),
            shapes: vec![s1.into(), s2.into()],
            groups: vec![],
            csgs: vec![],
        }
    }
}

impl World {
    pub fn intersect(&self, ray: Ray) -> Intersections {
        let xs = self
            .csgs
            .iter()
            .flat_map(|csg| csg.local_intersect(ray))
            .collect();

        let xs = self.groups.iter().fold(xs, |mut xs, group| {
            group.intersect(&mut xs, ray);
            xs
        });

        let xs = self.shapes.iter().fold(xs, |mut xs, shape| {
            shape.intersect(&mut xs, ray);
            xs
        });

        xs.into()
    }

    pub fn shade_hit(&self, comps: Computations, remaining: u8) -> Color {
        let light_intensity = self.light.intensity_at(comps.over_point, self);
        let surface = lighting(
            comps.shape.get_material(),
            comps.shape,
            self.light,
            comps.over_point,
            comps.eyev,
            comps.normalv,
            light_intensity,
        );

        let reflected = self.reflected_color(comps, remaining);
        let refracted = self.refracted_color(comps, remaining);

        let material = comps.shape.get_material();
        if material.reflective > 0. && material.transparency > 0. {
            let reflectance = comps.schlick();
            surface + reflected * reflectance + refracted * (1. - reflectance)
        } else {
            surface + reflected + refracted
        }
    }

    pub fn color_at(&self, ray: Ray, remaining: u8) -> Color {
        let xs = self.intersect(ray);
        if let Some(intersection) = xs.hit() {
            let comps = Computations::prepare(intersection, ray, &xs);
            self.shade_hit(comps, remaining)
        } else {
            Color::black()
        }
    }

    pub fn is_shadowed(&self, light_position: Point, p: Point) -> bool {
        let v = light_position - p;
        let distance = v.magnitude();
        let direction = v.normalize();
        let r = Ray::new(p, direction);
        if let Some(hit) = self.intersect(r).hit() {
            hit.t() < distance && hit.shape().shadow()
        } else {
            false
        }
    }

    pub fn reflected_color(&self, comps: Computations, remaining: u8) -> Color {
        if remaining == 0 {
            return Color::black();
        }

        if comps.shape.get_material().reflective == 0. {
            return Color::black();
        }

        let reflect_ray = Ray::new(comps.over_point, comps.reflectv);
        let color = self.color_at(reflect_ray, remaining - 1);
        color * comps.shape.get_material().reflective
    }

    pub fn refracted_color(&self, comps: Computations, remaining: u8) -> Color {
        if remaining == 0 {
            return Color::black();
        }

        if comps.shape.get_material().transparency == 0. {
            return Color::black();
        }

        // Check total internal reflection (Snell's Law)
        let n_ratio = comps.n1 / comps.n2;
        let cos_i = comps.eyev.dot(comps.normalv);
        let sin2_t = n_ratio.powi(2) * (1. - cos_i.powi(2));
        if sin2_t > 1. {
            return Color::black();
        }

        let cos_t = (1. - sin2_t).sqrt();
        let direction = (n_ratio * cos_i - cos_t) * comps.normalv - n_ratio * comps.eyev;
        let refracted_ray = Ray::new(comps.under_point, direction);
        self.color_at(refracted_ray, remaining - 1) * comps.shape.get_material().transparency
    }
}

impl Drop for World {
    fn drop(&mut self) {
        self.groups.drain(..).for_each(|group| group.delete());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::intersections::Intersection;
    use crate::patterns::XyzRgb;
    use crate::shapes::Plane;
    use crate::transformations::*;
    use crate::tuples::Vector;

    #[test]
    fn world_basics() {
        let world = World::default();
        let ray = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let xs = world.intersect(ray);
        assert!(xs.count() == 4);
        assert!(xs[0].t() == 4.);
        assert!(xs[1].t() == 4.5);
        assert!(xs[2].t() == 5.5);
        assert!(xs[3].t() == 6.);
    }

    #[test]
    fn world_shading() {
        let w = World::default();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let shape = w.shapes[0];
        let i = Intersection::new(4.0, shape);
        let comps = Computations::prepare(i, r, &vec![].into());
        let c = w.shade_hit(comps, 1);
        assert!(c.equal_approx(Color::new(0.38066, 0.47583, 0.2855)));

        let mut w = World::default();
        w.light = PointLight::new(Point::new(0., 0.25, 0.), Color::white()).into();
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let shape = w.shapes[1];
        let i = Intersection::new(0.5, shape);
        let comps = Computations::prepare(i, r, &vec![].into());
        let c = w.shade_hit(comps, 1);
        assert!(c.equal_approx(Color::new(0.90498, 0.90498, 0.90498)));

        let w = World::default();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 1., 0.));
        assert!(w.color_at(r, 1) == Color::black());

        let w = World::default();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        assert!(w
            .color_at(r, 1)
            .equal_approx(Color::new(0.38066, 0.47583, 0.2855)));

        let mut w = World::default();
        w.shapes[0].material_mut().ambient = 1.0;
        w.shapes[1].material_mut().ambient = 1.0;
        let r = Ray::new(Point::new(0., 0., 0.75), Vector::new(0., 0., -1.));
        assert!(w
            .color_at(r, 1)
            .equal_approx(w.shapes[1].get_material().color));

        let mut w = World::default();
        w.light = PointLight::new(Point::new(0., 0., -10.), Color::white()).into();
        let s1 = Sphere::default();
        w.shapes.push(s1.into());
        let s2 = Sphere::default().with_transform(translation(0., 0., 10.));
        w.shapes.push(s2.into());
        let r = Ray::new(Point::new(0., 0., 5.), Vector::new(0., 0., 1.));
        let i = Intersection::new(4., s2.into());
        let comps = Computations::prepare(i, r, &vec![].into());
        let c = w.shade_hit(comps, 1);
        assert!(c == Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn world_shadows() {
        let cases = vec![
            (Point::new(-10., -10., 10.), false),
            (Point::new(10., 10., 10.), true),
            (Point::new(-20., -20., -20.), false),
            (Point::new(-5., -5., -5.), false),
        ];
        let light_position = Point::new(-10., -10., -10.);
        let w = World::default();
        for (point, res) in cases {
            assert!(w.is_shadowed(light_position, point) == res);
        }
    }

    #[test]
    fn world_reflection() {
        let mut w = World::default();
        w.shapes[1].material_mut().ambient = 1.;
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let i = Intersection::new(1., w.shapes[1]);
        let comps = Computations::prepare(i, r, &vec![].into());
        assert!(w.reflected_color(comps, 1) == Color::black());

        let mut w = World::default();
        let shape = Plane::default()
            .with_transform(translation(0., -1., 0.))
            .with_material(Material::default().reflective(0.5));
        w.shapes.push(shape.into());
        let r = Ray::new(
            Point::new(0., 0., -3.),
            Vector::new(0., -f32::sqrt(2.) / 2., f32::sqrt(2.) / 2.),
        );
        let i = Intersection::new(f32::sqrt(2.), shape.into());
        let comps = Computations::prepare(i, r, &vec![].into());
        assert!(w
            .reflected_color(comps, 1)
            .equal_approx(Color::new(0.19032, 0.2379, 0.14274)));

        let mut w = World::default();
        let shape = Plane::default()
            .with_transform(translation(0., -1., 0.))
            .with_material(Material::default().reflective(0.5));
        w.shapes.push(shape.into());
        let r = Ray::new(
            Point::new(0., 0., -3.),
            Vector::new(0., -f32::sqrt(2.) / 2., f32::sqrt(2.) / 2.),
        );
        let i = Intersection::new(f32::sqrt(2.), shape.into());
        let comps = Computations::prepare(i, r, &vec![].into());
        assert!(w
            .shade_hit(comps, 1)
            .equal_approx(Color::new(0.87677, 0.92436, 0.82918)));

        let mut w = World::default();
        w.light = PointLight::new(Point::new(0., 0., 0.), Color::white()).into();
        let lower = Plane::default()
            .with_transform(translation(0., -1., 0.))
            .with_material(Material::default().reflective(1.));
        w.shapes.push(lower.into());
        let upper = Plane::default()
            .with_transform(translation(0., 1., 0.))
            .with_material(Material::default().reflective(1.));
        w.shapes.push(upper.into());
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 1., 0.));
        w.color_at(r, 1);
    }

    #[test]
    fn world_refraction() {
        let w = World::default();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let xs: Intersections = vec![
            Intersection::new(4., w.shapes[1]),
            Intersection::new(6., w.shapes[1]),
        ]
        .into();
        let comps = Computations::prepare(xs[0], r, &xs);
        assert!(w.refracted_color(comps, 1) == Color::black());

        let mut w = World::default();
        w.shapes[0].material_mut().transparency = 1.;
        w.shapes[0].material_mut().refractive_index = 1.5;
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let xs: Intersections = vec![
            Intersection::new(4., w.shapes[1]),
            Intersection::new(6., w.shapes[1]),
        ]
        .into();
        let comps = Computations::prepare(xs[0], r, &xs);
        assert!(w.refracted_color(comps, 0) == Color::black());

        let mut w = World::default();
        w.shapes[0].material_mut().transparency = 1.;
        w.shapes[0].material_mut().refractive_index = 1.5;
        let r = Ray::new(
            Point::new(0., 0., f32::sqrt(2.) / 2.),
            Vector::new(0., 1., 0.),
        );
        let xs: Intersections = vec![
            Intersection::new(-f32::sqrt(2.) / 2., w.shapes[1]),
            Intersection::new(f32::sqrt(2.) / 2., w.shapes[1]),
        ]
        .into();
        let comps = Computations::prepare(xs[1], r, &xs);
        assert!(w.refracted_color(comps, 1) == Color::black());

        let mut w = World::default();
        w.shapes[0].material_mut().ambient = 1.;
        w.shapes[0].material_mut().pattern = XyzRgb::new().into();
        w.shapes[1].material_mut().transparency = 1.;
        w.shapes[1].material_mut().refractive_index = 1.5;
        let r = Ray::new(Point::new(0., 0., 0.1), Vector::new(0., 1., 0.));
        let xs: Intersections = vec![
            Intersection::new(-0.9899, w.shapes[0]),
            Intersection::new(-0.4899, w.shapes[1]),
            Intersection::new(0.4899, w.shapes[1]),
            Intersection::new(0.9899, w.shapes[0]),
        ]
        .into();
        let comps = Computations::prepare(xs[2], r, &xs);
        assert!(w
            .refracted_color(comps, 1)
            .equal_approx(Color::new(0., 0.99888, 0.04725)));

        let mut w = World::default();
        let floor = Plane::default()
            .with_transform(translation(0., -1., 0.))
            .with_material(Material::default().transparency(0.5).refractive_index(1.5));
        w.shapes.push(floor.into());
        let ball = Sphere::default()
            .with_transform(translation(0., -3.5, -0.5))
            .with_material(
                Material::default()
                    .color(Color::new(1., 0., 0.))
                    .ambient(0.5),
            );
        w.shapes.push(ball.into());
        let r = Ray::new(
            Point::new(0., 0., -3.),
            Vector::new(0., -f32::sqrt(2.) / 2., f32::sqrt(2.) / 2.),
        );
        let xs: Intersections = vec![Intersection::new(f32::sqrt(2.), floor.into())].into();
        let comps = Computations::prepare(xs[0], r, &xs);
        assert!(w
            .shade_hit(comps, 1)
            .equal_approx(Color::new(0.93642, 0.68642, 0.68642)));

        let mut w = World::default();
        let floor = Plane::default()
            .with_transform(translation(0., -1., 0.))
            .with_material(
                Material::default()
                    .reflective(0.5)
                    .transparency(0.5)
                    .refractive_index(1.5),
            );
        w.shapes.push(floor.into());
        let ball = Sphere::default()
            .with_transform(translation(0., -3.5, -0.5))
            .with_material(
                Material::default()
                    .color(Color::new(1., 0., 0.))
                    .ambient(0.5),
            );
        w.shapes.push(ball.into());
        let r = Ray::new(
            Point::new(0., 0., -3.),
            Vector::new(0., -f32::sqrt(2.) / 2., f32::sqrt(2.) / 2.),
        );
        let xs: Intersections = vec![Intersection::new(f32::sqrt(2.), floor.into())].into();
        let comps = Computations::prepare(xs[0], r, &xs);
        assert!(w
            .shade_hit(comps, 1)
            .equal_approx(Color::new(0.93391, 0.69643, 0.69243)));
    }
}
