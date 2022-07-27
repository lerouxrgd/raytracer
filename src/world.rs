use crate::groups::Group;
use crate::intersections::{Computations, Intersections};
use crate::lights::PointLight;
use crate::materials::lighting;
use crate::rays::Ray;
use crate::shapes::{Shape, Sphere};
use crate::transformations::scaling;
use crate::tuples::{Color, Point};

#[derive(Debug, Clone, PartialEq)]
pub struct World {
    pub light: PointLight,
    pub objects: Vec<Shape>,
    pub groups: Vec<Group>,
}

impl Default for World {
    fn default() -> Self {
        let light = PointLight::new(Point::new(-10., 10., -10.), Color::white());

        let mut s1 = Sphere::new();
        s1.material.color = Color::new(0.8, 1., 0.6);
        s1.material.diffuse = 0.7;
        s1.material.specular = 0.2;

        let mut s2 = Sphere::new();
        s2.transform = scaling(0.5, 0.5, 0.5);

        Self {
            light,
            objects: vec![s1.into(), s2.into()],
            groups: vec![],
        }
    }
}

impl World {
    pub fn intersect(&self, ray: Ray) -> Intersections {
        let xs = self.groups.iter().fold(vec![], |mut xs, group| {
            group.intersect(&mut xs, ray);
            xs
        });

        self.objects
            .iter()
            .fold(xs, |mut xs, shape| {
                shape.intersect(&mut xs, ray);
                xs
            })
            .into()
    }

    pub fn shade_hit(&self, comps: Computations, remaining: u8) -> Color {
        let shadowed = self.is_shadowed(comps.over_point);
        let surface = lighting(
            comps.object.material(),
            comps.object,
            self.light,
            comps.over_point,
            comps.eyev,
            comps.normalv,
            shadowed,
        );

        let reflected = self.reflected_color(comps, remaining);
        let refracted = self.refracted_color(comps, remaining);

        let material = comps.object.material();
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

    pub fn is_shadowed(&self, p: Point) -> bool {
        let v = self.light.position - p;
        let distance = v.magnitude();
        let direction = v.normalize();
        let r = Ray::new(p, direction);
        if let Some(hit) = self.intersect(r).hit() {
            hit.t() < distance
        } else {
            false
        }
    }

    pub fn reflected_color(&self, comps: Computations, remaining: u8) -> Color {
        if remaining == 0 {
            return Color::black();
        }

        if comps.object.material().reflective == 0. {
            return Color::black();
        }

        let reflect_ray = Ray::new(comps.over_point, comps.reflectv);
        let color = self.color_at(reflect_ray, remaining - 1);
        color * comps.object.material().reflective
    }

    pub fn refracted_color(&self, comps: Computations, remaining: u8) -> Color {
        if remaining == 0 {
            return Color::black();
        }

        if comps.object.material().transparency == 0. {
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
        self.color_at(refracted_ray, remaining - 1) * comps.object.material().transparency
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
        let shape = w.objects[0];
        let i = Intersection::new(4.0, shape);
        let comps = Computations::prepare(i, r, &vec![].into());
        let c = w.shade_hit(comps, 1);
        assert!(c.equal_approx(Color::new(0.38066, 0.47583, 0.2855)));

        let mut w = World::default();
        w.light = PointLight::new(Point::new(0., 0.25, 0.), Color::white());
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let shape = w.objects[1];
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
        w.objects[0].material_mut().ambient = 1.0;
        w.objects[1].material_mut().ambient = 1.0;
        let r = Ray::new(Point::new(0., 0., 0.75), Vector::new(0., 0., -1.));
        assert!(w.color_at(r, 1).equal_approx(w.objects[1].material().color));

        let mut w = World::default();
        w.light = PointLight::new(Point::new(0., 0., -10.), Color::white());
        let s1 = Sphere::new();
        w.objects.push(s1.into());
        let mut s2 = Sphere::new();
        s2.transform = translation(0., 0., 10.);
        w.objects.push(s2.into());
        let r = Ray::new(Point::new(0., 0., 5.), Vector::new(0., 0., 1.));
        let i = Intersection::new(4., s2.into());
        let comps = Computations::prepare(i, r, &vec![].into());
        let c = w.shade_hit(comps, 1);
        assert!(c == Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn world_shadows() {
        let w = World::default();
        let p = Point::new(0., 10., 0.);
        assert!(w.is_shadowed(p) == false);

        let w = World::default();
        let p = Point::new(10., -10., 10.);
        assert!(w.is_shadowed(p) == true);

        let w = World::default();
        let p = Point::new(-20., 20., -20.);
        assert!(w.is_shadowed(p) == false);

        let w = World::default();
        let p = Point::new(-2., 2., -2.);
        assert!(w.is_shadowed(p) == false);
    }

    #[test]
    fn world_reflection() {
        let mut w = World::default();
        w.objects[1].material_mut().ambient = 1.;
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let i = Intersection::new(1., w.objects[1]);
        let comps = Computations::prepare(i, r, &vec![].into());
        assert!(w.reflected_color(comps, 1) == Color::black());

        let mut w = World::default();
        let mut shape = Plane::new();
        shape.material.reflective = 0.5;
        shape.transform = translation(0., -1., 0.);
        w.objects.push(shape.into());
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
        let mut shape = Plane::new();
        shape.material.reflective = 0.5;
        shape.transform = translation(0., -1., 0.);
        w.objects.push(shape.into());
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
        w.light = PointLight::new(Point::new(0., 0., 0.), Color::white());
        let mut lower = Plane::new();
        lower.material.reflective = 1.;
        lower.transform = translation(0., -1., 0.);
        w.objects.push(lower.into());
        let mut upper = Plane::new();
        upper.material.reflective = 1.;
        upper.transform = translation(0., 1., 0.);
        w.objects.push(upper.into());
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 1., 0.));
        w.color_at(r, 1);
    }

    #[test]
    fn world_refraction() {
        let w = World::default();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let xs: Intersections = vec![
            Intersection::new(4., w.objects[1]),
            Intersection::new(6., w.objects[1]),
        ]
        .into();
        let comps = Computations::prepare(xs[0], r, &xs);
        assert!(w.refracted_color(comps, 1) == Color::black());

        let mut w = World::default();
        w.objects[0].material_mut().transparency = 1.;
        w.objects[0].material_mut().refractive_index = 1.5;
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let xs: Intersections = vec![
            Intersection::new(4., w.objects[1]),
            Intersection::new(6., w.objects[1]),
        ]
        .into();
        let comps = Computations::prepare(xs[0], r, &xs);
        assert!(w.refracted_color(comps, 0) == Color::black());

        let mut w = World::default();
        w.objects[0].material_mut().transparency = 1.;
        w.objects[0].material_mut().refractive_index = 1.5;
        let r = Ray::new(
            Point::new(0., 0., f32::sqrt(2.) / 2.),
            Vector::new(0., 1., 0.),
        );
        let xs: Intersections = vec![
            Intersection::new(-f32::sqrt(2.) / 2., w.objects[1]),
            Intersection::new(f32::sqrt(2.) / 2., w.objects[1]),
        ]
        .into();
        let comps = Computations::prepare(xs[1], r, &xs);
        assert!(w.refracted_color(comps, 1) == Color::black());

        let mut w = World::default();
        w.objects[0].material_mut().ambient = 1.;
        w.objects[0].material_mut().pattern = XyzRgb::new().into();
        w.objects[1].material_mut().transparency = 1.;
        w.objects[1].material_mut().refractive_index = 1.5;
        let r = Ray::new(Point::new(0., 0., 0.1), Vector::new(0., 1., 0.));
        let xs: Intersections = vec![
            Intersection::new(-0.9899, w.objects[0]),
            Intersection::new(-0.4899, w.objects[1]),
            Intersection::new(0.4899, w.objects[1]),
            Intersection::new(0.9899, w.objects[0]),
        ]
        .into();
        let comps = Computations::prepare(xs[2], r, &xs);
        assert!(w
            .refracted_color(comps, 1)
            .equal_approx(Color::new(0., 0.99888, 0.04725)));

        let mut w = World::default();
        let mut floor = Plane::new();
        floor.transform = translation(0., -1., 0.);
        floor.material.transparency = 0.5;
        floor.material.refractive_index = 1.5;
        w.objects.push(floor.into());
        let mut ball = Sphere::new();
        ball.material.color = Color::new(1., 0., 0.);
        ball.material.ambient = 0.5;
        ball.transform = translation(0., -3.5, -0.5);
        w.objects.push(ball.into());
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
        let mut floor = Plane::new();
        floor.transform = translation(0., -1., 0.);
        floor.material.reflective = 0.5;
        floor.material.transparency = 0.5;
        floor.material.refractive_index = 1.5;
        w.objects.push(floor.into());
        let mut ball = Sphere::new();
        ball.material.color = Color::new(1., 0., 0.);
        ball.material.ambient = 0.5;
        ball.transform = translation(0., -3.5, -0.5);
        w.objects.push(ball.into());
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
