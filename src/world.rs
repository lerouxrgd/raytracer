use std::ops::{Index, IndexMut};

use crate::intersections::{Computations, Intersections, Shape};
use crate::lights::LightPoint;
use crate::materials::lighting;
use crate::rays::Ray;
use crate::spheres::Sphere;
use crate::transformations::scaling;
use crate::tuples::{Color, Point};

#[derive(Debug, Clone, PartialEq)]
pub struct World {
    pub light: LightPoint,
    pub objects: Vec<Shape>,
}

impl Default for World {
    fn default() -> Self {
        let light = LightPoint::new(Point::new(-10., 10., -10.), Color::white());

        let mut s1 = Sphere::new();
        s1.material.color = Color::new(0.8, 1., 0.6);
        s1.material.diffuse = 0.7;
        s1.material.specular = 0.2;

        let mut s2 = Sphere::new();
        s2.transform = scaling(0.5, 0.5, 0.5);

        Self {
            light,
            objects: vec![s1.into(), s2.into()],
        }
    }
}

impl World {
    pub fn intersect(&self, ray: Ray) -> Intersections {
        let mut intersections = vec![];
        for object in self.objects.iter() {
            match object {
                Shape::Sphere(s) => {
                    if let Some(xs) = s.intersect(ray) {
                        intersections.extend(xs)
                    }
                }
            }
        }
        intersections.into()
    }

    pub fn shade_hit(&self, comps: Computations) -> Color {
        let shadowed = self.is_shadowed(comps.over_point);
        lighting(
            comps.object.material(),
            self.light,
            comps.point,
            comps.eyev,
            comps.normalv,
            shadowed,
        )
    }

    pub fn color_at(&self, ray: Ray) -> Color {
        let xs = self.intersect(ray);
        if let Some(intersection) = xs.hit() {
            let comps = Computations::prepare(intersection, ray);
            self.shade_hit(comps)
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
}

impl Index<usize> for World {
    type Output = Shape;

    fn index(&self, i: usize) -> &Self::Output {
        &self.objects[i]
    }
}

impl IndexMut<usize> for World {
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        &mut self.objects[i]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::intersections::Intersection;
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
        let shape = w[0];
        let i = Intersection::new(4.0, shape);
        let comps = Computations::prepare(i, r);
        let c = w.shade_hit(comps);
        assert!(c.equal_approx(Color::new(0.38066, 0.47583, 0.2855)));

        let mut w = World::default();
        w.light = LightPoint::new(Point::new(0., 0.25, 0.), Color::white());
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let shape = w[1];
        let i = Intersection::new(0.5, shape);
        let comps = Computations::prepare(i, r);
        let c = w.shade_hit(comps);
        assert!(c.equal_approx(Color::new(0.90498, 0.90498, 0.90498)));

        let w = World::default();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 1., 0.));
        assert!(w.color_at(r) == Color::black());

        let w = World::default();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        assert!(w
            .color_at(r)
            .equal_approx(Color::new(0.38066, 0.47583, 0.2855)));

        let mut w = World::default();
        w[0].material_mut().ambient = 1.0;
        w[1].material_mut().ambient = 1.0;
        let r = Ray::new(Point::new(0., 0., 0.75), Vector::new(0., 0., -1.));
        assert!(w.color_at(r).equal_approx(w[1].material().color));

        let mut w = World::default();
        w.light = LightPoint::new(Point::new(0., 0., -10.), Color::white());
        let s1 = Sphere::new();
        w.objects.push(s1.into());
        let mut s2 = Sphere::new();
        s2.transform = translation(0., 0., 10.);
        w.objects.push(s2.into());
        let r = Ray::new(Point::new(0., 0., 5.), Vector::new(0., 0., 1.));
        let i = Intersection::new(4., s2.into());
        let comps = Computations::prepare(i, r);
        let c = w.shade_hit(comps);
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
}
