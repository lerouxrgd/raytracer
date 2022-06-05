use crate::intersections::Intersection;
use crate::materials::Material;
use crate::matrices::Matrix;
use crate::rays::Ray;
use crate::tuples::{Point, Tuple, Vector};

/// A unit sphere centered on the origin
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sphere {
    pub transform: Matrix<4, 4>,
    pub material: Material,
}

impl Sphere {
    pub fn new() -> Self {
        Self {
            transform: Matrix::identity(),
            material: Material::default(),
        }
    }

    pub fn intersect(&self, ray: Ray) -> Option<[Intersection; 2]> {
        let ray = ray.transform(self.transform.inverse().unwrap());

        let sphere_to_ray = ray.origin - Point::new(0., 0., 0.);
        let a = ray.direction.dot(ray.direction);
        let b = 2. * ray.direction.dot(sphere_to_ray);
        let c = sphere_to_ray.dot(sphere_to_ray) - 1.;
        let discriminant = b.powi(2) - 4. * a * c;

        if discriminant < 0. {
            None
        } else {
            let t1 = (-b - discriminant.sqrt()) / (2. * a);
            let t2 = (-b + discriminant.sqrt()) / (2. * a);
            Some([
                Intersection::new(t1, self.into()),
                Intersection::new(t2, self.into()),
            ])
        }
    }

    pub fn normal_at(&self, world_point: Point) -> Vector {
        let inverse = self.transform.inverse().unwrap();
        let object_point = inverse * world_point;
        let object_normal = object_point - Point::new(0., 0., 0.);
        let word_normal = inverse.transpose() * object_normal;
        let word_normal = Tuple::from(word_normal);
        let word_normal = Vector::new(word_normal[0], word_normal[1], word_normal[2]);
        word_normal.normalize()
    }
}

impl Default for Sphere {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::intersections::Shape;
    use crate::transformations::*;
    use crate::tuples::Vector;
    use std::f32::consts::PI;

    #[test]
    fn sphere_ray_intersect() {
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.intersect(r);
        assert!(xs.unwrap()[0].object() == Shape::Sphere(s));
        assert!(xs.unwrap()[1].object() == Shape::Sphere(s));

        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.intersect(r);
        assert!(xs.unwrap()[0].t() == 4.);
        assert!(xs.unwrap()[1].t() == 6.);

        let r = Ray::new(Point::new(0., 1., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.intersect(r);
        assert!(xs.unwrap()[0].t() == 5.);
        assert!(xs.unwrap()[1].t() == 5.);

        let r = Ray::new(Point::new(0., 2., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.intersect(r);
        assert!(xs.is_none());

        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.intersect(r);
        assert!(xs.unwrap()[0].t() == -1.);
        assert!(xs.unwrap()[1].t() == 1.);

        let r = Ray::new(Point::new(0., 0., 5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.intersect(r);
        assert!(xs.unwrap()[0].t() == -6.);
        assert!(xs.unwrap()[1].t() == -4.);

        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let mut s = Sphere::new();
        s.transform = scaling(2., 2., 2.);
        let xs = s.intersect(r);
        assert!(xs.unwrap()[0].t() == 3.);
        assert!(xs.unwrap()[1].t() == 7.);

        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let mut s = Sphere::new();
        s.transform = translation(5., 0., 0.);
        let xs = s.intersect(r);
        assert!(xs.is_none());
    }

    #[test]
    fn sphere_normal() {
        let s = Sphere::new();
        assert!(s.normal_at(Point::new(1., 0., 0.)) == Vector::new(1., 0., 0.));
        assert!(s.normal_at(Point::new(0., 1., 0.)) == Vector::new(0., 1., 0.));
        assert!(s.normal_at(Point::new(0., 0., 1.)) == Vector::new(0., 0., 1.));

        let s = Sphere::new();
        let p = Point::new(f32::sqrt(3.) / 3., f32::sqrt(3.) / 3., f32::sqrt(3.) / 3.);
        let n = Vector::new(f32::sqrt(3.) / 3., f32::sqrt(3.) / 3., f32::sqrt(3.) / 3.);
        assert!(s.normal_at(p).equal_approx(n));
        assert!(s.normal_at(p).equal_approx(s.normal_at(p).normalize()));

        let mut s = Sphere::new();
        s.transform = translation(0., 1., 0.);
        let n = s.normal_at(Point::new(0., 1.70711, -0.70711));
        assert!(n.equal_approx(Vector::new(0., 0.70711, -0.70711)));

        let mut s = Sphere::new();
        s.transform = Transform::new()
            .roation_z(PI / 5.)
            .scaling(1., 0.5, 1.)
            .into();
        let n = s.normal_at(Point::new(0., f32::sqrt(2.) / 2., -f32::sqrt(2.) / 2.));
        assert!(n.equal_approx(Vector::new(0., 0.97014, -0.24254)));
    }

    #[test]
    fn sphere_material() {
        let s = Sphere::new();
        let m = Material::default();
        assert!(s.material == m);
    }
}
