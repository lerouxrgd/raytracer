use derivative::Derivative;
use slotmap::Key;

use crate::groups::Group;
use crate::intersections::Intersection;
use crate::materials::Material;
use crate::matrices::Matrix;
use crate::rays::Ray;
use crate::shapes::Shape;
use crate::tuples::{Point, Vector};

/// A unit sphere centered on the origin
#[derive(Debug, Clone, Copy, Derivative)]
#[derivative(PartialEq)]
pub struct Sphere {
    pub(super) transform: Matrix<4, 4>,
    pub(super) material: Material,
    #[derivative(PartialEq = "ignore")]
    pub(super) parent: Group,
}

impl Sphere {
    pub fn with_transform<T: Into<Matrix<4, 4>>>(mut self, transform: T) -> Self {
        self.transform = transform.into();
        self
    }

    pub fn get_transform(&self) -> Matrix<4, 4> {
        self.transform
    }

    pub fn set_transform<T: Into<Matrix<4, 4>>>(&mut self, t: T) {
        self.transform = t.into()
    }

    pub fn with_material(mut self, material: Material) -> Self {
        self.material = material;
        self
    }

    pub fn get_material(&self) -> Material {
        self.material
    }

    pub fn set_material(&mut self, m: Material) {
        self.material = m
    }

    pub fn local_intersect(&self, local_ray: Ray) -> Option<[Intersection; 2]> {
        let sphere_to_ray = local_ray.origin - Point::new(0., 0., 0.);
        let a = local_ray.direction.dot(local_ray.direction);
        let b = 2. * local_ray.direction.dot(sphere_to_ray);
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

    pub fn local_normal_at(&self, local_point: Point) -> Vector {
        local_point - Point::new(0., 0., 0.)
    }
}

impl Default for Sphere {
    fn default() -> Self {
        Self {
            transform: Matrix::identity(),
            material: Material::default(),
            parent: Group::null(),
        }
    }
}

impl From<Sphere> for Shape {
    fn from(s: Sphere) -> Self {
        Self::Sphere(s)
    }
}

impl From<&Sphere> for Shape {
    fn from(s: &Sphere) -> Self {
        Self::Sphere(*s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transformations::*;
    use crate::tuples::Vector;
    use std::f32::consts::PI;

    #[test]
    fn sphere_ray_intersect() {
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::default();
        let xs = s.local_intersect(r);
        assert!(xs.unwrap()[0].shape() == Shape::Sphere(s));
        assert!(xs.unwrap()[1].shape() == Shape::Sphere(s));

        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::default();
        let xs = s.local_intersect(r);
        assert!(xs.unwrap()[0].t() == 4.);
        assert!(xs.unwrap()[1].t() == 6.);

        let r = Ray::new(Point::new(0., 1., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::default();
        let xs = s.local_intersect(r);
        assert!(xs.unwrap()[0].t() == 5.);
        assert!(xs.unwrap()[1].t() == 5.);

        let r = Ray::new(Point::new(0., 2., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::default();
        let xs = s.local_intersect(r);
        assert!(xs.is_none());

        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let s = Sphere::default();
        let xs = s.local_intersect(r);
        assert!(xs.unwrap()[0].t() == -1.);
        assert!(xs.unwrap()[1].t() == 1.);

        let r = Ray::new(Point::new(0., 0., 5.), Vector::new(0., 0., 1.));
        let s = Sphere::default();
        let xs = s.local_intersect(r);
        assert!(xs.unwrap()[0].t() == -6.);
        assert!(xs.unwrap()[1].t() == -4.);

        let s = Sphere::default().with_transform(scaling(2., 2., 2.));
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let mut xs = vec![];
        Shape::Sphere(s).intersect(&mut xs, r);
        assert!(xs[0].t() == 3.);
        assert!(xs[1].t() == 7.);

        let s = Sphere::default().with_transform(translation(5., 0., 0.));
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let mut xs = vec![];
        Shape::Sphere(s).intersect(&mut xs, r);
        assert!(xs.is_empty());
    }

    #[test]
    fn sphere_normal() {
        let s = Sphere::default();
        assert!(s.local_normal_at(Point::new(1., 0., 0.)) == Vector::new(1., 0., 0.));
        assert!(s.local_normal_at(Point::new(0., 1., 0.)) == Vector::new(0., 1., 0.));
        assert!(s.local_normal_at(Point::new(0., 0., 1.)) == Vector::new(0., 0., 1.));

        let s = Sphere::default();
        let p = Point::new(f32::sqrt(3.) / 3., f32::sqrt(3.) / 3., f32::sqrt(3.) / 3.);
        let n = Vector::new(f32::sqrt(3.) / 3., f32::sqrt(3.) / 3., f32::sqrt(3.) / 3.);
        assert!(s.local_normal_at(p).equal_approx(n));
        assert!(s
            .local_normal_at(p)
            .equal_approx(s.local_normal_at(p).normalize()));

        let s = Sphere::default().with_transform(translation(0., 1., 0.));
        let p = Point::new(0., 1.70711, -0.70711);
        let n = Shape::Sphere(s).normal_at(p);
        assert!(n.equal_approx(Vector::new(0., 0.70711, -0.70711)));

        let s = Sphere::default().with_transform(
            Transform::default()
                .rotation_z(PI / 5.)
                .scaling(1., 0.5, 1.),
        );
        let p = Point::new(0., f32::sqrt(2.) / 2., -f32::sqrt(2.) / 2.);
        let n = Shape::Sphere(s).normal_at(p);
        assert!(n.equal_approx(Vector::new(0., 0.97014, -0.24254)));
    }
}
