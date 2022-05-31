use crate::intersections::Intersection;
use crate::matrices::Matrix;
use crate::rays::Ray;
use crate::tuples::Point;

/// A unit sphere centered on the origin
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sphere {
    transform: Matrix<4, 4>,
}

impl Sphere {
    pub fn new() -> Self {
        Self {
            transform: Matrix::identity(),
        }
    }

    pub fn intersect(&self, ray: Ray) -> Option<[Intersection; 2]> {
        let ray = match self.transform.inverse() {
            Some(inverse) => ray.transform(inverse),
            None => return None,
        };

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

    pub fn set_transform(&mut self, transform: Matrix<4, 4>) {
        self.transform = transform;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::intersections::Intersectable;
    use crate::transformations::*;
    use crate::tuples::Vector;

    #[test]
    fn sphere_ray_intersect() {
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.intersect(r);
        assert!(xs.unwrap()[0].object() == Intersectable::Sphere(s));
        assert!(xs.unwrap()[1].object() == Intersectable::Sphere(s));

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
        s.set_transform(scaling(2., 2., 2.));
        let xs = s.intersect(r);
        assert!(xs.unwrap()[0].t() == 3.);
        assert!(xs.unwrap()[1].t() == 7.);

        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let mut s = Sphere::new();
        s.set_transform(translation(5., 0., 0.));
        let xs = s.intersect(r);
        assert!(xs.is_none());
    }
}
