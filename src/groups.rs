use std::f32::consts::PI;
use std::mem;
use std::sync::RwLock;

use lazy_static::lazy_static;
use slotmap::{Key, SlotMap};

use crate::intersections::Intersection;
use crate::matrices::Matrix;
use crate::rays::Ray;
use crate::shapes::{Cylinder, Shape, Sphere};
use crate::transformations::*;
use crate::tuples::{Point, Tuple, Vector};

lazy_static! {
    static ref GROUPS: RwLock<SlotMap<Group, GroupData>> = RwLock::new(SlotMap::with_key());
}

struct GroupData {
    transform: Matrix<4, 4>,
    shapes: Vec<Shape>,
    children: Vec<Group>,
    parent: Group,
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
#[repr(transparent)]
pub struct Group(slotmap::KeyData);

impl From<slotmap::KeyData> for Group {
    fn from(k: slotmap::KeyData) -> Self {
        Self(k)
    }
}

unsafe impl slotmap::Key for Group {
    fn data(&self) -> slotmap::KeyData {
        self.0
    }
}

impl Group {
    pub fn with_transform<T: Into<Matrix<4, 4>>>(self, transform: T) -> Self {
        let mut groups = GROUPS.write().unwrap();
        let group = groups.get_mut(self).unwrap();
        group.transform = transform.into();
        self
    }

    pub fn get_transform(&self) -> Matrix<4, 4> {
        let groups = GROUPS.read().unwrap();
        let group = groups.get(*self).unwrap();
        group.transform
    }

    pub fn set_transform<T: Into<Matrix<4, 4>>>(&mut self, t: T) {
        let mut groups = GROUPS.write().unwrap();
        let group = groups.get_mut(*self).unwrap();
        group.transform = t.into();
    }

    pub fn add_child(&mut self, child: Group) {
        let mut groups = GROUPS.write().unwrap();

        let group = groups.get_mut(*self).unwrap();
        group.children.push(child);

        let child = groups.get_mut(child).unwrap();
        child.parent = *self;
    }

    pub fn get_child(&self, i: usize) -> Group {
        let groups = GROUPS.read().unwrap();
        let group = groups.get(*self).unwrap();
        group.children[i]
    }

    pub fn add_shape(&mut self, mut shape: Shape) {
        let mut groups = GROUPS.write().unwrap();
        let group = groups.get_mut(*self).unwrap();
        *shape.parent_mut() = *self;
        group.shapes.push(shape);
    }

    pub fn get_shape(&self, i: usize) -> Shape {
        let groups = GROUPS.read().unwrap();
        let group = groups.get(*self).unwrap();
        group.shapes[i]
    }

    pub fn set_shape(&mut self, i: usize, mut shape: Shape) {
        let mut groups = GROUPS.write().unwrap();
        let group = groups.get_mut(*self).unwrap();
        *shape.parent_mut() = *self;
        group.shapes[i] = shape;
    }

    pub fn len(&self) -> usize {
        let groups = GROUPS.read().unwrap();
        let group = groups.get(*self).unwrap();
        group.shapes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn local_intersect(&self, local_ray: Ray) -> Vec<Intersection> {
        let groups = GROUPS.read().unwrap();
        let group = groups.get(*self).unwrap();
        let mut xs = vec![];
        for shape in &group.shapes {
            shape.intersect(&mut xs, local_ray);
        }
        for child in &group.children {
            child.intersect(&mut xs, local_ray);
        }
        xs
    }

    pub fn intersect(&self, intersections: &mut Vec<Intersection>, world_ray: Ray) {
        let groups = GROUPS.read().unwrap();
        let group = groups.get(*self).unwrap();
        let inverse = group.transform.inverse().unwrap();
        let local_ray = world_ray.transform(inverse);
        intersections.extend(self.local_intersect(local_ray));
    }

    pub fn parent(&self) -> Group {
        let groups = GROUPS.read().unwrap();
        let group = groups.get(*self).unwrap();
        group.parent
    }

    pub fn world_to_object(&self, point: Point) -> Point {
        let parent = self.parent();
        let point = if !parent.is_null() {
            parent.world_to_object(point)
        } else {
            point
        };
        self.get_transform().inverse().unwrap() * point
    }

    pub fn normal_to_world(&self, normal: Vector) -> Vector {
        let groups = GROUPS.read().unwrap();
        let group = groups.get(*self).unwrap();

        let inverse = group.transform.inverse().unwrap();
        let normal = Tuple::from(inverse.transpose() * normal);
        let normal = Vector::new(normal[0], normal[1], normal[2]);
        let normal = normal.normalize();

        let parent = self.parent();
        if !parent.is_null() {
            parent.normal_to_world(normal)
        } else {
            normal
        }
    }

    pub fn delete(self) {
        let mut groups = GROUPS.write().unwrap();
        let group = groups.get_mut(self).unwrap();
        let children = mem::take(&mut group.children);
        for child in children {
            groups.remove(child);
        }
        groups.remove(self);
    }
}

impl Default for Group {
    fn default() -> Self {
        let group = GroupData {
            transform: Matrix::identity(),
            shapes: vec![],
            children: vec![],
            parent: Group::null(),
        };
        let handle = GROUPS.write().unwrap().insert(group);
        handle
    }
}

pub fn hexagon<T: Into<Matrix<4, 4>>>(transform: T) -> Group {
    fn hexagon_corner() -> Shape {
        Sphere::default()
            .with_transform(
                Transform::default()
                    .scaling(0.25, 0.25, 0.25)
                    .translation(0., 0., -1.),
            )
            .into()
    }

    fn hexagon_edge() -> Shape {
        Cylinder::default()
            .with_transform(
                Transform::default()
                    .scaling(0.25, 1., 0.25)
                    .rotation_z(-PI / 2.)
                    .rotation_y(-PI / 6.)
                    .translation(0., 0., -1.),
            )
            .min(0.)
            .max(1.)
            .into()
    }

    fn hexagon_side(rot_y: Matrix<4, 4>) -> Group {
        let mut side = Group::default().with_transform(rot_y);
        side.add_shape(hexagon_corner());
        side.add_shape(hexagon_edge());
        side
    }

    let mut hexagon = Group::default().with_transform(transform);
    for n in 0..6 {
        let side = hexagon_side(rotation_y(n as f32 * PI / 3.));
        hexagon.add_child(side);
    }

    hexagon
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::intersections::Intersections;
    use crate::tuples::{Point, Vector};

    #[test]
    fn group_intersections() {
        let g = Group::default();
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let xs = g.local_intersect(r);
        assert!(xs.is_empty());

        let mut g = Group::default();
        let s1 = Sphere::default();
        let s2 = Sphere::default().with_transform(translation(0., 0., -3.));
        let s3 = Sphere::default().with_transform(translation(5., 0., 0.));
        g.add_shape(s1.into());
        g.add_shape(s2.into());
        g.add_shape(s3.into());
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let xs = g.local_intersect(r);
        let xs = Vec::<Intersection>::from(Intersections::from(xs)); // sorted
        assert!(xs.len() == 4);
        assert!(xs[0].shape() == s2.into());
        assert!(xs[1].shape() == s2.into());
        assert!(xs[2].shape() == s1.into());
        assert!(xs[3].shape() == s1.into());

        let mut g = Group::default().with_transform(scaling(2., 2., 2.));
        let s = Sphere::default().with_transform(translation(5., 0., 0.));
        g.add_shape(s.into());
        let r = Ray::new(Point::new(10., 0., -10.), Vector::new(0., 0., 1.));
        let mut xs = vec![];
        g.intersect(&mut xs, r);
        assert!(xs.len() == 2);
    }

    #[test]
    fn group_normals() {
        let mut g1 = Group::default().with_transform(rotation_y(PI / 2.));
        let mut g2 = Group::default().with_transform(scaling(2., 2., 2.));
        let s = Sphere::default().with_transform(translation(5., 0., 0.));
        g2.add_shape(s.into());
        g1.add_child(g2);
        let p = Point::new(-2., 0., -10.);
        assert!(g2
            .get_shape(0)
            .world_to_object(p)
            .equal_approx(Point::new(0., 0., -1.)));

        let mut g1 = Group::default().with_transform(rotation_y(PI / 2.));
        let mut g2 = Group::default().with_transform(scaling(1., 2., 3.));
        let s = Sphere::default().with_transform(translation(5., 0., 0.));
        g2.add_shape(s.into());
        g1.add_child(g2);
        let n = Vector::new(f32::sqrt(3.) / 3., f32::sqrt(3.) / 3., f32::sqrt(3.) / 3.);
        assert!(g2
            .get_shape(0)
            .normal_to_world(n)
            .equal_approx(Vector::new(0.2857, 0.4286, -0.8571)));

        let mut g1 = Group::default().with_transform(rotation_y(PI / 2.));
        let mut g2 = Group::default().with_transform(scaling(1., 2., 3.));
        let s = Sphere::default().with_transform(translation(5., 0., 0.));
        g2.add_shape(s.into());
        g1.add_child(g2);
        let p = Point::new(1.7321, 1.1547, -5.5774);
        assert!(g2
            .get_shape(0)
            .normal_at(p, None)
            .equal_approx(Vector::new(0.2857, 0.4286, -0.8571)));
    }
}
