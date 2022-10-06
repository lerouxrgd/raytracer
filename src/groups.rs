use std::f32::consts::PI;
use std::mem;

use lazy_static::lazy_static;
use parking_lot::RwLock;
use slotmap::{Key, SlotMap};

use crate::bounds::BoundingBox;
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
    bounds: Option<BoundingBox>,
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
        let mut groups = GROUPS.write();
        let group = groups.get_mut(self).unwrap();
        group.transform = transform.into();
        self
    }

    pub fn get_transform(&self) -> Matrix<4, 4> {
        let groups = GROUPS.read_recursive();
        let group = groups.get(*self).unwrap();
        group.transform
    }

    pub fn set_transform<T: Into<Matrix<4, 4>>>(&mut self, t: T) {
        let mut groups = GROUPS.write();
        let group = groups.get_mut(*self).unwrap();
        group.transform = t.into();
    }

    pub fn add_child(&mut self, child: Group) {
        let mut groups = GROUPS.write();

        let group = groups.get_mut(*self).unwrap();
        group.children.push(child);

        let child = groups.get_mut(child).unwrap();
        child.parent = *self;
    }

    pub fn get_child(&self, i: usize) -> Group {
        let groups = GROUPS.read_recursive();
        let group = groups.get(*self).unwrap();
        group.children[i]
    }

    pub fn add_shape(&mut self, mut shape: Shape) {
        let mut groups = GROUPS.write();
        let group = groups.get_mut(*self).unwrap();
        *shape.parent_mut() = *self;
        group.shapes.push(shape);
    }

    pub fn get_shape(&self, i: usize) -> Shape {
        let groups = GROUPS.read_recursive();
        let group = groups.get(*self).unwrap();
        group.shapes[i]
    }

    pub fn set_shape(&mut self, i: usize, mut shape: Shape) {
        let mut groups = GROUPS.write();
        let group = groups.get_mut(*self).unwrap();
        *shape.parent_mut() = *self;
        group.shapes[i] = shape;
    }

    pub fn len(&self) -> usize {
        let groups = GROUPS.read_recursive();
        let group = groups.get(*self).unwrap();
        group.shapes.len() + group.children.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn local_intersect(&self, local_ray: Ray) -> Vec<Intersection> {
        if !self.bounds().intersects(local_ray) {
            return vec![];
        }

        let groups = GROUPS.read_recursive();
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
        let groups = GROUPS.read_recursive();
        let group = groups.get(*self).unwrap();
        let inverse = group.transform.inverse().unwrap();
        let local_ray = world_ray.transform(inverse);
        intersections.extend(self.local_intersect(local_ray));
    }

    pub fn parent(&self) -> Group {
        GROUPS.read_recursive().get(*self).unwrap().parent
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
        let groups = GROUPS.read_recursive();
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

    pub fn bounds(&self) -> BoundingBox {
        let groups = GROUPS.read_recursive();
        let group = groups.get(*self).unwrap();
        if let Some(bb) = group.bounds {
            bb
        } else {
            Self::make_bounds(group)
        }
    }

    fn make_bounds(group: &GroupData) -> BoundingBox {
        let mut bb = BoundingBox::default();
        for shape in &group.shapes {
            bb.add_box(shape.parent_space_bounds());
        }
        for child in &group.children {
            if group.parent.is_null() {
                bb = bb.transform(group.transform);
            }
            bb.add_box(child.bounds());
        }
        bb
    }

    pub fn cache_bounds(&mut self) {
        let groups = GROUPS.read_recursive();
        let group = groups.get(*self).unwrap();
        let bb = Self::make_bounds(group);
        drop(groups);
        GROUPS.write().get_mut(*self).unwrap().bounds = Some(bb);
    }

    fn partition(&mut self) -> (Group, Group) {
        let (left_bb, right_bb) = self.bounds().split();

        let mut groups = GROUPS.write();
        let group = groups.get_mut(*self).unwrap();

        // Split shapes

        let mut remaining_shapes = vec![];
        mem::swap(&mut remaining_shapes, &mut group.shapes);

        let (left_shapes, remaining_shapes) = remaining_shapes
            .into_iter()
            .partition::<Vec<_>, _>(|shape| left_bb.contains_box(shape.parent_space_bounds()));
        let (right_shapes, remaining_shapes) = remaining_shapes
            .into_iter()
            .partition::<Vec<_>, _>(|shape| right_bb.contains_box(shape.parent_space_bounds()));
        group.shapes = remaining_shapes;

        // Split children

        let mut remaining_children = vec![];
        mem::swap(&mut remaining_children, &mut group.children);
        drop(groups);

        let (left_children, remaining_children) = remaining_children
            .into_iter()
            .partition::<Vec<_>, _>(|group| left_bb.contains_box(group.bounds()));
        let (right_children, remaining_children) = remaining_children
            .into_iter()
            .partition::<Vec<_>, _>(|group| right_bb.contains_box(group.bounds()));
        GROUPS.write().get_mut(*self).unwrap().children = remaining_children;

        // Make left group

        let left_group = GROUPS.write().insert(GroupData {
            transform: Matrix::identity(),
            shapes: left_shapes,
            children: left_children.clone(),
            parent: Group::null(),
            bounds: None,
        });

        GROUPS
            .write()
            .get_mut(left_group)
            .unwrap()
            .shapes
            .iter_mut()
            .for_each(|shape| *shape.parent_mut() = left_group);

        let mut groups = GROUPS.write();
        for left_child in left_children {
            groups.get_mut(left_child).unwrap().parent = left_group;
        }
        drop(groups);

        // Make right group

        let right_group = GROUPS.write().insert(GroupData {
            transform: Matrix::identity(),
            shapes: right_shapes,
            children: right_children.clone(),
            parent: Group::null(),
            bounds: None,
        });

        GROUPS
            .write()
            .get_mut(right_group)
            .unwrap()
            .shapes
            .iter_mut()
            .for_each(|shape| *shape.parent_mut() = right_group);

        let mut groups = GROUPS.write();
        for right_child in right_children {
            groups.get_mut(right_child).unwrap().parent = right_group;
        }
        drop(groups);

        // Return partitions

        (left_group, right_group)
    }

    pub fn divide(&mut self, threshold: usize) {
        if threshold <= self.len() {
            let (left, right) = self.partition();
            if !left.is_empty() {
                self.add_child(left);
            }
            if !right.is_empty() {
                self.add_child(right);
            }
        }

        let children = GROUPS.read().get(*self).unwrap().children.clone();
        for mut child in children {
            child.divide(threshold);
        }
    }

    pub fn delete(self) {
        let mut groups = GROUPS.write();
        let group = groups.get_mut(self).unwrap();
        let children = mem::take(&mut group.children);
        drop(groups);
        for child in children {
            child.delete();
        }
        let mut groups = GROUPS.write();
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
            bounds: None,
        };
        let handle = GROUPS.write().insert(group);
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

    #[test]
    fn group_bounding_box() {
        let s = Sphere::default().with_transform(
            Transform::default()
                .scaling(2., 2., 2.)
                .translation(2., 5., -3.),
        );
        let c = Cylinder::default()
            .with_transform(
                Transform::default()
                    .scaling(0.5, 1., 0.5)
                    .translation(-4., -1., 4.),
            )
            .min(-2.)
            .max(2.);
        let mut g = Group::default();
        g.add_shape(s.into());
        g.add_shape(c.into());
        let bb = g.bounds();
        assert!(bb.min == Point::new(-4.5, -3., -5.));
        assert!(bb.max == Point::new(4., 7., 4.5));
    }

    #[test]
    fn group_partitions() {
        let s1 = Shape::sphere().with_transform(translation(-2., 0., 0.));
        let s2 = Shape::sphere().with_transform(translation(2., 0., 0.));
        let s3 = Shape::sphere();
        let mut g = Group::default();
        for s in [s1, s2, s3] {
            g.add_shape(s);
        }
        let (left, right) = g.partition();
        assert!(g.get_shape(0) == s3);
        assert!(left.get_shape(0) == s1);
        assert!(right.get_shape(0) == s2);

        let s1 = Shape::sphere().with_transform(translation(-2., -2., 0.));
        let s2 = Shape::sphere().with_transform(translation(-2., 2., 0.));
        let s3 = Shape::sphere().with_transform(scaling(4., 4., 4.));
        let mut g = Group::default();
        for s in [s1, s2, s3] {
            g.add_shape(s);
        }
        g.divide(1);
        assert!(g.get_shape(0) == s3);
        let sub = g.get_child(0);
        assert!(sub.len() == 2);
        assert!(sub.get_child(0).get_shape(0) == s1);
        assert!(sub.get_child(1).get_shape(0) == s2);

        let s1 = Shape::sphere().with_transform(translation(-2., 0., 0.));
        let s2 = Shape::sphere().with_transform(translation(2., 1., 0.));
        let s3 = Shape::sphere().with_transform(translation(2., -1., 0.));
        let mut sub = Group::default();
        for s in [s1, s2, s3] {
            sub.add_shape(s);
        }
        let s4 = Shape::sphere();
        let mut g = Group::default();
        g.add_shape(s4);
        g.add_child(sub);
        g.divide(3);
        assert!(g.get_child(0) == sub);
        assert!(g.get_shape(0) == s4);
        assert!(sub.len() == 2);
        assert!(sub.get_child(0).get_shape(0) == s1);
        assert!(sub.get_child(1).get_shape(0) == s2);
        assert!(sub.get_child(1).get_shape(1) == s3);
    }
}
