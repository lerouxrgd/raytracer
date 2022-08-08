use std::io::Read;

use crate::groups::Group;
use crate::shapes::{SmoothTriangle, Triangle};
use crate::tuples::{Point, Vector};
use crate::wavefront;

pub fn parse_obj<R: Read>(reader: R) -> Result<Group, wavefront::Error> {
    let obj = wavefront::Obj::from_reader(reader)?;

    let mut groups = Group::default();
    for (_, g) in obj.groups() {
        let mut group = Group::default();
        for [a, b, c] in g.triangles() {
            let [a1, a2, a3] = a.position();
            let p1 = Point::new(a1, a2, a3);

            let [b1, b2, b3] = b.position();
            let p2 = Point::new(b1, b2, b3);

            let [c1, c2, c3] = c.position();
            let p3 = Point::new(c1, c2, c3);

            let tri = match (a.normal(), b.normal(), c.normal()) {
                (Some([a1, a2, a3]), Some([b1, b2, b3]), Some([c1, c2, c3])) => {
                    let n1 = Vector::new(a1, a2, a3);
                    let n2 = Vector::new(b1, b2, b3);
                    let n3 = Vector::new(c1, c2, c3);
                    SmoothTriangle::new(p1, p2, p3, n1, n2, n3).into()
                }
                _ => Triangle::new(p1, p2, p3).into(),
            };
            group.add_shape(tri);
        }
        groups.add_child(group);
    }

    Ok(groups)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn obj_triangles() {
        let content = "
v -1 1 0
v -1 0 0
v 1 0 0
v 1 1 0

f 1 2 3
f 1 3 4
";
        let reader = Cursor::new(content.as_bytes());
        let g = parse_obj(reader).unwrap();

        assert!(
            g.get_child(0).get_shape(0)
                == Triangle::new(
                    Point::new(-1., 1., 0.),
                    Point::new(-1., 0., 0.),
                    Point::new(1., 0., 0.),
                )
                .into()
        );
        assert!(
            g.get_child(0).get_shape(1)
                == Triangle::new(
                    Point::new(-1., 1., 0.),
                    Point::new(1., 0., 0.),
                    Point::new(1., 1., 0.),
                )
                .into()
        );
    }

    #[test]
    fn obj_polygon() {
        let content = "
v -1 1 0
v -1 0 0
v 1 0 0
v 1 1 0
v 0 2 0

f 1 2 3 4 5
";
        let reader = Cursor::new(content.as_bytes());
        let g = parse_obj(reader).unwrap();

        assert!(
            g.get_child(0).get_shape(0)
                == Triangle::new(
                    Point::new(-1., 1., 0.),
                    Point::new(-1., 0., 0.),
                    Point::new(1., 0., 0.),
                )
                .into()
        );
        assert!(
            g.get_child(0).get_shape(1)
                == Triangle::new(
                    Point::new(-1., 1., 0.),
                    Point::new(1., 0., 0.),
                    Point::new(1., 1., 0.),
                )
                .into()
        );
        assert!(
            g.get_child(0).get_shape(2)
                == Triangle::new(
                    Point::new(-1., 1., 0.),
                    Point::new(1., 1., 0.),
                    Point::new(0., 2., 0.),
                )
                .into()
        );
    }

    #[test]
    fn obj_groups() {
        let content = "
v -1 1 0
v -1 0 0
v 1 0 0
v 1 1 0

g FirstGroup
f 1 2 3

g SecondGroup
f 1 3 4
";
        let reader = Cursor::new(content.as_bytes());
        let g = parse_obj(reader).unwrap();

        let triangles = vec![g.get_child(0).get_shape(0), g.get_child(1).get_shape(0)];
        assert!(triangles.contains(
            &Triangle::new(
                Point::new(-1., 1., 0.),
                Point::new(-1., 0., 0.),
                Point::new(1., 0., 0.),
            )
            .into()
        ));
        assert!(triangles.contains(
            &Triangle::new(
                Point::new(-1., 1., 0.),
                Point::new(1., 0., 0.),
                Point::new(1., 1., 0.),
            )
            .into()
        ));
    }

    #[test]
    fn obj_smooth_triangles() {
        let content = "
v 0 1 0
v -1 0 0
v 1 0 0

vn -1 0 0
vn 1 0 0
vn 0 1 0

f 1//3 2//1 3//2
f 1//3 2//1 3//2
";
        let reader = Cursor::new(content.as_bytes());
        let g = parse_obj(reader).unwrap();

        assert!(
            g.get_child(0).get_shape(0)
                == SmoothTriangle::new(
                    Point::new(0., 1., 0.),
                    Point::new(-1., 0., 0.),
                    Point::new(1., 0., 0.),
                    Vector::new(0., 1., 0.),
                    Vector::new(-1., 0., 0.),
                    Vector::new(1., 0., 0.),
                )
                .into()
        );
        assert!(g.get_child(0).get_shape(0) == g.get_child(0).get_shape(1));
    }
}
