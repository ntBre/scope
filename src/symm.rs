#![allow(dead_code)]

use std::{collections::HashMap, str::FromStr, string::ParseError};

use nalgebra as na;

// TODO expand beyond cartesian axes. an alternative formulation of this is to
// align the geometry to a cartesian axis if it doesn't start like that. I think
// rotations pretty much assume you are along the cartesian axes

// restrict these to the cartesian axes for now
#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub enum Axis {
    X,
    Y,
    Z,
}

// restrict these to combinations of cartesian axes for now. a more general
// plane is described by (a, b, c) in the equation ax + by + cz = 0
#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub struct Plane(Axis, Axis);

#[derive(Debug, Clone)]
pub struct Atom {
    pub atomic_number: usize,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl PartialEq for Atom {
    fn eq(&self, other: &Self) -> bool {
        let eps = 1e-8;
        let close = |a: f64, b: f64| (a - b).abs() < eps;
        self.atomic_number == other.atomic_number
            && close(self.x, other.x)
            && close(self.y, other.y)
            && close(self.z, other.z)
    }
}

impl Atom {
    pub fn new(atomic_number: usize, x: f64, y: f64, z: f64) -> Self {
        Self {
            atomic_number,
            x,
            y,
            z,
        }
    }
}

#[derive(Debug, PartialEq)]
pub enum PointGroup {
    C1,
    C2 { axis: Axis },
    Cs { plane: Plane },
    C2v { axis: Axis, planes: Vec<Plane> },
}

pub enum Irrep {
    // C1
    A,
    // C2
    B,
    // Cs - p = prime
    Ap,
    App,
    // C2v
    A1,
    B1,
    B2,
    A2,
}

#[derive(Debug, Clone)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
}

impl PartialEq for Molecule {
    /// compare molecules for equality, irrespective of order. try to find an
    /// atom in other that equals the current atom in self. If found, remove it,
    /// so it can't be double-counted.
    fn eq(&self, other: &Self) -> bool {
        let mut theirs = other.atoms.clone();
        if self.atoms.len() != theirs.len() {
            return false;
        }
        for atom in &self.atoms {
            let mut pops = Vec::new();
            let mut found = false;
            for (i, btom) in theirs.iter().enumerate() {
                if *atom == *btom {
                    pops.push(i);
                    found = true;
                    break;
                }
            }
            if !found {
                return false;
            }
            // remove high indices first
            pops.sort();
            pops.reverse();
            for p in pops {
                theirs.remove(p);
            }
        }
        true
    }
}

impl FromStr for Molecule {
    type Err = ParseError;

    /// parse lines like
    ///      O           0.000000000    0.000000000   -0.124238453
    ///      H           0.000000000    1.431390207    0.986041184
    ///      H           0.000000000   -1.431390207    0.986041184
    /// into a molecule
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut ret = Self::default();
        let atomic_symbols = HashMap::from([("H", 1), ("C", 6), ("O", 8)]);
        for line in s.lines() {
            let fields = line.split_whitespace().collect::<Vec<_>>();
            if fields.len() == 4 {
                let sym = if let Some(&s) = atomic_symbols.get(fields[0]) {
                    s
                } else {
                    panic!(
                        "atomic symbol '{}' not found, tell Brent!",
                        fields[0]
                    );
                };
                ret.atoms.push(Atom::new(
                    sym,
                    fields[1].parse().unwrap(),
                    fields[2].parse().unwrap(),
                    fields[3].parse().unwrap(),
                ));
            }
        }
        Ok(ret)
    }
}

impl Molecule {
    pub fn default() -> Self {
        Self { atoms: Vec::new() }
    }

    fn to_vecs(&self) -> Vec<na::Vector3<f64>> {
        let mut ret = Vec::with_capacity(self.atoms.len());
        for atom in &self.atoms {
            ret.push(na::Vector3::new(atom.x, atom.y, atom.z));
        }
        ret
    }

    pub fn point_group(&self) -> PointGroup {
        use Axis::*;
        use PointGroup::*;
        let mut axes = Vec::new();
        let mut planes = Vec::new();
        for ax in vec![X, Y, Z] {
            if self.rotate(180.0, &ax) == *self {
                axes.push(ax);
            }
        }
        for plane in vec![Plane(X, Y), Plane(X, Z), Plane(Y, Z)] {
            if self.reflect(&plane) == *self {
                planes.push(plane);
            }
        }
        match (axes.len(), planes.len()) {
            (0, 1) => Cs { plane: planes[0] },
            (1, 0) => C2 { axis: axes[0] },
            (1, 2) => C2v {
                planes,
                axis: axes[0],
            },
            _ => C1,
        }
    }

    /// apply the transformation matrix `mat` to the atoms in `self` and return
    /// the new Molecule
    fn transform(&self, mat: na::Matrix3<f64>) -> Self {
        let mut ret = Self::default();
        for (i, atom) in self.to_vecs().iter().enumerate() {
            let v = mat * atom;
            ret.atoms.push(Atom::new(
                self.atoms[i].atomic_number,
                v[0],
                v[1],
                v[2],
            ));
        }
        ret
    }

    pub fn rotate(&self, deg: f64, axis: &Axis) -> Self {
        use Axis::*;
        let deg = deg.to_radians();
        let ct = deg.cos();
        let st = deg.sin();
        // from
        // https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
        #[rustfmt::skip]
	let rot_mat = match axis {
            X => {
		na::Matrix3::new(
		    1., 0., 0.,
		    0., ct, -st,
		    0., st, ct,
		)
            }
            Y => {
		na::Matrix3::new(
		    ct, 0., st,
		    0., 1., 0.,
		    -st, 0., ct,
		)
            }
            Z => {
		na::Matrix3::new(
		    ct, -st, 0.,
		    st, ct, 0.,
		    0., 0., 1.,
		)
            }
        };
        self.transform(rot_mat)
    }

    #[rustfmt::skip]
    /// return the special case of the Householder reflection in 3 dimensions
    /// described here:
    /// <https://en.wikipedia.org/wiki/Transformation_matrix#Reflection_2>
    fn householder(a: f64, b: f64, c: f64) -> na::Matrix3<f64> {
        na::Matrix3::new(
            1. - 2. * a * a, -2. * a * b, -2. * a * c,
            -2. * a * b, 1. - 2. * b * b, -2. * b * c,
            -2. * a * c, -2. * b * c, 1. - 2. * c * c,
        )
    }

    pub fn reflect(&self, plane: &Plane) -> Self {
        use Axis::*;
        let ref_mat = match plane {
            Plane(X, Y) | Plane(Y, X) => Self::householder(0.0, 0.0, 1.0),
            Plane(X, Z) | Plane(Z, X) => Self::householder(0.0, 1.0, 0.0),
            Plane(Y, Z) | Plane(Z, Y) => Self::householder(1.0, 0.0, 0.0),
            _ => panic!("unrecognized plane {:?}", plane),
        };
        self.transform(ref_mat)
    }

    pub fn irrep(&self, _pg: PointGroup) -> Irrep {
        // use Irrep::*;
        todo!();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use Axis::*;

    #[test]
    fn test_rotate() {
        let tests = vec![
            // X around all axes
            (
                vec![Atom::new(1, 1.0, 0.0, 0.0)],
                vec![Atom::new(1, 1.0, 0.0, 0.0)],
                180.0,
                X,
            ),
            (
                vec![Atom::new(1, 1.0, 0.0, 0.0)],
                vec![Atom::new(1, -1.0, 0.0, 0.0)],
                180.0,
                Y,
            ),
            (
                vec![Atom::new(1, 1.0, 0.0, 0.0)],
                vec![Atom::new(1, -1.0, 0.0, 0.0)],
                180.0,
                Z,
            ),
            // Y around all axes
            (
                vec![Atom::new(1, 0.0, 1.0, 0.0)],
                vec![Atom::new(1, 0.0, -1.0, 0.0)],
                180.0,
                X,
            ),
            (
                vec![Atom::new(1, 0.0, 1.0, 0.0)],
                vec![Atom::new(1, 0.0, 1.0, 0.0)],
                180.0,
                Y,
            ),
            (
                vec![Atom::new(1, 0.0, 1.0, 0.0)],
                vec![Atom::new(1, 0.0, -1.0, 0.0)],
                180.0,
                Z,
            ),
            // Z around all axes
            (
                vec![Atom::new(1, 0.0, 0.0, 1.0)],
                vec![Atom::new(1, 0.0, 0.0, -1.0)],
                180.0,
                X,
            ),
            (
                vec![Atom::new(1, 0.0, 0.0, 1.0)],
                vec![Atom::new(1, 0.0, 0.0, -1.0)],
                180.0,
                Y,
            ),
            (
                vec![Atom::new(1, 0.0, 0.0, 1.0)],
                vec![Atom::new(1, 0.0, 0.0, 1.0)],
                180.0,
                Z,
            ),
        ];
        for test in tests {
            let h = Molecule { atoms: test.0 };
            let want = Molecule { atoms: test.1 };
            let got = h.rotate(test.2, &test.3);
            assert_eq!(got, want);
        }
    }

    #[test]
    fn test_reflect() {
        let tests = vec![
            // X through all the planes
            (
                vec![Atom::new(1, 1.0, 0.0, 0.0)],
                vec![Atom::new(1, -1.0, 0.0, 0.0)],
                Plane(Y, Z),
            ),
            (
                vec![Atom::new(1, 1.0, 0.0, 0.0)],
                vec![Atom::new(1, 1.0, 0.0, 0.0)],
                Plane(X, Z),
            ),
            (
                vec![Atom::new(1, 1.0, 0.0, 0.0)],
                vec![Atom::new(1, 1.0, 0.0, 0.0)],
                Plane(X, Y),
            ),
            // Y through all the planes
            (
                vec![Atom::new(1, 0.0, 1.0, 0.0)],
                vec![Atom::new(1, 0.0, 1.0, 0.0)],
                Plane(Y, Z),
            ),
            (
                vec![Atom::new(1, 0.0, 1.0, 0.0)],
                vec![Atom::new(1, 0.0, -1.0, 0.0)],
                Plane(X, Z),
            ),
            (
                vec![Atom::new(1, 0.0, 1.0, 0.0)],
                vec![Atom::new(1, 0.0, 1.0, 0.0)],
                Plane(X, Y),
            ),
            // Z through all the planes
            (
                vec![Atom::new(1, 0.0, 0.0, 1.0)],
                vec![Atom::new(1, 0.0, 0.0, 1.0)],
                Plane(Y, Z),
            ),
            (
                vec![Atom::new(1, 0.0, 0.0, 1.0)],
                vec![Atom::new(1, 0.0, 0.0, 1.0)],
                Plane(X, Z),
            ),
            (
                vec![Atom::new(1, 0.0, 0.0, 1.0)],
                vec![Atom::new(1, 0.0, 0.0, -1.0)],
                Plane(X, Y),
            ),
        ];
        for test in tests {
            let h = Molecule { atoms: test.0 };
            let want = Molecule { atoms: test.1 };
            let got = h.reflect(&test.2);
            assert_eq!(got, want);
        }
    }

    #[test]
    fn test_point_group() {
        use PointGroup::*;
        let mol = Molecule::from_str(
            "
  O           0.000000000    0.000000000   -0.124238453
  H           0.000000000    1.431390207    0.986041184
  H           0.000000000   -1.431390207    0.986041184
",
        )
        .unwrap();
        assert_eq!(
            mol.point_group(),
            C2v {
                axis: Z,
                planes: vec![Plane(X, Z), Plane(Y, Z)]
            }
        );

        let mol = Molecule::from_str(
            "
    C        0.000000   -0.888844    0.000000
    C       -0.662697    0.368254    0.000000
    C        0.662697    0.368254    0.000000
    H       -1.595193    0.906925    0.000000
    H        1.595193    0.906925    0.000000
",
        )
        .unwrap();
        assert_eq!(
            mol.point_group(),
            C2v {
                axis: Y,
                planes: vec![Plane(X, Y), Plane(Y, Z)]
            }
        );
    }
}
