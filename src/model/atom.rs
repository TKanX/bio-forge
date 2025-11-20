use super::types::{Element, Point};
use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub name: String,
    pub element: Element,
    pub pos: Point,
}

impl Atom {
    pub fn new(name: &str, element: Element, pos: Point) -> Self {
        Self {
            name: name.to_string(),
            element,
            pos,
        }
    }

    pub fn distance_squared(&self, other: &Atom) -> f64 {
        nalgebra::distance_squared(&self.pos, &other.pos)
    }

    pub fn distance(&self, other: &Atom) -> f64 {
        nalgebra::distance(&self.pos, &other.pos)
    }

    pub fn translate_by(&mut self, vector: &nalgebra::Vector3<f64>) {
        self.pos += vector;
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Atom {{ name: \"{}\", element: {}, pos: [{:.3}, {:.3}, {:.3}] }}",
            self.name, self.element, self.pos.x, self.pos.y, self.pos.z
        )
    }
}
