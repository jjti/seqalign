use std::fmt::Debug;

use ordered_float::OrderedFloat;

/// Step is a single step through a pairwise alignment.
#[derive(Clone, Eq)]
pub struct Step {
    pub val: OrderedFloat<f32>,
    pub i: usize,
    pub j: usize,
    pub next: Option<(usize, usize)>,
}

impl Step {
    pub fn from(i: usize, j: usize, val: f32) -> Self {
        Step {
            val: OrderedFloat(val),
            i,
            j,
            next: None,
        }
    }
}

impl Default for Step {
    fn default() -> Self {
        Step {
            val: OrderedFloat(f32::MIN),
            i: usize::MIN,
            j: usize::MIN,
            next: None,
        }
    }
}

impl Debug for Step {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {}): {}", self.i, self.j, self.val)
    }
}

impl PartialEq for Step {
    fn eq(&self, other: &Self) -> bool {
        self.i == other.i && self.j == other.j && self.val == other.val
    }
}

impl Ord for Step {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.val.cmp(&other.val)
    }
}

impl PartialOrd for Step {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.val.partial_cmp(&other.val)
    }
}
