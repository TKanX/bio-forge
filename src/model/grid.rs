//! Spatial indexing primitives for accelerating geometric queries.
//!
//! This module provides a [`Grid`] structure that partitions 3D space into uniform cells,
//! enabling **O(1)** average-case lookups for neighbor searches, collision detection, and
//! range queries.

use super::types::Point;
use nalgebra::Vector3;

/// Sentinel value indicating the end of a linked list.
const SENTINEL: u32 = u32::MAX;

/// A uniform spatial grid that bins items into cubic cells.
///
/// The grid is defined by a bounding box and a cell size. Items are mapped to cells
/// based on their coordinates. This structure is optimized for "fixed-radius" queries,
/// where the search radius is comparable to the cell size.
///
/// # Performance
///
/// - Construction: **O(N)** where N is the number of items.
/// - Neighbor queries: **O(1)** average-case per query, assuming uniform distribution.
#[derive(Debug, Clone)]
pub struct Grid<T> {
    /// Side length of each cubic cell.
    cell_size: f64,
    /// Minimum coordinate of the grid's bounding box.
    origin: Point,
    /// Number of cells along each dimension (x, y, z).
    dims: Vector3<usize>,
    /// Index of the first item in each cell. Size = num_cells.
    head: Vec<u32>,
    /// Index of the next item in the linked list. Size = num_items.
    next: Vec<u32>,
    /// Stored items with their positions. Size = num_items.
    items: Vec<(Point, T)>,
}

impl<T> Grid<T> {
    /// Creates a new grid enclosing the provided points.
    ///
    /// The grid dimensions are automatically calculated to encompass all points with
    /// a small padding.
    ///
    /// # Arguments
    ///
    /// * `items` - Iterator yielding `(position, item)` pairs.
    /// * `cell_size` - The side length of each spatial bin.
    ///
    /// # Panics
    ///
    /// Panics if `cell_size` is non-positive.
    pub fn new(items: impl IntoIterator<Item = (Point, T)>, cell_size: f64) -> Self {
        assert!(cell_size > 0.0, "Cell size must be positive");

        let input_items: Vec<_> = items.into_iter().collect();
        let num_items = input_items.len();

        if num_items == 0 {
            return Self {
                cell_size,
                origin: Point::origin(),
                dims: Vector3::zeros(),
                head: Vec::new(),
                next: Vec::new(),
                items: Vec::new(),
            };
        }

        let mut min = Point::new(f64::MAX, f64::MAX, f64::MAX);
        let mut max = Point::new(f64::MIN, f64::MIN, f64::MIN);

        for (pos, _) in &input_items {
            min = min.inf(pos);
            max = max.sup(pos);
        }

        let epsilon = 1e-6;
        max += Vector3::new(epsilon, epsilon, epsilon);

        let extent = max - min;
        let dims = Vector3::new(
            (extent.x / cell_size).ceil() as usize,
            (extent.y / cell_size).ceil() as usize,
            (extent.z / cell_size).ceil() as usize,
        );

        let total_cells = dims.x * dims.y * dims.z;

        let mut head = vec![SENTINEL; total_cells];
        let mut next = vec![SENTINEL; num_items];
        let mut stored_items = Vec::with_capacity(num_items);

        for (i, (pos, item)) in input_items.into_iter().enumerate() {
            stored_items.push((pos, item));

            if let Some(cell_idx) = Self::get_cell_index_static(&pos, dims, min, cell_size) {
                next[i] = head[cell_idx];
                head[cell_idx] = i as u32;
            }
        }

        Self {
            cell_size,
            origin: min,
            dims,
            head,
            next,
            items: stored_items,
        }
    }

    /// Static helper to compute cell index without `self`.
    fn get_cell_index_static(
        pos: &Point,
        dims: Vector3<usize>,
        origin: Point,
        cell_size: f64,
    ) -> Option<usize> {
        if pos.x < origin.x || pos.y < origin.y || pos.z < origin.z {
            return None;
        }

        let offset = pos - origin;
        let x = (offset.x / cell_size).floor() as usize;
        let y = (offset.y / cell_size).floor() as usize;
        let z = (offset.z / cell_size).floor() as usize;

        if x >= dims.x || y >= dims.y || z >= dims.z {
            return None;
        }

        Some(x + y * dims.x + z * dims.x * dims.y)
    }

    /// Iterates over all items in cells overlapping with the query sphere.
    ///
    /// The returned iterator yields items from candidate cells. To filter strictly by
    /// Euclidean distance, use the `.exact()` method on the returned iterator.
    ///
    /// # Arguments
    ///
    /// * `center` - Center of the search sphere.
    /// * `radius` - Radius of the search sphere.
    pub fn neighbors<'a>(&'a self, center: &Point, radius: f64) -> GridNeighborhood<'a, T> {
        if self.items.is_empty() {
            return GridNeighborhood {
                grid: self,
                min_x: 0,
                max_x: 0,
                min_y: 0,
                max_y: 0,
                max_z: 0,
                curr_x: 0,
                curr_y: 0,
                curr_z: 1,
                curr_item_idx: SENTINEL,
                center: *center,
                radius_sq: radius * radius,
            };
        }

        let min_idx = self.get_grid_coords(&(center - Vector3::new(radius, radius, radius)));
        let max_idx = self.get_grid_coords(&(center + Vector3::new(radius, radius, radius)));

        let (min_x, min_y, min_z) = min_idx;
        let (max_x, max_y, max_z) = max_idx;

        GridNeighborhood {
            grid: self,
            min_x,
            max_x,
            min_y,
            max_y,
            max_z,
            curr_x: min_x,
            curr_y: min_y,
            curr_z: min_z,
            curr_item_idx: SENTINEL,
            center: *center,
            radius_sq: radius * radius,
        }
    }

    /// Helper to get clamped grid coordinates (x, y, z).
    fn get_grid_coords(&self, pos: &Point) -> (usize, usize, usize) {
        let offset = pos - self.origin;
        let x = (offset.x / self.cell_size).floor() as isize;
        let y = (offset.y / self.cell_size).floor() as isize;
        let z = (offset.z / self.cell_size).floor() as isize;

        (
            x.clamp(0, (self.dims.x as isize) - 1) as usize,
            y.clamp(0, (self.dims.y as isize) - 1) as usize,
            z.clamp(0, (self.dims.z as isize) - 1) as usize,
        )
    }

    /// Checks if any item in the grid is within `radius` of `point`.
    ///
    /// This is optimized to return early.
    ///
    /// # Arguments
    ///
    /// * `point` - The query point.
    /// * `radius` - The cutoff distance.
    /// * `predicate` - A closure to filter items (e.g., check exact distance).
    pub fn has_neighbor<F>(&self, point: &Point, radius: f64, mut predicate: F) -> bool
    where
        F: FnMut(&T) -> bool,
    {
        for item in self.neighbors(point, radius) {
            if predicate(item) {
                return true;
            }
        }
        false
    }
}

/// Iterator for traversing grid cells and their linked lists.
///
/// This iterator yields all items in the cells that overlap with the query sphere.
/// Use [`GridNeighborhood::exact`] to filter items strictly within the radius.
pub struct GridNeighborhood<'a, T> {
    grid: &'a Grid<T>,
    min_x: usize,
    max_x: usize,
    min_y: usize,
    max_y: usize,
    max_z: usize,
    curr_x: usize,
    curr_y: usize,
    curr_z: usize,
    curr_item_idx: u32,
    center: Point,
    radius_sq: f64,
}

impl<'a, T> GridNeighborhood<'a, T> {
    /// Returns an iterator that yields only items strictly within the search radius.
    ///
    /// This method uses the internally stored positions to perform the distance check,
    /// avoiding the need for external lookups.
    pub fn exact(self) -> impl Iterator<Item = &'a T> + 'a {
        ExactGridNeighborhood { inner: self }
    }
}

/// Iterator that yields items strictly within the Euclidean radius.
///
/// This iterator filters candidates from [`GridNeighborhood`] using the stored positions
/// and the query center/radius, returning only items whose distance to the center is
/// less than or equal to the specified radius.
pub struct ExactGridNeighborhood<'a, T> {
    inner: GridNeighborhood<'a, T>,
}

impl<'a, T> Iterator for ExactGridNeighborhood<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.inner.curr_item_idx != SENTINEL {
                let (pos, item) = &self.inner.grid.items[self.inner.curr_item_idx as usize];
                self.inner.curr_item_idx = self.inner.grid.next[self.inner.curr_item_idx as usize];

                if nalgebra::distance_squared(pos, &self.inner.center) <= self.inner.radius_sq {
                    return Some(item);
                }
                continue;
            }

            if self.inner.curr_x > self.inner.max_x {
                self.inner.curr_x = self.inner.min_x;
                self.inner.curr_y += 1;
            }
            if self.inner.curr_y > self.inner.max_y {
                self.inner.curr_y = self.inner.min_y;
                self.inner.curr_z += 1;
            }
            if self.inner.curr_z > self.inner.max_z {
                return None;
            }

            let cell_idx = self.inner.curr_x
                + self.inner.curr_y * self.inner.grid.dims.x
                + self.inner.curr_z * self.inner.grid.dims.x * self.inner.grid.dims.y;

            self.inner.curr_x += 1;

            if cell_idx < self.inner.grid.head.len() {
                self.inner.curr_item_idx = self.inner.grid.head[cell_idx];
            }
        }
    }
}

impl<'a, T> Iterator for GridNeighborhood<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.curr_item_idx != SENTINEL {
                let (_, item) = &self.grid.items[self.curr_item_idx as usize];
                self.curr_item_idx = self.grid.next[self.curr_item_idx as usize];
                return Some(item);
            }

            if self.curr_x > self.max_x {
                self.curr_x = self.min_x;
                self.curr_y += 1;
            }
            if self.curr_y > self.max_y {
                self.curr_y = self.min_y;
                self.curr_z += 1;
            }
            if self.curr_z > self.max_z {
                return None;
            }

            let cell_idx = self.curr_x
                + self.curr_y * self.grid.dims.x
                + self.curr_z * self.grid.dims.x * self.grid.dims.y;

            self.curr_x += 1;

            if cell_idx < self.grid.head.len() {
                self.curr_item_idx = self.grid.head[cell_idx];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::types::Point;

    #[test]
    fn grid_bins_points_correctly() {
        let points = vec![
            (Point::new(0.5, 0.5, 0.5), 1),
            (Point::new(1.5, 0.5, 0.5), 2),
            (Point::new(0.5, 1.5, 0.5), 3),
        ];

        let grid = Grid::new(points, 1.0);

        assert_eq!(grid.dims, Vector3::new(2, 2, 1));

        let center = Point::new(0.5, 0.5, 0.5);
        let neighbors: Vec<_> = grid.neighbors(&center, 0.1).collect();
        assert!(neighbors.contains(&&1));
        assert!(!neighbors.contains(&&2));
    }

    #[test]
    fn grid_neighbors_returns_nearby_items() {
        let points = vec![
            (Point::new(0.0, 0.0, 0.0), "A"),
            (Point::new(10.0, 0.0, 0.0), "B"),
        ];
        let grid = Grid::new(points, 2.0);

        let center = Point::new(0.1, 0.1, 0.1);
        let neighbors: Vec<_> = grid.neighbors(&center, 1.0).collect();
        assert_eq!(neighbors.len(), 1);
        assert_eq!(*neighbors[0], "A");
    }

    #[test]
    fn grid_handles_empty_input() {
        let points: Vec<(Point, i32)> = vec![];
        let grid = Grid::new(points, 1.0);
        assert_eq!(grid.items.len(), 0);
        assert_eq!(grid.neighbors(&Point::origin(), 1.0).count(), 0);
    }

    #[test]
    fn grid_handles_dense_packing() {
        let mut points = Vec::new();
        for i in 0..100 {
            points.push((Point::new(0.1, 0.1, 0.1), i));
        }
        let grid = Grid::new(points, 1.0);

        let center = Point::new(0.1, 0.1, 0.1);
        let count = grid.neighbors(&center, 0.5).count();
        assert_eq!(count, 100);
    }

    #[test]
    fn grid_handles_boundary_conditions() {
        let points = vec![
            (Point::new(0.0, 0.0, 0.0), 1),
            (Point::new(10.0, 10.0, 10.0), 2),
        ];
        let grid = Grid::new(points, 1.0);

        let center = Point::new(0.0, 0.0, 0.0);
        assert!(grid.has_neighbor(&center, 0.1, |&i| i == 1));

        let center = Point::new(10.0, 10.0, 10.0);
        assert!(grid.has_neighbor(&center, 0.1, |&i| i == 2));
    }

    #[test]
    fn grid_exact_filtering_works() {
        let points = vec![
            (Point::new(0.0, 0.0, 0.0), "Center"),
            (Point::new(0.9, 0.0, 0.0), "Inside"),
            (Point::new(1.1, 0.0, 0.0), "Outside"),
        ];
        let grid = Grid::new(points, 2.0);

        let center = Point::new(0.0, 0.0, 0.0);
        let radius = 1.0;

        let coarse_count = grid.neighbors(&center, radius).count();
        assert_eq!(coarse_count, 3);

        let exact_neighbors: Vec<_> = grid.neighbors(&center, radius).exact().collect();
        assert_eq!(exact_neighbors.len(), 2);
        assert!(exact_neighbors.contains(&&"Center"));
        assert!(exact_neighbors.contains(&&"Inside"));
        assert!(!exact_neighbors.contains(&&"Outside"));
    }

    #[test]
    fn grid_exact_filtering_handles_empty_grid() {
        let points: Vec<(Point, i32)> = vec![];
        let grid = Grid::new(points, 1.0);

        let count = grid.neighbors(&Point::origin(), 1.0).exact().count();
        assert_eq!(count, 0);
    }

    #[test]
    fn grid_exact_filtering_handles_boundary_points() {
        let points = vec![(Point::new(1.0, 0.0, 0.0), "OnBoundary")];
        let grid = Grid::new(points, 2.0);

        let center = Point::new(0.0, 0.0, 0.0);
        let count = grid.neighbors(&center, 1.0).exact().count();
        assert_eq!(count, 1);

        let count = grid.neighbors(&center, 0.99).exact().count();
        assert_eq!(count, 0);
    }
}
