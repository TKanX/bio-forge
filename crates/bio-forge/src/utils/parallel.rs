//! Abstraction layer for parallel iteration.
//!
//! This module provides conditional compilation for parallel processing.
//! When the `parallel` feature is enabled, it exports Rayon's parallelism primitives.
//! When disabled, it provides serial fallbacks that mimic the parallel API,
//! allowing internal code to be written once.

#[cfg(feature = "parallel")]
pub use rayon::prelude::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator,
    IntoParallelRefMutIterator, ParallelBridge, ParallelIterator, ParallelSliceMut,
};

#[cfg(not(feature = "parallel"))]
pub use self::fallback::*;

#[cfg(not(feature = "parallel"))]
mod fallback {
    pub use std::iter::Iterator as ParallelIterator;
    pub use std::iter::Iterator as IndexedParallelIterator;

    /// Shim trait to allow `into_par_iter()` on types that implement `IntoIterator`.
    pub trait IntoParallelIterator {
        type Item;
        type Iter: Iterator<Item = Self::Item>;
        fn into_par_iter(self) -> Self::Iter;
    }

    impl<I: IntoIterator> IntoParallelIterator for I {
        type Item = I::Item;
        type Iter = I::IntoIter;
        fn into_par_iter(self) -> Self::Iter {
            self.into_iter()
        }
    }

    /// Shim trait to allow `par_iter()` on types that implement `IntoIterator` for `&T`.
    pub trait IntoParallelRefIterator<'data> {
        type Item;
        type Iter: Iterator<Item = Self::Item>;
        fn par_iter(&'data self) -> Self::Iter;
    }

    impl<'data, I: 'data + ?Sized> IntoParallelRefIterator<'data> for I
    where
        &'data I: IntoIterator,
    {
        type Item = <&'data I as IntoIterator>::Item;
        type Iter = <&'data I as IntoIterator>::IntoIter;
        fn par_iter(&'data self) -> Self::Iter {
            self.into_iter()
        }
    }

    /// Shim trait to allow `par_iter_mut()` on types that implement `IntoIterator` for `&mut T`.
    pub trait IntoParallelRefMutIterator<'data> {
        type Item;
        type Iter: Iterator<Item = Self::Item>;
        fn par_iter_mut(&'data mut self) -> Self::Iter;
    }

    impl<'data, I: 'data + ?Sized> IntoParallelRefMutIterator<'data> for I
    where
        &'data mut I: IntoIterator,
    {
        type Item = <&'data mut I as IntoIterator>::Item;
        type Iter = <&'data mut I as IntoIterator>::IntoIter;
        fn par_iter_mut(&'data mut self) -> Self::Iter {
            self.into_iter()
        }
    }

    /// Shim trait to allow `par_bridge()` on Iterators.
    pub trait ParallelBridge: Iterator {
        fn par_bridge(self) -> Self;
    }

    impl<T: Iterator> ParallelBridge for T {
        fn par_bridge(self) -> Self {
            self
        }
    }

    /// Shim trait to allow parallel sorting on slices.
    pub trait ParallelSliceMut<T> {
        /// Sorts the slice using the default ordering, potentially in parallel.
        fn par_sort_unstable(&mut self)
        where
            T: Ord;

        /// Sorts the slice using a key extraction function, potentially in parallel.
        fn par_sort_unstable_by_key<K, F>(&mut self, f: F)
        where
            K: Ord,
            F: Fn(&T) -> K;
    }

    impl<T> ParallelSliceMut<T> for [T] {
        fn par_sort_unstable(&mut self)
        where
            T: Ord,
        {
            self.sort_unstable();
        }

        fn par_sort_unstable_by_key<K, F>(&mut self, f: F)
        where
            K: Ord,
            F: Fn(&T) -> K,
        {
            self.sort_unstable_by_key(f)
        }
    }

    /// Extension trait to add Rayon-like methods to standard Iterators.
    pub trait ParallelIteratorExt: Iterator {
        fn flat_map_iter<U, F>(self, f: F) -> std::iter::FlatMap<Self, U, F>
        where
            Self: Sized,
            U: IntoIterator,
            F: FnMut(Self::Item) -> U,
        {
            self.flat_map(f)
        }

        fn try_reduce<T, E, ID, OP>(mut self, identity: ID, mut op: OP) -> Result<T, E>
        where
            Self: Sized + Iterator<Item = Result<T, E>>,
            ID: Fn() -> T,
            OP: FnMut(T, T) -> Result<T, E>,
        {
            let mut acc = identity();
            for item in self {
                match item {
                    Ok(v) => acc = op(acc, v)?,
                    Err(e) => return Err(e),
                }
            }
            Ok(acc)
        }
    }

    impl<I: Iterator> ParallelIteratorExt for I {}
}
