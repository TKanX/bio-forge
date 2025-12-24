//! Abstraction layer for parallel iteration.
//!
//! This module provides conditional compilation for parallel processing.
//! When the `parallel` feature is enabled, it exports Rayon's parallelism primitives.
//! When disabled, it provides serial fallbacks that mimic the parallel API,
//! allowing internal code to be written once.

#[cfg(feature = "parallel")]
pub use rayon::prelude::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator,
    IntoParallelRefMutIterator, ParallelBridge, ParallelIterator,
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
}
