//! Read-only access to residue templates.

use crate::db;

pub use crate::db::TemplateView;

/// Looks up a template by its canonical name (for example `"ALA"` or `"HOH"`).
pub fn get(name: &str) -> Option<TemplateView<'_>> {
    db::get_template(name)
}
