use anyhow::{Result, anyhow, bail};
use clap::Args;

use bio_forge::Structure;
use bio_forge::ops::Transform;

use crate::commands::run_with_spinner;

/// Applies geometric transforms to the structure.
#[derive(Debug, Args)]
pub struct TransformArgs {
    /// Move the geometric center to the origin.
    #[arg(long)]
    pub center: bool,
    /// Move the center of mass to the origin.
    #[arg(long = "center-mass")]
    pub center_mass: bool,
    /// Translate by the provided vector (format: "x,y,z").
    #[arg(long, value_name = "X,Y,Z")]
    pub translate: Option<String>,
    /// Rotate around the X axis (degrees).
    #[arg(long = "rotate-x", value_name = "DEG")]
    pub rotate_x: Option<f64>,
    /// Rotate around the Y axis (degrees).
    #[arg(long = "rotate-y", value_name = "DEG")]
    pub rotate_y: Option<f64>,
    /// Rotate around the Z axis (degrees).
    #[arg(long = "rotate-z", value_name = "DEG")]
    pub rotate_z: Option<f64>,
}

/// Executes the configured transformations in the prescribed order.
pub fn run(structure: &mut Structure, args: &TransformArgs) -> Result<()> {
    run_with_spinner("Applying transforms", || {
        if args.center {
            Transform::center_geometry(structure, None);
        }

        if args.center_mass {
            Transform::center_mass(structure, None);
        }

        if let Some(deg) = args.rotate_x {
            Transform::rotate_x(structure, deg.to_radians());
        }
        if let Some(deg) = args.rotate_y {
            Transform::rotate_y(structure, deg.to_radians());
        }
        if let Some(deg) = args.rotate_z {
            Transform::rotate_z(structure, deg.to_radians());
        }

        if let Some(vector) = &args.translate {
            let (x, y, z) = parse_vector(vector)?;
            Transform::translate(structure, x, y, z);
        }

        Ok(())
    })
}

fn parse_vector(value: &str) -> Result<(f64, f64, f64)> {
    let parts: Vec<_> = value
        .split(',')
        .map(|p| p.trim())
        .filter(|p| !p.is_empty())
        .collect();
    if parts.len() != 3 {
        bail!("Invalid translation '{}'. Expected format 'x,y,z'.", value);
    }

    let x = parts[0]
        .parse::<f64>()
        .map_err(|_| anyhow!("Invalid X component in translation: '{}'", parts[0]))?;
    let y = parts[1]
        .parse::<f64>()
        .map_err(|_| anyhow!("Invalid Y component in translation: '{}'", parts[1]))?;
    let z = parts[2]
        .parse::<f64>()
        .map_err(|_| anyhow!("Invalid Z component in translation: '{}'", parts[2]))?;
    Ok((x, y, z))
}
