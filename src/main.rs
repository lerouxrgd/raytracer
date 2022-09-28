use std::error::Error;
use std::path::PathBuf;
use std::{fs, io};

use clap::Parser;
use raytracer::scene::Scene;

/// The Ray Tracer Challenge CLI
#[derive(Debug, Parser)]
struct Args {
    /// A yaml description of the scene to render
    #[clap(long, value_name = "FILE", display_order = 1)]
    scene: PathBuf,
    /// Optional obj models to add to the scene
    #[clap(long, value_name = "FILE", display_order = 2)]
    obj: Option<Vec<PathBuf>>,
    /// Optional ppm textures to use as material
    #[clap(long, value_name = "FILE", display_order = 3)]
    ppm: Option<Vec<PathBuf>>,
    /// Optional output ppm file, defaults to stdout
    #[clap(long, display_order = 3)]
    out: Option<PathBuf>,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let scene = fs::read_to_string(args.scene)?;
    let scene: Scene = serde_yaml::from_str(&scene)?;

    let obj = args.obj.unwrap_or_default();
    let ppm = args.ppm.unwrap_or_default();
    if let Some(out) = args.out {
        let out = fs::File::create(out)?;
        scene.render(out, &obj, &ppm)?;
    } else {
        scene.render(io::stdout(), &obj, &ppm)?;
    }

    Ok(())
}
