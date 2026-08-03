use std::path::PathBuf;

const PIXI_TOML: &str = include_str!("../pixi.toml");
const PIXI_LOCK: &str = include_str!("../pixi.lock");

/// Write the bundled pixi manifest to a cache directory and return the path to pixi.toml.
/// This allows `pixi run --manifest-path <path> -e <env> <tool>` to work from any
/// directory without the user needing their own pixi workspace.
/// Tries `~/.local/share/galah/` first, falls back to a temp directory.
pub fn galah_manifest_path() -> PathBuf {
    let candidates: Vec<PathBuf> = std::env::var("HOME")
        .map(|h| vec![PathBuf::from(h).join(".local").join("share").join("galah")])
        .unwrap_or_default()
        .into_iter()
        .chain(std::iter::once(std::env::temp_dir().join("galah")))
        .collect();

    for dir in &candidates {
        if std::fs::create_dir_all(dir).is_err() {
            continue;
        }
        let toml_path = dir.join("pixi.toml");
        let lock_path = dir.join("pixi.lock");
        if std::fs::write(&toml_path, PIXI_TOML).is_err() {
            continue;
        }
        if std::fs::write(&lock_path, PIXI_LOCK).is_err() {
            continue;
        }
        return toml_path;
    }

    panic!(
        "Could not write galah pixi manifest to any of: {:?}",
        candidates
    )
}

pub fn is_on_path(name: &str) -> bool {
    if let Some(path) = std::env::var_os("PATH") {
        std::env::split_paths(&path).any(|dir| dir.join(name).is_file())
    } else {
        false
    }
}
