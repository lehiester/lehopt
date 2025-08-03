//! Cargo build script
//! (see https://doc.rust-lang.org/cargo/reference/build-scripts.html)

use cc;
use std::fs::read_dir;


fn main() {
    // get all the .c files for CSparse
    let dot_c_files: Vec<String> =
        read_dir("deps/CSparse_modified/Source")
        .unwrap()
        .filter(|e| e.as_ref().unwrap().path().is_file())
        .map(|e| e.unwrap().file_name().to_str().unwrap().to_string())
        .filter(|s| s.ends_with(".c"))
        .collect();

    // build CSparse and make it available for linking
    cc::Build::new()
        .include("deps/CSparse_modified/Include")
        .files(dot_c_files.iter().map(|name| format!("deps/CSparse_modified/Source/{0}", name)))
        .compile("CSparse");
}
