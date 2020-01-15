pub mod ani_correction;
pub mod minhash_clusterer;

#[macro_use]
extern crate log;
extern crate finch;
extern crate rayon;
extern crate partitions;

use std::io::Read;

pub fn finish_command_safely(
    mut process: std::process::Child, process_name: &str)
-> std::process::Child {
    let es = process.wait()
        .expect(&format!("Failed to glean exitstatus from failing {} process", process_name));
    debug!("Process {} finished", process_name);
    if !es.success() {
        error!("Error when running {} process.", process_name);
        let mut err = String::new();
        process.stderr.expect(&format!("Failed to grab stderr from failed {} process", process_name))
            .read_to_string(&mut err).expect("Failed to read stderr into string");
        error!("The STDERR was: {:?}", err);
        let mut out = String::new();
        process.stdout.expect(&format!("Failed to grab stdout from failed {} process", process_name))
            .read_to_string(&mut out).expect("Failed to read stdout into string");
        error!("The STDOUT was: {:?}", out);
        error!("Cannot continue after {} failed.", process_name);
        std::process::exit(1);
    }
    return process;
}
