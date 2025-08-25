use bird_tool_utils::external_command_checker::*;

pub fn check_for_dependencies() {
    check_for_fastani();
}

pub fn check_for_fastani() {
    self::check_for_external_command_presence("fastaANI", "which fastANI")
        .expect("Failed to find installed fastANI");
    self::default_version_check("fastANI", "1.31", false, None)
        .expect("Failed to find sufficient version of fastANI");
}
