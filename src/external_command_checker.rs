use bird_tool_utils::external_command_checker::*;

pub fn check_for_dependencies() {
    check_for_fastani();
    check_for_dashing();
}

pub fn check_for_fastani() {
    self::check_for_external_command_presence("fastaANI", "which fastANI")
        .expect("Failed to find installed fastANI");
    self::default_version_check("fastANI", "1.31", false, None)
        .expect("Failed to find sufficient version of fastANI");
}

pub fn check_for_dashing() {
    self::check_for_external_command_presence("dashing", "which dashing")
        .expect("Failed to find installed dashing. You may wish to use the finch precluster method if you are having problems with dashing.");
    self::default_version_check("dashing", "0.4.0", true, None)
        .expect("Failed to find sufficient version of dashing. You may wish to use the finch precluster method if you are having problems with dashing.");
}
