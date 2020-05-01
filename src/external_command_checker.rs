use bird_tool_utils::external_command_checker::*;

pub fn check_for_dependencies() {
    check_for_fastani();
    check_for_dashing();
}

pub fn check_for_fastani() {
    self::check_for_external_command_presence("fastaANI", "which fastANI");
    self::default_version_check("fastANI", "1.3", false, None);
}

pub fn check_for_dashing() {
    self::check_for_external_command_presence("dashing", "which dashing");
    self::default_version_check("dashing", "0.4.0", true, None);
}
