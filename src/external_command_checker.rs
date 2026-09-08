use bird_tool_utils::external_command_checker::*;

/// Check for the external tools required by the chosen preclustering and
/// clustering methods. The finch preclusterer is a library, so requires nothing.
pub fn check_for_clustering_dependencies(precluster_method: &str, cluster_method: &str) {
    if precluster_method == "skani" || cluster_method == "skani" {
        check_for_skani();
    }
    if cluster_method == "fastani" {
        check_for_fastani();
    }
}

pub fn check_for_fastani() {
    self::check_for_external_command_presence_with_which("fastANI")
        .expect("Failed to find installed fastANI");
    self::default_version_check("fastANI", "1.31", false, None)
        .expect("Failed to find sufficient version of fastANI");
}

pub fn check_for_skani() {
    self::check_for_external_command_presence_with_which("skani")
        .expect("Failed to find installed skani");
    self::default_version_check("skani", "0.2.2", false, None)
        .expect("Failed to find sufficient version of skani");
}
