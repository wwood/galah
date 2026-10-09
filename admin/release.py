#!/usr/bin/env python3

import argparse
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def run(cmd: list[str], *, check: bool = True) -> subprocess.CompletedProcess:
    print(f"+ {' '.join(cmd)}")
    return subprocess.run(cmd, check=check)


def capture(cmd: list[str]) -> str:
    return subprocess.check_output(cmd, text=True).strip()


def die(message: str) -> None:
    print(f"error: {message}", file=sys.stderr)
    sys.exit(1)


def validate_version(version: str) -> None:
    if not re.fullmatch(r"\d+\.\d+\.\d+(?:[-+][A-Za-z0-9.-]+)?", version):
        die(
            f"invalid version {version!r}; expected something like "
            "1.2.3, 1.2.3-beta.1, or 1.2.3+build.1"
        )


def check_clean_git() -> None:
    status = capture(["git", "status", "--porcelain"])
    if status:
        die("git working tree is not clean; commit or stash changes first")


# The package version is the first top-level `version = "..."` line, under [package]
CARGO_VERSION_RE = r'(?m)^version\s*=\s*"([^"]+)"'


def get_current_version() -> str:
    cargo_path = Path("Cargo.toml")
    if m := re.search(CARGO_VERSION_RE, cargo_path.read_text()):
        return m.group(1)
    die(f"could not find version in {cargo_path}")


def update_cargo_version(version: str) -> None:
    cargo_path = Path("Cargo.toml")
    text = cargo_path.read_text()
    new_text, count = re.subn(
        CARGO_VERSION_RE,
        f'version = "{version}"',
        text,
        count=1,
    )
    if count != 1:
        die(f"could not find a version line in {cargo_path}")
    cargo_path.write_text(new_text)
    print(f"Updated {cargo_path}: version = {version!r}")

    # Only galah's own entry in Cargo.lock changes, dependencies are left as-is
    run(["cargo", "update", "--workspace"])


def changelog_contains_version(version: str, path: Path) -> bool:
    if not path.exists():
        return False
    text = path.read_text()
    return bool(re.search(rf"(?m)^## \[{re.escape(version)}\]", text))


def update_changelog(version: str, path: Path) -> None:
    if not path.exists():
        die(f"{path} not found — create it with an '## [Unreleased]' section first")

    text = path.read_text()

    if changelog_contains_version(version, path):
        print(f"{path} already contains a section for {version}; skipping.")
        return

    unreleased_match = re.search(r"(?m)^## \[Unreleased\]\s*$", text)
    if not unreleased_match:
        die(f"could not find '## [Unreleased]' section in {path}")

    next_heading_match = re.search(r"(?m)^## \[", text[unreleased_match.end():])
    if not next_heading_match:
        die(f"could not find next '## [' release section after Unreleased in {path}")

    unreleased_start = unreleased_match.end()
    next_heading_start = unreleased_start + next_heading_match.start()
    unreleased_body = text[unreleased_start:next_heading_start].strip()

    if not unreleased_body:
        die("CHANGELOG.md [Unreleased] section is empty; add entries before releasing")

    # Leave an empty Unreleased section in place for the next release
    date_str = datetime.today().strftime("%Y-%m-%d")
    new_section = f"\n\n## [{version}] - {date_str}\n\n{unreleased_body}\n\n"

    before = text[:unreleased_start].rstrip()
    after = text[next_heading_start:].lstrip()
    path.write_text(before + new_section + after)
    print(f"Moved [Unreleased] entries to [{version}] - {date_str} in {path}")


def update_citation(version: str) -> None:
    citation_path = Path("CITATION.cff")
    if not citation_path.exists():
        print("No CITATION.cff found; skipping.")
        return
    lines = []
    r_version = re.compile(r"( *version: )")
    r_date = re.compile(r"( *date-released: )")
    for line in citation_path.read_text().splitlines(keepends=True):
        if m := r_version.match(line):
            line = m.group(1) + version + "\n"
        elif m := r_date.match(line):
            line = m.group(1) + datetime.today().strftime("%Y-%m-%d") + "\n"
        lines.append(line)
    citation_path.write_text("".join(lines))
    print(f"Updated {citation_path}: version = {version!r}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare a galah release."
    )
    parser.add_argument(
        "--version",
        help="Release version, e.g. 0.6.0 (default: read from Cargo.toml)",
    )
    parser.add_argument(
        "--no-commit",
        action="store_true",
        help="Modify files but do not commit or tag",
    )
    parser.add_argument(
        "--tag-prefix",
        default="v",
        help='Git tag prefix; default is "v", producing tags like v0.6.0',
    )
    parser.add_argument(
        "--allow-dirty",
        action="store_true",
        help="Allow running with an already-dirty git working tree",
    )

    args = parser.parse_args()

    version = args.version or get_current_version()
    validate_version(version)
    tag = f"{args.tag_prefix}{version}"

    if not args.allow_dirty:
        check_clean_git()

    existing_tags = capture(["git", "tag", "--list", tag])
    if existing_tags:
        die(f"git tag {tag!r} already exists")

    yes_no = input(
        "\nDid you run the non-CI tests first, to make sure everything is OK? [y/n]\n\n"
        "  bash tests/run_expensive_tests_at_cmr.sh\n\n"
        "or equivalently:\n\n"
        "  CHECKM2DB=/work/microbiome/db/CheckM2_database/uniref100.KO.1.dmnd "
        "ISITEUK_METAPACKAGE_PATH=/work/microbiome/db/isiteuk/isiteuk-0.0.1.smpkg "
        "EUKCC2_DB=/work/microbiome/db/eukcc/eukcc2_db_ver_1.1 "
        "pixi run cargo test -- --ignored\n\n"
    ).strip().lower()
    if yes_no not in {"y", "yes"}:
        die("run the non-CI tests first")

    print(f"Version is {version}")

    # Update version in Cargo.toml (and Cargo.lock) if a new version was specified
    current_version = get_current_version()
    if version != current_version:
        update_cargo_version(version)

    # Move Unreleased → versioned section in CHANGELOG.md
    update_changelog(version, Path("CHANGELOG.md"))

    # Update CITATION.cff
    update_citation(version)

    # Build docs
    print("\nBuilding docs ...")
    run(["pixi", "run", "-e", "dev", "python3", "admin/build_docs.py", "--version", version])

    if args.no_commit:
        print("\nStopped before commit/tag because --no-commit was supplied.")
        print("Review changes with: git diff")
        return

    print(f"\nTagging the release as {tag}")
    run(["git", "add", "Cargo.toml", "Cargo.lock", "CHANGELOG.md", "CITATION.cff", "docs/"])
    run(["git", "commit", "-m", f"Release {tag}"])
    run(["git", "tag", tag])
    print("\nNow run: git push && git push --tags")
    print("CI will build binaries and create the GitHub Release.")
    print("Then run `cargo publish` to publish to crates.io.")


if __name__ == "__main__":
    main()
