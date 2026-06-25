#!/usr/bin/env python3

import io
import re
from datetime import datetime
from os.path import dirname, join
import extern


def get_version(relpath):
    """Read version info from a file without importing it"""
    for line in io.open(join(dirname(__file__), relpath), encoding="utf-8"):
        if "version" in line:
            if '"' in line:
                return line.split('"')[1]
            elif "'" in line:
                return line.split("'")[1]


if __name__ == "__main__":
    yes_no = input(
        "Did you run the non-CI tests first, to make sure everything is OK (y/n)? \n\n"
        "CHECKM2DB=/work/microbiome/db/CheckM2_database/uniref100.KO.1.dmnd pixi run cargo test -- --ignored\n\n"
    )
    if yes_no != "y":
        raise Exception("Please run the non-CI tests first")

    version = get_version('../Cargo.toml')
    print(f"Version is {version}")

    print("Updating CITATION.cff ..")
    citations_lines = []
    with open("CITATION.cff", "r") as f:
        r = re.compile(r"( *version: )")
        r2 = re.compile(r"( *date-released: )")
        for line in f:
            if matches := r.match(line):
                line = matches.group(1) + version + "\n"
            elif matches := r2.match(line):
                line = matches.group(1) + datetime.today().strftime('%Y-%m-%d') + "\n"
            citations_lines.append(line)
    with open("CITATION.cff", "w") as f:
        f.writelines(citations_lines)

    print("Building docs ..")
    extern.run(f"pixi run -e dev python3 admin/build_docs.py --version {version}")

    print("Committing release files ..")
    extern.run("git add CITATION.cff docs/")
    extern.run(f'git commit -m "Release v{version}"')

    print("Checking repo is clean ..")
    extern.run('if [[ $(git diff --shortstat 2> /dev/null | tail -n1) != "" ]]; then exit 1; fi')

    extern.run(f"git tag v{version}")
    print(f"\nRelease prepared. Now run:\n\n  git push && git push --tags\n\nCI will build binaries and create the GitHub Release.")
    print("Run: `cargo publish` to crates.io.")
