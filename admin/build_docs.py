#!/usr/bin/env python3

import extern
import logging
import argparse
import io
from os.path import dirname, join
import os


def remove_before(marker, string_to_process):
    splitter = "\n# " + marker + "\n"
    if splitter not in string_to_process:
        raise Exception(f"Marker '{marker}' not found in string")
    return splitter + string_to_process.split(splitter)[1]


def get_version(relpath):
    """Read version info from a file without importing it"""
    for line in io.open(join(dirname(__file__), relpath), encoding="cp437"):
        if "__version__" in line:
            if '"' in line:
                return line.split('"')[1]
            elif "'" in line:
                return line.split("'")[1]


if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")
    parent_parser.add_argument('--version', help='not with v e.g. 0.12.11', required=True)

    args = parent_parser.parse_args()

    # Setup logging
    debug = True
    if args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.DEBUG
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    subdir_and_commands = [
        ["tools", ["cluster", "analyse", "process"]],
    ]

    for subdir, commands in subdir_and_commands:
        for subcommand in commands:
            cmd_stub = f"cargo run -- {subcommand} --full-help-roff |pandoc - -t markdown-multiline_tables-simple_tables-grid_tables -f man |sed 's/\\\\\\[/[/g; s/\\\\\\]/]/g; s/^: //'"
            logging.debug("May need to install pandoc (e.g. conda install pandoc)")
            man_usage = extern.run(cmd_stub)

            subcommand_prelude = f"docs/preludes/{subcommand}_prelude.md"
            if os.path.exists(subcommand_prelude):
                # Remove everything before the options section
                splitters = {
                    "cluster": "GENOME INPUT",
                    "analyse": "GENOME INPUT",
                    "process": "GENOME INPUT",
                }
                logging.info(f"For ROFF for command {subcommand}, removing everything before '{splitters[subcommand]}'")
                man_usage = remove_before(splitters[subcommand], man_usage)

                with open(f"docs/{subdir}/{subcommand}.md",'w') as f:
                    f.write("---\n")
                    f.write(f"title: Galah {subcommand}\n")
                    f.write("---\n")
                    f.write(f"# galah {subcommand}\n")

                    with open(subcommand_prelude) as f2:
                        f.write(f2.read())

                    f.write(man_usage)
            else:
                man_usage = remove_before("DESCRIPTION", man_usage)
                with open("docs/{}/{}.md".format(subdir, subcommand),'w') as f:
                    f.write("---\n")
                    f.write(f"title: Galah {subcommand}\n")
                    f.write("---\n")
                    f.write(f"# galah {subcommand}\n")

                    f.write(man_usage)

    extern.run("doctave build")
