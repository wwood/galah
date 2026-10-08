#!/bin/bash
# Submit expensive tests to the PBS queue.
# Run from the repo root: bash tests/run_expensive_tests_at_cmr.sh
mqsub --no-email -t 8 -m 64 --hours 2 --segregated-log-files -- \
    bash "$(dirname "$0")/run_expensive_tests.sh"
