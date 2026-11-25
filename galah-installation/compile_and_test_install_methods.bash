#!/bin/bash
# Build and test each Galah install method in Docker.
#   1. build: install Galah, then smoke-test `galah --version` and `galah cluster`
#   2. run:   download reference databases (once, into ./databases) and test `galah analyse` / `process`
# Usage: bash compile_and_test_install_methods.bash GALAH_VERSION [--skip-quality]
#   e.g. bash compile_and_test_install_methods.bash 0.7.0
# Databases need roughly 10+ GB of disk and a long first download; --skip-quality skips step 2.
set -u

if [ -z "$1" ]; then
    echo "Usage: $0 GALAH_VERSION [--skip-quality] (e.g. 0.7.0)" >&2
    exit 1
fi
GALAH_VERSION=$1
SKIP_QUALITY=${2:-}
cd "$(dirname "$0")"
mkdir -p databases

FAILED=()
for method in pixi conda cargo; do
    echo "=== $method ==="
    if [ "$method" = pixi ]; then
        sed "s/GALAH_VERSION/$GALAH_VERSION/g" pixi.toml.in > pixi.toml
    fi
    sed "s/GALAH_VERSION/$GALAH_VERSION/g" Dockerfile.$method.in > Dockerfile.$method
    if docker build --no-cache --progress=plain -t galah-install-$method -f Dockerfile.$method . &> $method.build.log; then
        echo "$method: build + cluster OK (log: $method.build.log)"
        if [ "$SKIP_QUALITY" != "--skip-quality" ]; then
            if docker run --rm -v "$PWD/databases:/db" galah-install-$method \
                bash /app/test_quality.sh $method &> $method.quality.log; then
                echo "$method: analyse/process OK (log: $method.quality.log)"
            else
                echo "$method: analyse/process FAILED (log: $method.quality.log)"
                FAILED+=("$method:quality")
            fi
        fi
    else
        echo "$method: build FAILED (log: $method.build.log)"
        FAILED+=("$method:build")
    fi
    docker system prune -af
done

if [ ${#FAILED[@]} -gt 0 ]; then
    echo "FAILED: ${FAILED[*]}"
    exit 1
fi
echo "All methods passed"
