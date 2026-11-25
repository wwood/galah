#!/bin/bash
# Build and smoke-test each Galah install method in Docker.
# Usage: bash compile_and_test_install_methods.bash GALAH_VERSION   (e.g. 0.7.0)
set -e

if [ -z "$1" ]; then
    echo "Usage: $0 GALAH_VERSION (e.g. 0.7.0)" >&2
    exit 1
fi
GALAH_VERSION=$1
cd "$(dirname "$0")"

for method in pixi conda cargo; do
    echo "=== $method ==="
    if [ "$method" = pixi ]; then
        sed "s/GALAH_VERSION/$GALAH_VERSION/g" pixi.toml.in > pixi.toml
    fi
    sed "s/GALAH_VERSION/$GALAH_VERSION/g" Dockerfile.$method.in > Dockerfile.$method
    docker build --no-cache --progress=plain -f Dockerfile.$method . &> $method.build.log
    echo "$method: OK (log: $method.build.log)"
    docker system prune -af
done
