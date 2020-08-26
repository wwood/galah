#!/bin/bash -e

set -o pipefail

export VERSION=`cargo run -- --version |awk '{print $2}'`

# For dashing which does not currently work on desktop via conda.
export PATH=~/bioinfo/dashing:$PATH

echo "Found version $VERSION .."

echo "Building normally .."
cargo build --release

echo "Testing release version .."
cargo test --release

echo "Building musl static binary .."
cargo build --target x86_64-unknown-linux-musl --release

echo "Making static dist .."
mkdir -p dist/galah-x86_64-unknown-linux-musl-$VERSION
cp \
 target/x86_64-unknown-linux-musl/release/galah \
 dist/galah-x86_64-unknown-linux-musl-$VERSION/
cd dist
tar czf galah-x86_64-unknown-linux-musl-$VERSION.tar.gz galah-x86_64-unknown-linux-musl-$VERSION
cd ..

echo "Building HTML versions of man pages .."
for SUBCOMMAND in cluster
do
    cargo run -- $SUBCOMMAND --full-help-roff |pandoc - -t markdown -f man |sed 's/\\\[/[/g; s/\\\]/]/g' |cat <(sed s/SUBCOMMAND/$SUBCOMMAND/ prelude) - >docs/galah-$SUBCOMMAND.Rmd
    echo "library(prettydoc); setwd('docs'); rmarkdown::render('galah-$SUBCOMMAND.Rmd','prettydoc::html_pretty','galah-$SUBCOMMAND.html')" |R --no-save </dev/stdin
    rm docs/galah-$SUBCOMMAND.Rmd
done

echo "Now make sure git is up to date, and run cargo publish"
