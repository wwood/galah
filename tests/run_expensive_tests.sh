#!/bin/bash

export CHECKM2DB=/work/microbiome/db/CheckM2_database/uniref100.KO.1.dmnd
export ISITEUK_METAPACKAGE_PATH=/work/microbiome/db/isiteuk/isiteuk-0.0.1.smpkg
export EUKCC2_DB=/work/microbiome/db/eukcc/eukcc2_db_ver_1.1

pixi run cargo test --no-fail-fast -- --ignored 2>&1
