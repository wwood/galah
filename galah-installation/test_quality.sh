#!/bin/bash
# Download reference databases (if absent) and test `galah analyse` and `galah process`.
# Runs inside an install-method container. Expects:
#   GALAH_CMD  command prefix to run galah (e.g. "galah", "pixi run galah")
#   PIXI_CMD   command to run pixi (used to fetch checkm2/isiteuk just to download their databases)
#   /db        (mounted) database directory, persisted between runs/methods
# Usage: bash test_quality.sh METHOD
set -euo pipefail

METHOD=${1:?Usage: $0 METHOD}
DB=/db
DATA=/app/galah_src/tests/data/domain_examples
OUT=/app/test_quality_$METHOD
mkdir -p "$DB" "$OUT"

# --- Databases (large; downloaded once into the mounted /db) ---
if ! find "$DB/checkm2" -name 'uniref100.KO.1.dmnd' 2>/dev/null | grep -q .; then
    $PIXI_CMD exec -c conda-forge -c bioconda -s checkm2 checkm2 database --download --path "$DB/checkm2"
fi
CHECKM2DB=$(find "$DB/checkm2" -name 'uniref100.KO.1.dmnd' | head -n1)

if ! ls -d "$DB"/isiteuk/isiteuk-* >/dev/null 2>&1; then
    $PIXI_CMD exec -c conda-forge -c bioconda -s isiteuk isiteuk data --output-directory "$DB/isiteuk"
fi
ISITEUK_METAPACKAGE_PATH=$(ls -d "$DB"/isiteuk/isiteuk-* | head -n1)

if [ ! -d "$DB/eukcc/eukcc2_db_ver_1.1" ]; then
    mkdir -p "$DB/eukcc"
    wget -q -O "$DB/eukcc/eukcc2_db_ver_1.1.tar.gz" \
        http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.1.tar.gz
    tar -xzf "$DB/eukcc/eukcc2_db_ver_1.1.tar.gz" -C "$DB/eukcc"
    rm "$DB/eukcc/eukcc2_db_ver_1.1.tar.gz"
fi
EUKCC2_DB="$DB/eukcc/eukcc2_db_ver_1.1"
export CHECKM2DB ISITEUK_METAPACKAGE_PATH EUKCC2_DB
echo "CHECKM2DB=$CHECKM2DB ISITEUK_METAPACKAGE_PATH=$ISITEUK_METAPACKAGE_PATH EUKCC2_DB=$EUKCC2_DB"

GENOMES="$DATA/GCF_002008365.1_genomic.fna.gz $DATA/GCA_003139855.1_genomic.fna.gz $DATA/binchicken_co8412.34_euk.fna.gz"

check_mimag() {
    local f=$1
    [ "$(wc -l < "$f")" -eq 4 ] || { echo "FAIL: $f should have header + 3 genomes"; exit 1; }
    grep -P 'GCF_002008365.1\S*\tBacteria\t' "$f" >/dev/null || { echo "FAIL: GCF_002008365.1 not Bacteria"; exit 1; }
    grep -P 'GCA_003139855.1\S*\tArchaea\t' "$f" >/dev/null || { echo "FAIL: GCA_003139855.1 not Archaea"; exit 1; }
    grep -P 'binchicken_co8412\S*\tEukaryota\t' "$f" >/dev/null || { echo "FAIL: euk genome not Eukaryota"; exit 1; }
}

# --- analyse ---
$GALAH_CMD analyse --genome-fasta-files $GENOMES --output-mimag-summary "$OUT/analyse_mimag.tsv"
cat "$OUT/analyse_mimag.tsv"
check_mimag "$OUT/analyse_mimag.tsv"

# --- process (analyse + cluster) ---
$GALAH_CMD process --genome-fasta-files $GENOMES \
    --output-mimag-summary "$OUT/process_mimag.tsv" \
    --output-cluster-definition "$OUT/process_clusters.tsv"
cat "$OUT/process_mimag.tsv" "$OUT/process_clusters.tsv"
check_mimag "$OUT/process_mimag.tsv"

echo "$METHOD: analyse/process OK"
