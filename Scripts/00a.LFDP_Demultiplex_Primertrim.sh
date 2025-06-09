#!/bin/bash

set -e  # Exit if any command fails

# ======= USER INPUT ========
BASE_DIR="/Users/glennd/Documents/GitHub/eDNA-LFDP/Raw_data"
BARCODE_DATA="/Users/glennd/Documents/GitHub/eDNA-LFDP/Raw_data/LFDP_barcode_data.txt"
TAR_FILE="X204SC24022146-Z01-F001.tar"
# ============================

echo ">>> Navigating to base directory: $BASE_DIR"
cd "$BASE_DIR"

echo ">>> Extracting TAR file: $TAR_FILE"
tar -xf "$TAR_FILE"

DATA_DIR="$BASE_DIR/X204SC24022146-Z01-F001/01.RawData"
cd "$DATA_DIR"

echo ">>> Creating directories for plate1 and plate2"
mkdir -p plate1 plate2

echo ">>> Moving raw data into plate directories"
mv *P1 plate1/
mv *P2 plate2/

# Function for processing each plate
process_plate () {
    local PLATE=$1
    echo ">>> Processing $PLATE"

    cd "$DATA_DIR/$PLATE"
    echo ">>> Flattening *.fq.gz files into $PLATE directory"
    mv **/*.fq.gz . || true  # Avoid error if already moved

    echo ">>> Running sabre for $PLATE"
    cat <<EOF > sabre_${PLATE}.sh
for i in *1.fq.gz; do
    bn=\${i/1.fq.gz}
    sabre pe -f \${bn}1.fq.gz -r \${bn}2.fq.gz -b "$BARCODE_DATA" -u \${bn}unassigned1.fq -w \${bn}unassigned.fq
    for r in {1..4}; do
        mv rep\${r}f \${bn}_rep\${r}f.fq
        mv rep\${r}r \${bn}_rep\${r}r.fq
    done
done
EOF
    bash sabre_${PLATE}.sh

    echo ">>> Running cutadapt for $PLATE"
    for rep in 1 2 3 4; do
        case $rep in
            1) rbarcode="AGGAA"; rcRv="GATAGGTGCAGAGACTCAATGGTTCCT"; Rv="AGGAACCATTGAGTCTCTGCACCTATC";;
            2) rbarcode="GAGTGG"; rcRv="GATAGGTGCAGAGACTCAATGGCCACTC"; Rv="GAGTGGCCATTGAGTCTCTGCACCTATC";;
            3) rbarcode="CCACGTC"; rcRv="GATAGGTGCAGAGACTCAATGGGACGTGG"; Rv="CCACGTCCCATTGAGTCTCTGCACCTATC";;
            4) rbarcode="TTCTCAGC"; rcRv="GATAGGTGCAGAGACTCAATGGGCTGAGAA"; Rv="TTCTCAGCCCATTGAGTCTCTGCACCTATC";;
        esac

        echo ">>> Running cutadapt for $PLATE replicate $rep"
        cat <<EOF > ${PLATE}_rep${rep}.sh
for i in *rep${rep}f.fq; do
    bn=\${i/rep${rep}f.fq}
    cutadapt -a ^GGGCAATCCTGAGCCAA...$rcRv -A $Rv...TTGGCTCAGGATTGCCC \\
        --untrimmed-output \${bn}.rep${rep}out1.fq.gz \\
        --untrimmed-paired-output \${bn}.rep${rep}out2.fq.gz \\
        -o \${bn}.rep${rep}.trim1.fq.gz -p \${bn}.rep${rep}.trim2.fq.gz \\
        \${bn}rep${rep}f.fq \${bn}rep${rep}r.fq
done
EOF
        bash ${PLATE}_rep${rep}.sh
    done

    echo ">>> Moving trimmed files for $PLATE to ${PLATE}trimmed/"
    mkdir -p ${PLATE}trimmed
    mv *trim1.fq.gz *trim2.fq.gz ${PLATE}trimmed/

    echo ">>> Cleaning up intermediate files for $PLATE"
    rm -f *_rep*f.fq *_rep*r.fq *unassigned*.fq *.out*.fq.gz sabre_${PLATE}.sh ${PLATE}_rep*.sh

    cd "$DATA_DIR"
}

process_plate "plate1"
process_plate "plate2"

echo ">>> Cleanup complete. Final files: raw fastq + trimmed outputs."
echo ">>> Pipeline complete! Ready for R."
