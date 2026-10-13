#!/bin/bash

# Regenerates the simulated datasets used by the simulated validation tests (./gradlew simulatedTest).
# Requires python3 with pysam, bwa and samtools (e.g. conda activate genomics_analysis).
# The outputs are committed under src/test/resources/simulated, so this script only needs to be run when the
# simulation itself changes. After regenerating, the baselines must be refreshed:
#   ./gradlew simulatedTest -PupdateBaseline

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESOURCES_DIR="${SCRIPT_DIR}/../resources/simulated"
WELLFORMED_DIR="${RESOURCES_DIR}/wellformed"
MALFORMED_DIR="${RESOURCES_DIR}/malformed"
SEED="${1:-20261009}"

for tool in python3 bwa samtools; do
    if ! command -v ${tool} &> /dev/null; then
        echo "Error: ${tool} is not installed or not in PATH"
        exit 1
    fi
done

WORK_DIR="$(mktemp -d)"
trap 'rm -rf "${WORK_DIR}"' EXIT

echo "Simulating reads with seed ${SEED}..."
python3 "${SCRIPT_DIR}/simulate_tumor_normal.py" -o "${WORK_DIR}" -s "${SEED}"

echo "Aligning reads..."
# Relative paths are used so the program records of the BAM headers do not change between executions
pushd "${WORK_DIR}" > /dev/null
bwa index reference.fa 2> /dev/null
samtools faidx reference.fa
for sample in normal tumor; do
    suffix="N"
    if [ "${sample}" == "tumor" ]; then suffix="T"; fi
    read_group="@RG\tID:HSIM${suffix}.1\tSM:PATIENT_0042_${suffix}\tLB:PATIENT_0042_${suffix}_lib1\tPL:ILLUMINA\tPU:HSIM${suffix}.1"
    bwa mem -t 1 -R "${read_group}" reference.fa "${sample}_R1.fastq" "${sample}_R2.fastq" 2> /dev/null \
        | samtools fixmate -m - - \
        | samtools sort -o "${sample}.sorted.bam" -
    samtools markdup "${sample}.sorted.bam" "${sample}.bam" 2> /dev/null
    samtools index "${sample}.bam"
done
popd > /dev/null

rm -rf "${WELLFORMED_DIR}" "${MALFORMED_DIR}"
mkdir -p "${WELLFORMED_DIR}" "${MALFORMED_DIR}"
for file in reference.fa reference.fa.fai normal.bam normal.bam.bai tumor.bam tumor.bam.bai truth_variants.tsv regions.bed; do
    cp "${WORK_DIR}/${file}" "${WELLFORMED_DIR}/${file}"
done

echo "Building malformed datasets..."
python3 "${SCRIPT_DIR}/make_malformed.py" -i "${WELLFORMED_DIR}" -o "${MALFORMED_DIR}"

echo "Datasets written to ${RESOURCES_DIR}"
for sample in normal tumor; do
    echo "${sample}: $(samtools view -c "${WELLFORMED_DIR}/${sample}.bam") alignments, mean depth" \
        "$(samtools depth -a "${WELLFORMED_DIR}/${sample}.bam" | awk '{s+=$3} END {printf "%.1f", s/NR}')"
done
