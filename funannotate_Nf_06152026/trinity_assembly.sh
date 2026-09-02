#!/bin/bash
#############################################################################
# trinity_assembly.slurm
#
# Read QC + de novo transcriptome assembly of paired-end fungal RNA-seq
# reads with Trinity, for submission to a SLURM-managed HPC cluster.
#
# Pipeline:
#   1. FastQC on raw reads
#   2. fastp adapter/quality trimming (per sample)
#   3. FastQC on trimmed reads
#   4. MultiQC report aggregating all of the above
#   5. Trinity assembly on the TRIMMED reads (pooled across samples)
#   6. Assembly stats (TrinityStats.pl)
#
# Usage:
#   1. Edit every value in the "USER-EDITABLE SETTINGS" block below.
#   2. Submit with:  sbatch trinity_assembly.slurm
#   3. Monitor with: squeue -u $USER   /   sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS
#
# Notes on the approach:
#   - Reads for all samples are QC'd individually, then pooled into ONE
#     Trinity run, which is the standard way to build a single reference
#     transcriptome across conditions/replicates. If you want a separate
#     assembly per sample, see the comment near the Trinity command.
#   - Trinity writes a very large number of small intermediate files during
#     the Inchworm/Chrysalis/Butterfly phases. On parallel filesystems
#     (Lustre/GPFS, typical for /scratch or /project on HPC) this can be
#     slow and can even trip metadata-server quotas. This script therefore
#     runs Trinity's (and fastp/FastQC's) working files on node-local
#     scratch ($SLURM_TMPDIR or $TMPDIR, whichever the cluster provides)
#     and rsyncs only the QC reports and final assembly back to persistent
#     storage at the end.
#   - This QC is generic adapter/quality trimming, appropriate for any
#     Illumina RNA-seq data. It does NOT screen for rRNA contamination; if
#     rRNA reads are a known issue for your samples, consider adding a
#     SortMeRNA step between trimming and Trinity (commented pointer left
#     below where trimmed reads are finalized).
#############################################################################

#SBATCH --job-name=trinity_fungal
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=250G
#SBATCH --output=trinity_%x_%j.out
#SBATCH --error=trinity_%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eric.morrison@unh.edu

set -euo pipefail

#############################################################################
# USER-EDITABLE SETTINGS
#############################################################################

# --- Input reads -----------------------------------------------------------
# Paired-end fastq(.gz) files. SAMPLE_NAMES, LEFT_READS and RIGHT_READS are
# comma-separated with NO spaces, and must all be the same length and in
# the same order (Nth name <-> Nth left file <-> Nth right file).
#
READS_DIR="${HOME}/neonectria_faginata_draft_genome/Nf_annotate/funannotate_09022026/Nf_rnaseq_reads"

SAMPLE_NAMES="C-A1,C-B1,C-C1,C-D1,N-E2,N-F2,N-G1,N-H1"
LEFT_READS="${READS_DIR}/C-A1_S224_L002_R1_001.fastq.gz,${READS_DIR}/C-B1_S225_L002_R1_001.fastq.gz,${READS_DIR}/C-C1_S226_L002_R1_001.fastq.gz,${READS_DIR}/C-D1_S227_L002_R1_001.fastq.gz,${READS_DIR}/N-E2_S228_L002_R1_001.fastq.gz,${READS_DIR}/N-F2_S229_L002_R1_001.fastq.gz,${READS_DIR}/N-G1_S230_L002_R1_001.fastq.gz,${READS_DIR}/N-H1_S231_L002_R1_001.fastq.gz"
RIGHT_READS="${READS_DIR}/C-A1_S224_L002_R2_001.fastq.gz,${READS_DIR}/C-B1_S225_L002_R2_001.fastq.gz,${READS_DIR}/C-C1_S226_L002_R2_001.fastq.gz,${READS_DIR}/C-D1_S227_L002_R2_001.fastq.gz,${READS_DIR}/N-E2_S228_L002_R2_001.fastq.gz,${READS_DIR}/N-F2_S229_L002_R2_001.fastq.gz,${READS_DIR}/N-G1_S230_L002_R2_001.fastq.gz,${READS_DIR}/N-H1_S231_L002_R2_001.fastq.gz"

# Strand-specific library? Set to the appropriate Trinity code
# (RF or FR for paired-end dUTP-style kits) or leave empty ("") if your
# library prep is NOT strand-specific (most standard fungal RNA-seq kits are
# not, unless you specifically used a stranded kit -- check your library prep
# protocol / sequencing facility report if unsure).
SS_LIB_TYPE=""     # e.g. "RF"

# --- QC / trimming parameters ------------------------------------------------
FASTP_QUALIFIED_QUALITY=20      # -q : phred score threshold for "qualified" bases
FASTP_MIN_LENGTH=36             # -l : discard reads shorter than this after trimming
FASTP_EXTRA_ARGS=""             # any extra fastp flags, e.g. "--trim_poly_g" (auto-detected for NextSeq/NovaSeq anyway)

# --- Output ------------------------------------------------------------------
# Defaults to a sibling of the funannotate directory that holds the reads,
# i.e. .../Nf_annotate/trinity_out and .../Nf_annotate/qc_reports. Change if
# you'd rather it land somewhere else.
PROJECT_DIR="${HOME}/neonectria_faginata_draft_genome/Nf_annotate/funannotate_09022026"
FINAL_OUTDIR="${PROJECT_DIR}/trinity_out"          # final assembly lands here, on persistent storage
QC_OUTDIR="${PROJECT_DIR}/rnaseq_qc_reports"               # FastQC/fastp/MultiQC reports land here, on persistent storage

# --- Resources (must match the #SBATCH lines above) --------------------------
NUM_CPUS=${SLURM_CPUS_PER_TASK:-24}
MAX_MEM_GB=250                                     # keep a little headroom below --mem for OS/Java overhead

# --- Software environment -----------------------------------------------------
# This cluster splits the toolchain across a linuxbrew module and two conda
# envs reached through the anaconda module:
#   - Trinity + FastQC:  `module load linuxbrew/colsa`
#   - fastp:              conda env `fastp-1.0.1`, via `module load anaconda/colsa`
#   - MultiQC:             conda env `multiqc-1.10.1`, via the same anaconda module
#
# linuxbrew/colsa is loaded once, up front, and stays loaded for the whole
# job -- Trinity and FastQC are called directly throughout the script. The
# anaconda module is also loaded up front (it just provides the `conda`/
# `activate` machinery), but the fastp and multiqc conda envs are activated
# only around the specific steps that need them, and deactivated right after
# -- conda envs are not additive like modules, so only one is "active" PATH-
# wise at a time, and mixing them can shadow each other's dependencies.
module purge
module load linuxbrew/colsa
#module load anaconda/colsa


# Conda's `activate`/`deactivate` shell functions sometimes reference unset
# variables, which trips `set -u` above -- temporarily relax it around calls
# to them.
conda_activate()   { set +u; conda activate "$1";   set -u; }
conda_deactivate() { set +u; conda deactivate;      set -u; }

#############################################################################
# SCRATCH SPACE SETUP
#############################################################################

# Prefer a SLURM-provided node-local scratch dir if the cluster sets one;
# fall back to $TMPDIR, then to a subdirectory under FINAL_OUTDIR as a last
# resort. Check your cluster's documentation for the correct variable
# (common ones: $SLURM_TMPDIR, $LOCAL_SCRATCH, $TMPDIR).
if [ -n "${SLURM_TMPDIR:-}" ]; then
    SCRATCH_DIR="${SLURM_TMPDIR}/trinity_work_${SLURM_JOB_ID}"
elif [ -n "${TMPDIR:-}" ]; then
    SCRATCH_DIR="${TMPDIR}/trinity_work_${SLURM_JOB_ID}"
else
    SCRATCH_DIR="${FINAL_OUTDIR}/trinity_work_${SLURM_JOB_ID}"
fi

RAW_QC_DIR="${SCRATCH_DIR}/fastqc_raw"
TRIMMED_DIR="${SCRATCH_DIR}/trimmed"
TRIMMED_QC_DIR="${SCRATCH_DIR}/fastqc_trimmed"
FASTP_REPORT_DIR="${SCRATCH_DIR}/fastp_reports"

mkdir -p "${SCRATCH_DIR}" "${FINAL_OUTDIR}" "${QC_OUTDIR}" \
         "${RAW_QC_DIR}" "${TRIMMED_DIR}" "${TRIMMED_QC_DIR}" "${FASTP_REPORT_DIR}"

# Clean up scratch space on exit, whether the job succeeds, fails, or is
# cancelled -- avoids leaving large temp directories behind on the node.
cleanup() {
    echo "[$(date)] Cleaning up scratch directory: ${SCRATCH_DIR}"
    rm -rf "${SCRATCH_DIR}"
}
trap cleanup EXIT

#############################################################################
# JOB INFO (useful for the log file / troubleshooting)
#############################################################################
echo "=========================================="
echo "Job started:      $(date)"
echo "Job ID:            ${SLURM_JOB_ID:-N/A}"
echo "Node(s):            ${SLURM_JOB_NODELIST:-N/A}"
echo "CPUs allocated:    ${NUM_CPUS}"
echo "Memory requested:  ${MAX_MEM_GB}G"
echo "Scratch directory: ${SCRATCH_DIR}"
echo "QC output dir:     ${QC_OUTDIR}"
echo "Final output dir:  ${FINAL_OUTDIR}"
echo "Trinity version:   $(Trinity --version 2>&1 | head -n1)"
echo "FastQC version:    $(fastqc --version 2>&1 | head -n1)"
echo "=========================================="

# Parse the comma-separated sample lists into bash arrays.
IFS=',' read -ra SAMPLE_ARR  <<< "${SAMPLE_NAMES}"
IFS=',' read -ra LEFT_ARR    <<< "${LEFT_READS}"
IFS=',' read -ra RIGHT_ARR   <<< "${RIGHT_READS}"

if [ "${#SAMPLE_ARR[@]}" -ne "${#LEFT_ARR[@]}" ] || [ "${#SAMPLE_ARR[@]}" -ne "${#RIGHT_ARR[@]}" ]; then
    echo "ERROR: SAMPLE_NAMES, LEFT_READS, and RIGHT_READS must have the same number of comma-separated entries." >&2
    exit 1
fi

#############################################################################
# STEP 1: FASTQC ON RAW READS
#############################################################################
echo "[$(date)] Running FastQC on raw reads..."
fastqc --threads "${NUM_CPUS}" --outdir "${RAW_QC_DIR}" "${LEFT_ARR[@]}" "${RIGHT_ARR[@]}"

#############################################################################
# STEP 2: FASTP TRIMMING (per sample)
#############################################################################
# fastp handles adapter detection/removal, quality trimming, and low-
# complexity/short-read filtering in one pass, and writes its own per-sample
# HTML/JSON QC report -- useful in itself and also consumed by MultiQC below.
TRIMMED_LEFT_ARR=()
TRIMMED_RIGHT_ARR=()

module purge
module load anaconda/colsa
conda_activate fastp-1.0.1
echo "fastp version:     $(fastp --version 2>&1 | head -n1)"

for i in "${!SAMPLE_ARR[@]}"; do
    sample="${SAMPLE_ARR[$i]}"
    left_in="${LEFT_ARR[$i]}"
    right_in="${RIGHT_ARR[$i]}"
    left_out="${TRIMMED_DIR}/${sample}_R1.trimmed.fastq.gz"
    right_out="${TRIMMED_DIR}/${sample}_R2.trimmed.fastq.gz"

    echo "[$(date)] fastp: ${sample}"
    fastp \
        --in1 "${left_in}" --in2 "${right_in}" \
        --out1 "${left_out}" --out2 "${right_out}" \
        --thread "${NUM_CPUS}" \
        --qualified_quality_phred "${FASTP_QUALIFIED_QUALITY}" \
        --length_required "${FASTP_MIN_LENGTH}" \
        --detect_adapter_for_pe \
        --correction \
        ${FASTP_EXTRA_ARGS} \
        --json "${FASTP_REPORT_DIR}/${sample}.fastp.json" \
        --html "${FASTP_REPORT_DIR}/${sample}.fastp.html"

    TRIMMED_LEFT_ARR+=("${left_out}")
    TRIMMED_RIGHT_ARR+=("${right_out}")
done

conda_deactivate
module purge
module load linuxbrew/colsa

# Optional: screen out rRNA reads here with SortMeRNA before assembly, e.g.
#   sortmerna --ref <rRNA_db.fasta> --reads "${left_out}" --reads "${right_out}" \
#             --paired_in --fastx --other "${TRIMMED_DIR}/${sample}_norRNA" ...
# and point TRIMMED_LEFT_ARR/TRIMMED_RIGHT_ARR at the --other output instead.
# Not enabled by default since it requires a pre-built rRNA reference DB.

# Build comma-separated lists of trimmed reads for Trinity.
TRIMMED_LEFT_READS=$(IFS=,; echo "${TRIMMED_LEFT_ARR[*]}")
TRIMMED_RIGHT_READS=$(IFS=,; echo "${TRIMMED_RIGHT_ARR[*]}")

#############################################################################
# STEP 3: FASTQC ON TRIMMED READS
#############################################################################
echo "[$(date)] Running FastQC on trimmed reads..."
fastqc --threads "${NUM_CPUS}" --outdir "${TRIMMED_QC_DIR}" "${TRIMMED_LEFT_ARR[@]}" "${TRIMMED_RIGHT_ARR[@]}"

#############################################################################
# STEP 4: MULTIQC SUMMARY
#############################################################################
echo "[$(date)] Aggregating QC reports with MultiQC..."
module purge
module load anaconda/colsa
conda_activate multiqc-1.10.1
echo "MultiQC version:   $(multiqc --version 2>&1 | head -n1)"
multiqc "${RAW_QC_DIR}" "${TRIMMED_QC_DIR}" "${FASTP_REPORT_DIR}" \
    --outdir "${SCRATCH_DIR}/multiqc" --filename multiqc_report.html || true
conda_deactivate
module purge
module load linuxbrew/colsa
#############################################################################
# RUN TRINITY (on trimmed reads)
#############################################################################
# For a SEPARATE assembly per sample instead of one pooled assembly, replace
# this single call with a loop over TRIMMED_LEFT_ARR/TRIMMED_RIGHT_ARR, e.g.:
#   for i in "${!SAMPLE_ARR[@]}"; do
#       Trinity --seqType fq --max_memory ${MAX_MEM_GB}G --CPU ${NUM_CPUS} \
#           --left "${TRIMMED_LEFT_ARR[$i]}" --right "${TRIMMED_RIGHT_ARR[$i]}" \
#           --output "${SCRATCH_DIR}/trinity_${SAMPLE_ARR[$i]}" --full_cleanup
#   done

TRINITY_CMD=(Trinity
    --seqType fq
    --max_memory "${MAX_MEM_GB}G"
    --CPU "${NUM_CPUS}"
    --left "${TRIMMED_LEFT_READS}"
    --right "${TRIMMED_RIGHT_READS}"
    --output "${SCRATCH_DIR}/trinity_out"
    --full_cleanup
)

if [ -n "${SS_LIB_TYPE}" ]; then
    TRINITY_CMD+=(--SS_lib_type "${SS_LIB_TYPE}")
fi

echo "[$(date)] Running: ${TRINITY_CMD[*]}"
"${TRINITY_CMD[@]}"

echo "[$(date)] Trinity assembly finished."

#############################################################################
# ASSEMBLY STATS
#############################################################################
# --full_cleanup above renames the output to <output_dir>.Trinity.fasta and
# removes the bulky intermediate working directory automatically.
ASSEMBLY_FASTA="${SCRATCH_DIR}/trinity_out.Trinity.fasta"

if [ -f "${ASSEMBLY_FASTA}" ]; then
    echo "[$(date)] Generating assembly stats..."
    # TrinityStats.pl ships alongside Trinity's util/ scripts. It's usually
    # on PATH once the Trinity module is loaded, but linuxbrew-style
    # installs sometimes symlink only the main `Trinity` binary -- fall back
    # to locating it relative to that binary if the plain command isn't found.
    if command -v TrinityStats.pl &>/dev/null; then
        STATS_CMD="TrinityStats.pl"
    else
        TRINITY_HOME="$(dirname "$(readlink -f "$(command -v Trinity)")")"
        STATS_CMD="${TRINITY_HOME}/util/TrinityStats.pl"
    fi

    if [ -x "${STATS_CMD}" ] || command -v "${STATS_CMD}" &>/dev/null; then
        "${STATS_CMD}" "${ASSEMBLY_FASTA}" > "${SCRATCH_DIR}/trinity_out.Trinity.fasta.stats.txt" 2>&1 || true
    else
        echo "WARNING: could not locate TrinityStats.pl (tried PATH and ${TRINITY_HOME}/util/) -- skipping stats." >&2
    fi
else
    echo "WARNING: expected assembly file not found at ${ASSEMBLY_FASTA}" >&2
fi

#############################################################################
# COPY RESULTS BACK TO PERSISTENT STORAGE
#############################################################################
echo "[$(date)] Copying QC reports to ${QC_OUTDIR} ..."
rsync -avh "${RAW_QC_DIR}/"      "${QC_OUTDIR}/fastqc_raw/"      2>/dev/null || true
rsync -avh "${TRIMMED_QC_DIR}/"  "${QC_OUTDIR}/fastqc_trimmed/"  2>/dev/null || true
rsync -avh "${FASTP_REPORT_DIR}/" "${QC_OUTDIR}/fastp/"          2>/dev/null || true
rsync -avh "${SCRATCH_DIR}/multiqc/" "${QC_OUTDIR}/multiqc/"     2>/dev/null || true

echo "[$(date)] Copying assembly results to ${FINAL_OUTDIR} ..."
rsync -avh "${SCRATCH_DIR}/trinity_out.Trinity.fasta"                "${FINAL_OUTDIR}/" 2>/dev/null || true
rsync -avh "${SCRATCH_DIR}/trinity_out.Trinity.fasta.gene_trans_map"  "${FINAL_OUTDIR}/" 2>/dev/null || true
rsync -avh "${SCRATCH_DIR}/trinity_out.Trinity.fasta.stats.txt"       "${FINAL_OUTDIR}/" 2>/dev/null || true

echo "=========================================="
echo "Job finished: $(date)"
echo "QC reports:     ${QC_OUTDIR}/multiqc/multiqc_report.html"
echo "Final assembly: ${FINAL_OUTDIR}/trinity_out.Trinity.fasta"
echo "=========================================="

# (Scratch cleanup happens automatically via the trap above.)
