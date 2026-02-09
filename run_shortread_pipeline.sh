#!/bin/bash
set -euo pipefail

# ==============================
# Defaults: tools, refs, knobs
# ==============================
DM6_REF="/rd/caiya/dm6.fa"
THREADS_ALIGN=32
THREADS_SORT=16
THREADS_SAMTOOLS=32
THREADS_DEPTH=32

# Create symlinks with provided short names in workdir
LINK_FASTQ="yes"

# Toggle heavy steps
DO_ALIGN="yes"
DO_DEPTH="yes"
DO_CHIMERA="yes"

# Concurrency cap for background chimera Python jobs
CHIMERA_JOBS=4

# ==============================
# Usage
# ==============================
usage() {
  cat <<'EOF'
Usage:
  run_shortread_pipeline.sh \
    --workdir /path/to/workdir \
    [--link-fastq yes|no] \
    [--do-align yes|no] \
    [--do-depth yes|no] \
    [--do-chimera yes|no] \
    [--dm6-ref /path/to/dm6.fa] \
    [--chimera-jobs N] \
    --sample \
      --fq1-abs /abs/path/to/R1.fq.gz --fq1-short Short_R1.fq.gz \
      --fq2-abs /abs/path/to/R2.fq.gz --fq2-short Short_R2.fq.gz \
      --label SampleLabel \
    [--sample ...]  # repeat for more samples

Notes:
- --label is the logical sample name (e.g., EQ-42). Determines BAM prefixes and summary labels.
- If --link-fastq=yes (default), symlinks are created in workdir with given short names.
- Base counts always run; alignment/depth/chimera run only if the respective flags are yes.

Examples:

1) Total bases only:
  run_shortread_pipeline.sh \
    --workdir /rd/caiya/wmx/MDA_20251223 \
    --do-align no \
    --sample \
      --fq1-abs /rd/caiya/wmx/MDA_20251223/.../RawData/EQ-42/E251218002_L01_EQ-42_1.fq.gz \
      --fq1-short EQ-42_1.fq.gz \
      --fq2-abs /rd/caiya/wmx/MDA_20251223/.../RawData/EQ-42/E251218002_L01_EQ-42_2.fq.gz \
      --fq2-short EQ-42_2.fq.gz \
      --label EQ-42

2) Full run with chimera in background (4 concurrent):
  run_shortread_pipeline.sh \
    --workdir /rd/caiya/wmx/MDA_20251223 \
    --do-align yes --do-depth yes --do-chimera yes --chimera-jobs 4 \
    --sample \
      --fq1-abs /.../EQ-42_1.fq.gz --fq1-short EQ-42_1.fq.gz \
      --fq2-abs /.../EQ-42_2.fq.gz --fq2-short EQ-42_2.fq.gz \
      --label EQ-42 \
    --sample \
      --fq1-abs /.../NEB-phi_1.fq.gz --fq1-short NEB-phi_1.fq.gz \
      --fq2-abs /.../NEB-phi_2.fq.gz --fq2-short NEB-phi_2.fq.gz \
      --label NEB-phi
EOF
}

# ==============================
# Helpers
# ==============================

need_tool() {
  command -v "$1" >/dev/null 2>&1 || { echo "Error: Required tool '$1' not found in PATH." >&2; exit 1; }
}

# Safe append using flock to avoid interleaved lines
flock_append() {
  # Usage: flock_append /path/to/file "text to append"
  local file="$1"; shift
  {
    flock -x 9
    printf "%s\n" "$*"
  } 9>>"$file"
}

# Gate to limit number of background jobs
wait_for_bg_slots() {
  local max_jobs="$1"
  while (( $(jobs -rp | wc -l) >= max_jobs )); do
    sleep 1
  done
}

# Compute total bases in a gz FASTQ via pigz + awk
total_bases() {
  local label="$1" fq="$2" mate="$3"
  pigz -dc "$fq" | awk -v sample="$label" -v mate="$mate" '
    BEGIN{all_base=0}
    NR%4==2{all_base += length($0)}
    END{printf("%s\t%s\t%d\n", sample, mate, all_base)}
  '
}

# ==============================
# Parse arguments
# ==============================

parse_sample_block() {
  local fq1_abs="" fq1_short="" fq2_abs="" fq2_short="" label=""
  while (( $# )); do
    case "$1" in
      --fq1-abs) shift; fq1_abs="${1:?--fq1-abs requires a value}";;
      --fq1-short) shift; fq1_short="${1:?--fq1-short requires a value}";;
      --fq2-abs) shift; fq2_abs="${1:?--fq2-abs requires a value}";;
      --fq2-short) shift; fq2_short="${1:?--fq2-short requires a value}";;
      --label) shift; label="${1:?--label requires a value}";;
      --sample) break ;;
      --*) echo "Unknown sample option: $1" >&2; exit 1;;
      *) break;;
    esac
    shift
  done
  [[ -n "$fq1_abs" && -n "$fq1_short" && -n "$fq2_abs" && -n "$fq2_short" && -n "$label" ]] || {
    echo "Error: Incomplete --sample block." >&2; exit 1;
  }
  S_FQ1_ABS+=("$fq1_abs"); S_FQ1_SHORT+=("$fq1_short")
  S_FQ2_ABS+=("$fq2_abs"); S_FQ2_SHORT+=("$fq2_short")
  S_LABEL+=("$label")
  printf '%s\0' "$@"
}



WORKDIR=""
LINK_FASTQ="${LINK_FASTQ:-yes}"
DO_ALIGN="${DO_ALIGN:-yes}"
DO_DEPTH="${DO_DEPTH:-yes}"
DO_CHIMERA="${DO_CHIMERA:-yes}"
CHIMERA_JOBS="${CHIMERA_JOBS:-4}"
DM6_REF="${DM6_REF:-/rd/caiya/dm6.fa}"

declare -a S_FQ1_ABS S_FQ1_SHORT S_FQ2_ABS S_FQ2_SHORT S_LABEL


i=1
argc=$#
while (( i <= argc )); do
  arg="${!i}"
  case "$arg" in
    -h|--help)usage; exit 0;;
    --workdir)      ((i++)); (( i <= argc )) || { echo "Error: --workdir requires a value" >&2; exit 1; }; WORKDIR="${!i}" ;;
    --link-fastq)   ((i++)); (( i <= argc )) || { echo "Error: --link-fastq requires yes|no" >&2; exit 1; }; LINK_FASTQ="${!i}";;
    --do-align)     ((i++)); (( i <= argc )) || { echo "Error: --do-align requires yes|no" >&2; exit 1; }; DO_ALIGN="${!i}";;
    --do-depth)     ((i++)); (( i <= argc )) || { echo "Error: --do-depth requires yes|no" >&2; exit 1; }; DO_DEPTH="${!i}";;
    --do-chimera)   ((i++)); (( i <= argc )) || { echo "Error: --do-chimera requires yes|no" >&2; exit 1; }; DO_CHIMERA="${!i}";;
    --dm6-ref)      ((i++)); (( i <= argc )) || { echo "Error: --dm6-ref requires a value" >&2; exit 1; }; DM6_REF="${!i}";;
    --chimera-jobs) ((i++)); (( i <= argc )) || { echo "Error: --chimera-jobs requires a value" >&2; exit 1; }; CHIMERA_JOBS="${!i}";;
    --sample)
      # Expect a block: --fq1-abs X --fq1-short X --fq2-abs X --fq2-short X --label X
      fq1_abs=""; fq1_short=""; fq2_abs=""; fq2_short=""; label=""
      while (( i < argc )); do
        ((i++))
        key="${!i:-}"
        case "$key" in
          --fq1-abs)   ((i++)); fq1_abs="${!i:-}";;
          --fq1-short) ((i++)); fq1_short="${!i:-}";;
          --fq2-abs)   ((i++)); fq2_abs="${!i:-}";;
          --fq2-short) ((i++)); fq2_short="${!i:-}";;
          --label)     ((i++)); label="${!i:-}";;
          --sample|--workdir|--link-fastq|--do-align|--do-depth|--do-chimera|--dm6-ref|--chimera-jobs|-h|--help)
            # Step back one token so outer loop re-reads this flag
            ((i--))
            break;;
          "")
            break;;
          *)
            echo "Unknown option in --sample block: $key" >&2; exit 1;;
        esac
      done
      # Validate this sample block
      if [[ -z "$fq1_abs" || -z "$fq1_short" || -z "$fq2_abs" || -z "$fq2_short" || -z "$label" ]]; then
        echo "Error: Incomplete --sample block (need --fq1-abs/--fq1-short/--fq2-abs/--fq2-short/--label)." >&2
        exit 1
      fi
      S_FQ1_ABS+=("$fq1_abs"); S_FQ1_SHORT+=("$fq1_short")
      S_FQ2_ABS+=("$fq2_abs"); S_FQ2_SHORT+=("$fq2_short")
      S_LABEL+=("$label")
      ;;

    *)
      echo "Unknown option: $arg" >&2
      usage
      exit 1
      ;;
  esac

  ((i++))
done

# Validate
[[ -n "$WORKDIR" ]] || { echo "Error: --workdir is required." >&2; exit 1; }
[[ "${#S_LABEL[@]}" -gt 0 ]] || { echo "Error: At least one --sample block is required." >&2; exit 1; }
[[ "$LINK_FASTQ" =~ ^(yes|no)$ ]] || { echo "Error: --link-fastq must be yes|no." >&2; exit 1; }
[[ "$DO_ALIGN" =~ ^(yes|no)$ ]] || { echo "Error: --do-align must be yes|no." >&2; exit 1; }
[[ "$DO_DEPTH" =~ ^(yes|no)$ ]] || { echo "Error: --do-depth must be yes|no." >&2; exit 1; }
[[ "$DO_CHIMERA" =~ ^(yes|no)$ ]] || { echo "Error: --do-chimera must be yes|no." >&2; exit 1; }

mkdir -p "$WORKDIR"
cd "$WORKDIR"

# ==============================
# Initialize outputs
# ==============================
STAMP="$(date +%Y%m%d)"
COVERAGE_FILE="MDA_${STAMP}.coverage.count"
DISCORDANT_FILE="MDA_${STAMP}.discordant.count"
SPLITMAP_FILE="MDA_${STAMP}.splitmapping.count"
CHIMERA_FILE="MDA_${STAMP}.Interchromosomal_Inverted_Outward_Large_Insert_Unclassified_Normal.count"
TOTALBASE_FILE="MDA_${STAMP}.total_base.count"

: > "$TOTALBASE_FILE"
if [[ "$DO_ALIGN" == "yes" ]]; then
  : > "$COVERAGE_FILE"
  : > "$DISCORDANT_FILE"
  : > "$SPLITMAP_FILE"
  [[ "$DO_CHIMERA" == "yes" ]] && : > "$CHIMERA_FILE"
fi

# ==============================
# Tool checks
# ==============================
need_tool pigz
need_tool awk
need_tool flock
# need_tool test
if [[ "$DO_ALIGN" == "yes" ]]; then
  need_tool bwa
  need_tool samtools
  need_tool python3
fi

# ==============================
# Prepare inputs and launch base counts
# ==============================
for idx in "${!S_LABEL[@]}"; do
  fq1_abs="${S_FQ1_ABS[$idx]}"; fq1_short="${S_FQ1_SHORT[$idx]}"
  fq2_abs="${S_FQ2_ABS[$idx]}"; fq2_short="${S_FQ2_SHORT[$idx]}"
  label="${S_LABEL[$idx]}"

  [[ -r "$fq1_abs" ]] || { echo "Error: Cannot read $fq1_abs" >&2; exit 1; }
  [[ -r "$fq2_abs" ]] || { echo "Error: Cannot read $fq2_abs" >&2; exit 1; }

  if [[ "$LINK_FASTQ" == "yes" ]]; then
    [[ -e "$fq1_short" ]] || ln -s "$fq1_abs" "$fq1_short"
    [[ -e "$fq2_short" ]] || ln -s "$fq2_abs" "$fq2_short"
    fq1="$fq1_short"
    fq2="$fq2_short"
  else
    fq1="$fq1_abs"
    fq2="$fq2_abs"
  fi

  # Launch base counts in background (both mates)
  total_bases "$label" "$fq1" "read_1" >> "$TOTALBASE_FILE" &
  total_bases "$label" "$fq2" "read_2" >> "$TOTALBASE_FILE" &
done
# Do not wait here; allow overlap with alignment/chimera work.
# We'll use a final global 'wait' at the end.

# ==============================
# Alignment, stats, chimera (bg)
# ==============================
if [[ "$DO_ALIGN" == "yes" ]]; then
  for idx in "${!S_LABEL[@]}"; do
    label="${S_LABEL[$idx]}"
    fq1="${S_FQ1_ABS[$idx]}"; fq2="${S_FQ2_ABS[$idx]}"
    if [[ "$LINK_FASTQ" == "yes" ]]; then
      fq1="${S_FQ1_SHORT[$idx]}"
      fq2="${S_FQ2_SHORT[$idx]}"
    fi

    # 1) Align -> keep mapped (-F 0x4) -> BAM
    bwa mem -t "$THREADS_ALIGN" "$DM6_REF" "$fq1" "$fq2" \
      | samtools view -@ "$THREADS_ALIGN" -F 0x4 -bS - \
      > "${label}.pe.F4.bam"

    # 2) Sort and index
    samtools sort -@ "$THREADS_SORT" -o "${label}.pe.F4.s.bam" "${label}.pe.F4.bam"
    samtools index -@ "$THREADS_SORT" "${label}.pe.F4.s.bam"
    rm -f "${label}.pe.F4.bam"

    # 3) Filter -F 0x904 and index
    samtools view -@ "$THREADS_SAMTOOLS" -bS -F 0x904 "${label}.pe.F4.s.bam" > "${label}.pe.F904.s.bam"
    samtools index -@ "$THREADS_SAMTOOLS" "${label}.pe.F904.s.bam"

    # 4) Depth + coverage summary (optional)
    if [[ "$DO_DEPTH" == "yes" ]]; then
      samtools depth -@ "$THREADS_DEPTH" -a "${label}.pe.F904.s.bam" > "${label}.pe.F904.s.bam.depth"
      awk -v sample="$label" '
        BEGIN{total=0; over0=0}
        {total++; if($3>0) over0++}
        END{printf("%s\t%d\t%d\n", sample, over0, total)}
      ' "${label}.pe.F904.s.bam.depth" >> "$COVERAGE_FILE"
    fi

    # 5) Discordant pairs summary
    samtools view -@ "$THREADS_SAMTOOLS" "${label}.pe.F904.s.bam" \
      | awk -v sample="$label" '
        BEGIN{all=0; disc=0}
        {
          if($7 != "=" || $9 > 1000 || $9 < -1000) { disc++ }
          all++
        }
        END{printf("%s\t%d\t%d\n", sample, disc, all)}
      ' >> "$DISCORDANT_FILE"

    # 6) Split-mapping summary (SA:Z)
    samtools view -@ "$THREADS_SAMTOOLS" "${label}.pe.F904.s.bam" \
      | awk -v sample="$label" '
        BEGIN{total=0; sa=0}
        { if($0 ~ /SA:Z:/) sa++; total++ }
        END{printf("%s\t%d\t%d\n", sample, sa, total)}
      ' >> "$SPLITMAP_FILE"

    # 7) Chimera detection in background, with concurrency and safe append
    if [[ "$DO_CHIMERA" == "yes" ]]; then
      chim_bam="${label}.pe.F904.s.bam"  # choose F904 (filtered). Use F4 if you prefer.
      wait_for_bg_slots "$CHIMERA_JOBS"
      (
        python3 /rd/caiya/wmx/get_chimera.Interchromosomal_Inverted_Outward_Large_Insert_Unclassified_Normal.pe.quick.py \
          "$chim_bam" "$label" \
        | while IFS= read -r line; do
            flock_append "$CHIMERA_FILE" "$line"
          done
      ) &
    fi
  done
fi

# ==============================
# Final synchronization
# ==============================
# Wait for all background jobs (base counts + chimera) to finish
wait

echo "Done. Outputs in: $WORKDIR"
echo "Generated:"
echo " - $TOTALBASE_FILE"
if [[ "$DO_ALIGN" == "yes" ]]; then
  echo " - $COVERAGE_FILE"
  echo " - $DISCORDANT_FILE"
  echo " - $SPLITMAP_FILE"
  [[ "$DO_CHIMERA" == "yes" ]] && echo " - $CHIMERA_FILE"
fi

