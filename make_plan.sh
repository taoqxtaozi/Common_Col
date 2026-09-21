#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:

  1) Default mode: fixed global reference
     bash make_plan.sh \
       --tree '<rooted Newick tree>' \
       --reference H \
       --paths path.txt \
       --threads 60 \
       --common_workers 15 \
       [options]

     bash make_plan.sh \
       --tree-file input.tre \
       --reference H \
       --paths path.txt \
       --threads 60 \
       --common_workers 15 \
       [options]

  2) No-fixed-reference mode
     bash make_plan.sh \
       --tree '<rooted Newick tree>' \
       --noFixedRef 1 \
       --final_reference H \
       --paths path.txt \
       --threads 60 \
       --common_workers 15 \
       [options]

     bash make_plan.sh \
       --tree-file input.tre \
       --noFixedRef 1 \
       --final_reference H \
       --paths path.txt \
       --threads 60 \
       --common_workers 15 \
       [options]

Required in all modes:
  --paths             Two-column file: taxon and genome path
  --threads           Passed through to the selected planner
  --common_workers    Passed through to the selected planner

One of:
  --tree              Rooted Newick tree string
  --tree-file         File containing one rooted Newick tree

Reference mode:
  --noFixedRef        0 or 1. Default: 0.

  --noFixedRef 0:
      Uses plan_GlobalRef.py.
      A global reference is required.
      Specify it with:
        --reference TAXON
      or equivalently:
        --global_reference TAXON

  --noFixedRef 1:
      Uses plan_noFixedRef.py automatically.
      No global reference is used for planning, guide-tree generation,
      Cactus alignment, or per-polytomy consensus extraction.
      Each consensus task selects the longest directly participating genome
      as its local reference.
      --final_reference is required only for the final complete HAL -> MAF export.

Reference arguments:
  --reference         Global reference used when --noFixedRef 0.
                      Alias: --global_reference
  --global_reference  Alias of --reference.
  --final_reference   Final HAL -> MAF export reference used only when
                      --noFixedRef 1.

Tree requirement:
  The input tree must be rooted.
  Partially resolved trees and rooted star trees are supported.

Optional arguments:
  --outdir            Output directory (default: current working directory)
  --python            Python executable used to run the planner (default: python3)

Guide-tree generator options:
  --generator         Path to generate_random_guidetrees_2models_2modes_finalver.py
                      (default: same directory as the selected planner)
  --generator-python  Python executable used to run the generator
                      (default: same as --python)
  --guide-num-trees   Pass through --num_trees
  --guide-model       Pass through --model
  --guide-rf-threshold
                      Pass through --rf-threshold
  --guide-t-threshold
                      Pass through --t-threshold
  --guide-max-tries   Pass through --max-tries
  --guide-restarts    Pass through --restarts
  --guide-seed        Pass through --seed

Consensus pipeline options:
  --ModelFile         Path to tryMLstartree.mod
                      (default: same directory as this shell script)
  --sep_length        Pass through --sep_length
  --chrom_length_threshold
                      Pass through --chrom_length_threshold
  --global_num        Pass through --global_num
  --separate_workers  Pass through --separate_workers

No-fixed-reference-only option:
  --reference-selector
                      Path to select_longest_reference.py.
                      Used only when --noFixedRef 1.
                      Default: same directory as plan_noFixedRef.py.

Notes:
  * The selected planner is automatic:
      --noFixedRef 0 -> plan_GlobalRef.py
      --noFixedRef 1 -> plan_noFixedRef.py
  * Both planners are expected to use RunPipelineUseThis.sh as their own
    default consensus pipeline. make_plan.sh does not pass --pipeline.

Examples:

  # Default: one fixed global reference
  bash make_plan.sh \
    --tree '(A,B,C,D)i1;' \
    --reference A \
    --paths path.txt \
    --threads 60 \
    --common_workers 15

  # Same default mode, using the explicit alias --global_reference
  bash make_plan.sh \
    --tree-file input.tre \
    --global_reference H \
    --paths path.txt \
    --threads 60 \
    --common_workers 15 \
    --outdir .

  # No fixed reference:
  # local reference of each consensus task = longest participating genome
  # final_reference is only used for final HAL -> MAF export
  bash make_plan.sh \
    --tree-file input.tre \
    --noFixedRef 1 \
    --final_reference H \
    --paths path.txt \
    --threads 60 \
    --common_workers 15 \
    --outdir .

  # No fixed reference + guide-tree options
  bash make_plan.sh \
    --tree-file input.tre \
    --noFixedRef 1 \
    --final_reference H \
    --paths path.txt \
    --threads 60 \
    --common_workers 15 \
    --guide-model yule \
    --guide-rf-threshold 1 \
    --guide-t-threshold 2/3
USAGE
}


########################################
# 0. Default parameter values
########################################

TREE_STRING=""
TREE_FILE=""

NO_FIXED_REF=0
GLOBAL_REFERENCE=""
FINAL_REFERENCE=""

PATHS_FILE=""
OUTDIR="."
PYTHON_BIN="python3"

GENERATOR=""
GENERATOR_PYTHON=""

GUIDE_NUM_TREES=""
GUIDE_MODEL=""
GUIDE_RF_THRESHOLD=""
GUIDE_T_THRESHOLD=""
GUIDE_MAX_TRIES=""
GUIDE_RESTARTS=""
GUIDE_SEED=""

REFERENCE_SELECTOR=""
PIPELINE_THREADS=""
PIPELINE_COMMON_WORKERS=""
MODEL_FILE=""
PIPELINE_SEP_LENGTH=""
PIPELINE_CHROM_LENGTH_THRESHOLD=""
PIPELINE_GLOBAL_NUM=""
PIPELINE_SEPARATE_WORKERS=""


########################################
# 1. Parse command-line arguments
########################################

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tree)
      TREE_STRING="${2:-}"
      shift 2
      ;;

    --tree-file)
      TREE_FILE="${2:-}"
      shift 2
      ;;

    --noFixedRef)
      NO_FIXED_REF="${2:-}"
      shift 2
      ;;

    --reference|--global_reference)
      GLOBAL_REFERENCE="${2:-}"
      shift 2
      ;;

    --final_reference)
      FINAL_REFERENCE="${2:-}"
      shift 2
      ;;

    --paths)
      PATHS_FILE="${2:-}"
      shift 2
      ;;

    --outdir)
      OUTDIR="${2:-}"
      shift 2
      ;;

    --python)
      PYTHON_BIN="${2:-}"
      shift 2
      ;;

    --generator)
      GENERATOR="${2:-}"
      shift 2
      ;;

    --generator-python)
      GENERATOR_PYTHON="${2:-}"
      shift 2
      ;;

    --guide-num-trees)
      GUIDE_NUM_TREES="${2:-}"
      shift 2
      ;;

    --guide-model)
      GUIDE_MODEL="${2:-}"
      shift 2
      ;;

    --guide-rf-threshold)
      GUIDE_RF_THRESHOLD="${2:-}"
      shift 2
      ;;

    --guide-t-threshold)
      GUIDE_T_THRESHOLD="${2:-}"
      shift 2
      ;;

    --guide-max-tries)
      GUIDE_MAX_TRIES="${2:-}"
      shift 2
      ;;

    --guide-restarts)
      GUIDE_RESTARTS="${2:-}"
      shift 2
      ;;

    --guide-seed)
      GUIDE_SEED="${2:-}"
      shift 2
      ;;

    --reference-selector)
      REFERENCE_SELECTOR="${2:-}"
      shift 2
      ;;

    --threads)
      PIPELINE_THREADS="${2:-}"
      shift 2
      ;;

    --common_workers)
      PIPELINE_COMMON_WORKERS="${2:-}"
      shift 2
      ;;

    --ModelFile)
      MODEL_FILE="${2:-}"
      shift 2
      ;;

    --sep_length)
      PIPELINE_SEP_LENGTH="${2:-}"
      shift 2
      ;;

    --chrom_length_threshold)
      PIPELINE_CHROM_LENGTH_THRESHOLD="${2:-}"
      shift 2
      ;;

    --global_num)
      PIPELINE_GLOBAL_NUM="${2:-}"
      shift 2
      ;;

    --separate_workers)
      PIPELINE_SEPARATE_WORKERS="${2:-}"
      shift 2
      ;;

    -h|--help)
      usage
      exit 0
      ;;

    *)
      echo "Error: unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done


########################################
# 2. Check common arguments
########################################

if [[ -n "$TREE_STRING" && -n "$TREE_FILE" ]]; then
  echo "Error: use only one of --tree or --tree-file." >&2
  exit 1
fi

if [[ -z "$TREE_STRING" && -z "$TREE_FILE" ]]; then
  echo "Error: you must provide either --tree or --tree-file." >&2
  exit 1
fi

if [[ "$NO_FIXED_REF" != "0" && "$NO_FIXED_REF" != "1" ]]; then
  echo "Error: --noFixedRef must be 0 or 1." >&2
  exit 1
fi

if [[ -z "$PATHS_FILE" ]]; then
  echo "Error: --paths is required." >&2
  exit 1
fi

if [[ -z "$PIPELINE_THREADS" ]]; then
  echo "Error: --threads is required." >&2
  exit 1
fi

if [[ -z "$PIPELINE_COMMON_WORKERS" ]]; then
  echo "Error: --common_workers is required." >&2
  exit 1
fi

if [[ ! -f "$PATHS_FILE" ]]; then
  echo "Error: path file not found: $PATHS_FILE" >&2
  exit 1
fi


########################################
# 3. Check reference arguments by mode
########################################

if [[ "$NO_FIXED_REF" == "0" ]]; then
  if [[ -z "$GLOBAL_REFERENCE" ]]; then
    echo "Error: --reference (or --global_reference) is required when --noFixedRef 0." >&2
    exit 1
  fi

  if [[ -n "$FINAL_REFERENCE" ]]; then
    echo "Error: --final_reference is only used when --noFixedRef 1." >&2
    echo "       In the default fixed-reference mode, use --reference instead." >&2
    exit 1
  fi
else
  if [[ -n "$GLOBAL_REFERENCE" ]]; then
    echo "Error: --reference/--global_reference must not be used when --noFixedRef 1." >&2
    echo "       This mode has no global reference." >&2
    exit 1
  fi

  if [[ -z "$FINAL_REFERENCE" ]]; then
    echo "Error: --final_reference is required when --noFixedRef 1." >&2
    echo "       It is used only for the final complete HAL -> MAF export." >&2
    exit 1
  fi
fi


########################################
# 4. Read tree file if supplied
########################################

if [[ -n "$TREE_FILE" ]]; then
  if [[ ! -f "$TREE_FILE" ]]; then
    echo "Error: tree file not found: $TREE_FILE" >&2
    exit 1
  fi

  TREE_STRING="$(
    tr -d '\r' < "$TREE_FILE" \
      | tr '\n' ' ' \
      | sed 's/[[:space:]]\+/ /g; s/^ //; s/ $//'
  )"
fi


########################################
# 5. Resolve program paths
########################################

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "Error: Python executable not found: $PYTHON_BIN" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


# Planner is selected automatically from --noFixedRef.
if [[ "$NO_FIXED_REF" == "1" ]]; then
  PLANNER="$SCRIPT_DIR/plan_noFixedRef.py"
else
  PLANNER="$SCRIPT_DIR/plan_GlobalRef.py"
fi

if [[ ! -f "$PLANNER" ]]; then
  echo "Error: planner script not found: $PLANNER" >&2

  if [[ "$NO_FIXED_REF" == "1" ]]; then
    echo "Expected plan_noFixedRef.py for --noFixedRef 1." >&2
  else
    echo "Expected plan_GlobalRef.py for --noFixedRef 0." >&2
  fi

  exit 1
fi

PLANNER_DIR="$(cd "$(dirname "$PLANNER")" && pwd)"


if [[ -z "$GENERATOR" ]]; then
  GENERATOR="$PLANNER_DIR/generate_random_guidetrees_2models_2modes_finalver.py"
fi

if [[ ! -f "$GENERATOR" ]]; then
  echo "Error: guide-tree generator script not found: $GENERATOR" >&2
  exit 1
fi


if [[ -z "$GENERATOR_PYTHON" ]]; then
  GENERATOR_PYTHON="$PYTHON_BIN"
fi

if ! command -v "$GENERATOR_PYTHON" >/dev/null 2>&1; then
  echo "Error: generator Python executable not found: $GENERATOR_PYTHON" >&2
  exit 1
fi


# Only the no-fixed-reference planner needs the local-reference selector.
if [[ "$NO_FIXED_REF" == "1" ]]; then
  if [[ -z "$REFERENCE_SELECTOR" ]]; then
    REFERENCE_SELECTOR="$PLANNER_DIR/select_longest_reference.py"
  fi

  if [[ ! -f "$REFERENCE_SELECTOR" ]]; then
    echo "Error: local-reference selector script not found: $REFERENCE_SELECTOR" >&2
    exit 1
  fi
else
  if [[ -n "$REFERENCE_SELECTOR" ]]; then
    echo "Error: --reference-selector is only valid when --noFixedRef 1." >&2
    exit 1
  fi
fi


if [[ -z "$MODEL_FILE" ]]; then
  MODEL_FILE="$SCRIPT_DIR/tryMLstartree.mod"
fi

if [[ ! -f "$MODEL_FILE" ]]; then
  echo "Error: model file not found: $MODEL_FILE" >&2
  exit 1
fi


########################################
# 6. Prepare output directory
########################################

mkdir -p "$OUTDIR"


########################################
# 7. Build planner command
########################################

CMD=(
  "$PYTHON_BIN" "$PLANNER"
  --tree "$TREE_STRING"
  --paths "$PATHS_FILE"
  --outdir "$OUTDIR"
  --generator "$GENERATOR"
  --generator-python "$GENERATOR_PYTHON"
  --threads "$PIPELINE_THREADS"
  --common_workers "$PIPELINE_COMMON_WORKERS"
  --ModelFile "$MODEL_FILE"
)


# Mode-specific planner arguments.
if [[ "$NO_FIXED_REF" == "1" ]]; then
  CMD+=(
    --final_reference "$FINAL_REFERENCE"
    --reference-selector "$REFERENCE_SELECTOR"
  )
else
  CMD+=(
    --reference "$GLOBAL_REFERENCE"
  )
fi


########################################
# 8. Optional guide-tree arguments
########################################

if [[ -n "$GUIDE_NUM_TREES" ]]; then
  CMD+=( --guide-num-trees "$GUIDE_NUM_TREES" )
fi

if [[ -n "$GUIDE_MODEL" ]]; then
  CMD+=( --guide-model "$GUIDE_MODEL" )
fi

if [[ -n "$GUIDE_RF_THRESHOLD" ]]; then
  CMD+=( --guide-rf-threshold "$GUIDE_RF_THRESHOLD" )
fi

if [[ -n "$GUIDE_T_THRESHOLD" ]]; then
  CMD+=( --guide-t-threshold "$GUIDE_T_THRESHOLD" )
fi

if [[ -n "$GUIDE_MAX_TRIES" ]]; then
  CMD+=( --guide-max-tries "$GUIDE_MAX_TRIES" )
fi

if [[ -n "$GUIDE_RESTARTS" ]]; then
  CMD+=( --guide-restarts "$GUIDE_RESTARTS" )
fi

if [[ -n "$GUIDE_SEED" ]]; then
  CMD+=( --guide-seed "$GUIDE_SEED" )
fi


########################################
# 9. Optional consensus-pipeline arguments
########################################

if [[ -n "$PIPELINE_SEP_LENGTH" ]]; then
  CMD+=( --sep_length "$PIPELINE_SEP_LENGTH" )
fi

if [[ -n "$PIPELINE_CHROM_LENGTH_THRESHOLD" ]]; then
  CMD+=( --chrom_length_threshold "$PIPELINE_CHROM_LENGTH_THRESHOLD" )
fi

if [[ -n "$PIPELINE_GLOBAL_NUM" ]]; then
  CMD+=( --global_num "$PIPELINE_GLOBAL_NUM" )
fi

if [[ -n "$PIPELINE_SEPARATE_WORKERS" ]]; then
  CMD+=( --separate_workers "$PIPELINE_SEPARATE_WORKERS" )
fi


########################################
# 10. Run selected planner
########################################

if [[ "$NO_FIXED_REF" == "1" ]]; then
  echo "Planning mode: no fixed global reference"
  echo "Planner: $PLANNER"
  echo "Final HAL-to-MAF reference: $FINAL_REFERENCE"
else
  echo "Planning mode: fixed global reference"
  echo "Planner: $PLANNER"
  echo "Global reference: $GLOBAL_REFERENCE"
fi

exec "${CMD[@]}"
