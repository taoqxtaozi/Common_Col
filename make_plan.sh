#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash make_plan.sh --tree '<rooted Newick tree>' --reference H --paths path.txt --threads 60 --common_workers 15 [options]
  bash make_plan.sh --tree-file input.tre --reference H --paths path.txt --threads 60 --common_workers 15 [options]

Required arguments:
  --reference         Reference / outgroup taxon name
  --paths             Two-column file: taxon and genome path
  --threads           Passed through to RunPipelineUseThis.sh
  --common_workers    Passed through to RunPipelineUseThis.sh

One of:
  --tree              Rooted Newick tree string
  --tree-file         File containing one rooted Newick tree

Tree requirement:
  The input tree must be rooted. Partially resolved trees and rooted star trees are supported.

Optional arguments:
  --outdir            Output directory (default: current working directory)
  --python            Python executable to use for the planner (default: python3)
  --planner           Path to plan_partially_resolved_cactus.py
                      (default: same directory as this shell script)

Guide-tree generator options (passed through to generate_random_guidetrees_2models_2modes_finalver.py):
  --generator         Path to the generator script
                      (default: same directory as the planner script)
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

Consensus pipeline options (passed through to RunPipelineUseThis.sh):
  --pipeline          Path to RunPipelineUseThis.sh
                      (default: same directory as this shell script)
  --ModelFile         Path to tryMLstartree.mod
                      (default: same directory as this shell script)
  --sep_length        Pass through --sep_length
  --chrom_length_threshold
                      Pass through --chrom_length_threshold
  --global_num        Pass through --global_num
  --separate_workers  Pass through --separate_workers


Examples:
  bash make_plan2.sh --tree '(A,B,C,D)i1;' --reference A --paths path.txt --threads 60 --common_workers 15
  bash make_plan2.sh --tree-file input.tre --reference H --paths path.txt --threads 60 --common_workers 15 --outdir .
  bash make_plan2.sh --tree-file input.tre --reference H --paths path.txt --threads 60 --common_workers 15 \
    --guide-model yule --guide-rf-threshold 1 --guide-t-threshold 2/3
USAGE
}

TREE_STRING=""
TREE_FILE=""
REFERENCE=""
PATHS_FILE=""
OUTDIR="."
PYTHON_BIN="python3"
PLANNER=""
GENERATOR=""
GENERATOR_PYTHON=""
GUIDE_NUM_TREES=""
GUIDE_MODEL=""
GUIDE_RF_THRESHOLD=""
GUIDE_T_THRESHOLD=""
GUIDE_MAX_TRIES=""
GUIDE_RESTARTS=""
GUIDE_SEED=""
PIPELINE=""
PIPELINE_THREADS=""
PIPELINE_COMMON_WORKERS=""
MODEL_FILE=""
PIPELINE_SEP_LENGTH=""
PIPELINE_CHROM_LENGTH_THRESHOLD=""
PIPELINE_GLOBAL_NUM=""
PIPELINE_SEPARATE_WORKERS=""

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
    --reference)
      REFERENCE="${2:-}"
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
    --planner)
      PLANNER="${2:-}"
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
    --pipeline)
      PIPELINE="${2:-}"
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

if [[ -n "$TREE_STRING" && -n "$TREE_FILE" ]]; then
  echo "Error: use only one of --tree or --tree-file." >&2
  exit 1
fi

if [[ -z "$TREE_STRING" && -z "$TREE_FILE" ]]; then
  echo "Error: you must provide either --tree or --tree-file." >&2
  exit 1
fi

if [[ -z "$REFERENCE" ]]; then
  echo "Error: --reference is required." >&2
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

if [[ -n "$TREE_FILE" ]]; then
  if [[ ! -f "$TREE_FILE" ]]; then
    echo "Error: tree file not found: $TREE_FILE" >&2
    exit 1
  fi
  TREE_STRING="$(tr -d '\r' < "$TREE_FILE" | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g; s/^ //; s/ $//')"
fi

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "Error: Python executable not found: $PYTHON_BIN" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -z "$PLANNER" ]]; then
  PLANNER="$SCRIPT_DIR/plan_partially_resolved_cactus.py"
fi

if [[ ! -f "$PLANNER" ]]; then
  echo "Error: planner script not found: $PLANNER" >&2
  exit 1
fi

if [[ -z "$GENERATOR" ]]; then
  GENERATOR="$(cd "$(dirname "$PLANNER")" && pwd)/generate_random_guidetrees_2models_2modes_finalver.py"
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

if [[ -z "$PIPELINE" ]]; then
  PIPELINE="$SCRIPT_DIR/RunPipelineUseThis.sh"
fi

if [[ ! -f "$PIPELINE" ]]; then
  echo "Error: consensus pipeline script not found: $PIPELINE" >&2
  exit 1
fi

if [[ -z "$MODEL_FILE" ]]; then
  MODEL_FILE="$SCRIPT_DIR/tryMLstartree.mod"
fi

if [[ ! -f "$MODEL_FILE" ]]; then
  echo "Error: model file not found: $MODEL_FILE" >&2
  exit 1
fi

mkdir -p "$OUTDIR"

CMD=( "$PYTHON_BIN" "$PLANNER"
  --tree "$TREE_STRING"
  --reference "$REFERENCE"
  --paths "$PATHS_FILE"
  --outdir "$OUTDIR"
  --generator "$GENERATOR"
  --generator-python "$GENERATOR_PYTHON"
  --pipeline "$PIPELINE"
  --threads "$PIPELINE_THREADS"
  --common_workers "$PIPELINE_COMMON_WORKERS"
  --ModelFile "$MODEL_FILE"
)

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

exec "${CMD[@]}"
