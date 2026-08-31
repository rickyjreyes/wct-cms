#!/usr/bin/env bash
set -euo pipefail

INPUT_MANIFEST="data/files_run2016g.txt"
GOLDEN_JSON="data/cert/14221/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt"
OUTPUT_DIR="results/run2016g_frozen"

if [[ ! -f "$INPUT_MANIFEST" ]]; then
  echo "Missing input manifest: $INPUT_MANIFEST" >&2
  exit 1
fi

if [[ ! -f "$GOLDEN_JSON" ]]; then
  echo "Missing golden JSON: $GOLDEN_JSON" >&2
  exit 1
fi

ROOT_FILE=$(grep -v '^[[:space:]]*#' "$INPUT_MANIFEST" | sed '/^[[:space:]]*$/d' | head -n 1)
if [[ -z "$ROOT_FILE" ]]; then
  echo "Input manifest contains no ROOT file" >&2
  exit 1
fi

if [[ "$ROOT_FILE" != root://* && ! -f "$ROOT_FILE" ]]; then
  echo "Run2016G ROOT file is not present: $ROOT_FILE" >&2
  echo "Expected first-file sample from CERN Open Data record 30522." >&2
  exit 1
fi

cms-wct \
  --input "$INPUT_MANIFEST" \
  --golden-json "$GOLDEN_JSON" \
  --output-dir "$OUTPUT_DIR" \
  --mass-min 2 \
  --mass-max 120 \
  --bins 350 \
  --log-bins \
  --muon-pt-min 4 \
  --muon-eta-max 2.4 \
  --tight-id \
  --max-events 100000 \
  --fit-degree 7 \
  --mask-window 2.9:3.3 \
  --mask-window 3.55:3.85 \
  --mask-window 8.5:11.5 \
  --mask-window 80:100 \
  --omega-min 3.1 \
  --omega-max 80 \
  --omega-steps 1000 \
  --frozen-omega 7.025825825825827 \
  --permutations 1000 \
  --parametric-bootstrap 500 \
  --seed 20260827

echo
echo "Frozen Run2016G result written to: $OUTPUT_DIR/summary.json"
