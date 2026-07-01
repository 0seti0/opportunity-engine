#!/bin/bash
# Covalent Radar — bi-weekly run (launched by com.axiom.covalent-radar on the 1st & 15th).
# Runs the deterministic legs, unifies, and posts the fact-checked digest to Slack.
# Continue-on-error (set +e): a slow/failed radar step must not block the digest, which runs on existing data.
export PATH="$HOME/.local/bin:/opt/homebrew/bin:/usr/local/bin:/usr/bin:/bin"   # launchd PATH is minimal; add claude/uv/java
export HF_HUB_OFFLINE=1                                                          # MolNexTR checkpoint is local — no web revalidation

cd "/Users/seti/Desktop/sov ai/opportunity-engine/layer1" || exit 1
PY="./.venv-engine/bin/python"
mkdir -p runs
LOG="runs/biweekly_$(date +%Y%m%d_%H%M).log"

{
  echo "=== Covalent Radar bi-weekly run: $(date) ==="
  echo "--- census_topup: Day-1 covalent patents on the watchlist (incl. RAF1) ---"; "$PY" census_topup.py
  echo "--- clinical_radar: fresh CT.gov + EU CTIS trials ---";                       "$PY" clinical_radar.py
  echo "--- unify: rebuild the covalent-move log ---";                                "$PY" unify.py
  echo "--- notify: compose -> web fact-check -> post to Slack ---";                  "$PY" notify.py
  echo "=== done: $(date) ==="
} >> "$LOG" 2>&1
