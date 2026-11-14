#!/usr/bin/env bash
# ==============================================================================
# Phase 5 Google Colab / Cloud GPU Execution Script
# ==============================================================================
# This script sets up OpenMM on a cloud GPU instance (e.g. Google Colab T4/V100/A100)
# and executes the 6 independent 20 ns production MD replicates (120 ns total).
#
# Usage on Colab:
#   bash scripts/run_colab_production.sh [wt|cua|all] [replicate_num (1-3)]
# ==============================================================================

set -e

echo "================================================================================"
echo "Phase 5 Production MD: Cloud GPU Setup & Runner"
echo "================================================================================"

# 1. Verify GPU
nvidia-smi

# 2. Install dependencies
pip install -q openmm pyyaml mdtraj numpy

# 3. Detect Platform
python -c "
import openmm as mm
platforms = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
print('Available OpenMM platforms:', platforms)
"

TARGET_CONSTRUCT=${1:-all}
TARGET_REP=${2:-all}

echo "Target construct: $TARGET_CONSTRUCT | Target replicate: $TARGET_REP"

if [ "$TARGET_REP" = "all" ]; then
    python scripts/07_run_production.py --construct "$TARGET_CONSTRUCT" --duration-ns 20.0 --platform CUDA
else
    python scripts/07_run_production.py --construct "$TARGET_CONSTRUCT" --replicate "$TARGET_REP" --duration-ns 20.0 --platform CUDA
fi

echo "================================================================================"
echo "Production execution finished. Running QC validation..."
echo "================================================================================"
python scripts/07_run_production.py --validate-qc --duration-ns 20.0
