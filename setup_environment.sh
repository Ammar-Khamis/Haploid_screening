#!/usr/bin/env bash
set -euo pipefail

# Install micromamba into ~/.local/bin if it is not already present
MICROMAMBA="$HOME/.local/bin/micromamba"
if [ ! -x "$MICROMAMBA" ]; then
  mkdir -p "$HOME/.local/bin"
  curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xvj --strip-components=1 bin/micromamba -C /tmp
  install -m755 /tmp/micromamba "$MICROMAMBA"
fi

# Enable micromamba shell support
source <($MICROMAMBA shell hook -s bash)

# Link the conda command to micromamba for compatibility
ln -sf "$MICROMAMBA" "$HOME/.local/bin/conda"
if ! grep -q 'alias conda=micromamba' ~/.bashrc 2>/dev/null; then
  echo 'alias conda=micromamba' >> ~/.bashrc
fi

# Configure default channels
micromamba config set always_yes yes
micromamba config append channels conda-forge
micromamba config append channels bioconda

# Create main Python environment
micromamba create -n haploid python=3.12 -y
micromamba activate haploid
python -m pip install --upgrade pip
pip install -r requirements.txt --no-cache-dir

# Additional environments required by the workflow
micromamba create -n biotools samtools bedtools minimap2 -y
micromamba create -n last_env last perl perl-bioperl -y

# Convenience aliases
if ! grep -q 'use_biotools' ~/.bashrc 2>/dev/null; then
  echo 'alias use_biotools="conda activate biotools"' >> ~/.bashrc
fi
if ! grep -q 'use_last' ~/.bashrc 2>/dev/null; then
  echo 'alias use_last="conda activate last_env"' >> ~/.bashrc
fi
