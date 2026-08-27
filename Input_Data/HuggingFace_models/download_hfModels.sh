#!/bin/bash
set -e

# Optional: helps avoid download stalls on some HPC/shared filesystems
export HF_HUB_DISABLE_XET=1

echo "Downloading ChemBERTa zinc base..."
hf download seyonec/ChemBERTa-zinc-base-v1 \
  --local-dir ./ChemBERTa_zinc_base

echo "Downloading ChemBERTa 77M MLM..."
hf download DeepChem/ChemBERTa-77M-MLM \
  --local-dir ./ChemBERTa_77M_MLM

echo "Downloading ChemBERTa 77M MTR..."
hf download DeepChem/ChemBERTa-77M-MTR \
  --local-dir ./ChemBERTa_77M_MTR

echo "Downloading ChemBERTa 10M MLM..."
hf download DeepChem/ChemBERTa-10M-MLM \
  --local-dir ./ChemBERTa_10M_MLM

echo "Downloading ChemBERTa 100M MLM..."
hf download DeepChem/ChemBERTa-100M-MLM \
  --local-dir ./ChemBERTa_100M_MLM

echo "Downloading MoLFormer XL..."
hf download ibm-research/MoLFormer-XL-both-10pct \
  --local-dir ./MoLFormer_XL_both_10pct

echo "Downloading MIST-28M base model..."
hf download mist-models/mist-28M-ti624ev1 \
  --local-dir ./MIST_28M

echo "All models downloaded successfully."
