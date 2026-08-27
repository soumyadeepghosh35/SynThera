---
license: mit
library_name: transformers
pipeline_tag: fill-mask
tags:
- cheminformatics
- ChemBERTa
- masked-lm
- roberta
---

# ChemBERTa-100M-MLM

ChemBERTa model pretrained on a subset of 100M molecules from ZINC20 dataset using masked language modeling (MLM).

## Usage

```python
from transformers import AutoTokenizer, AutoModelForMaskedLM

tokenizer = AutoTokenizer.from_pretrained("DeepChem/ChemBERTa-100M-MLM")
model = AutoModelForMaskedLM.from_pretrained("DeepChem/ChemBERTa-100M-MLM")


