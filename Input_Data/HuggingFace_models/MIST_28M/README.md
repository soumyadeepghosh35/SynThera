---
language: en
library_name: transformers
license: apache-2.0
tags:
- mist
- chemistry
- molecular-property-prediction
---

# MIST: Molecular Insight SMILES Transformers


MIST is a family of molecular foundation models for molecular property prediction.
The models were pre-trained on SMILES strings from the [Enamine REAL Space](https://enamine.net/compound-collections/real-compounds/real-space-navigator) dataset using the Masked Language Modeling (MLM) objective, then fine-tuned for downstream prediction tasks.
Further information is available in [our pre-print on arVix](https://arxiv.org/abs/2510.18900).

## Model Details

### Model Description

This is a pre-trained MIST encoder with 28M parameters trained on ~14B molecular SMILES tokens.

- **Developed by:** [Electrochemical Energy Group](https://eeg.engin.umich.edu/), University of Michigan, Ann Arbor.
- **Model type:** Self-supervised pre-trained MIST encoder with supervised finetuning.

### Model Sources

- **Repository:** [Full MIST Code](https://github.com/BattModels/mist)
- **Paper:** [arVix Preprint](https://arxiv.org/abs/2510.18900)
- **Demo:** [Finetuning and Inference Demo](https://github.com/BattModels/mist-demo)


## Uses

<!-- Address questions around how the model is intended to be used, including the foreseeable users of the model and those affected by the model. -->


### How to Get Started with the Model

Use the code below to get started with the model.

#### Setting Up Your Environment

Create a virtual environment and install dependencies:

```bash
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

> **Note**: SMIRK tokenizers require Rust to be installed. See the [Rust installation guide](https://www.rust-lang.org/tools/install) for details.

or install [uv](https://docs.astral.sh/uv/getting-started/installation/#installation-methods) and the run scripts below with ``uv run ...``.

#### Zero-Shot Similarity
Evaluating molecular similarity using MIST embeddings.
<details>
  <summary>Demo Code</summary>

  ```python
    """
    Minimal example for getting embeddings from a pretrained model and calculating similarity.

    Usage:
        uv run embeddings_similarity.py

    """
    # /// script
    # requires-python = ">=3.9"
    # dependencies = [
    #     "torch>=2.0.0",
    #     "transformers>=4.30.0",
    #     "smirk>=0.1.0",
    #     "numpy>=1.20.0",
    #     "rdkit>=2022.0.0",
    # ]
    # ///
    
  
    import torch
    import numpy as np
    from rdkit import Chem
    from smirk import SmirkTokenizerFast
    from transformers import AutoModel, AutoTokenizer
  
  
    def kekulize_smiles(smiles):
        """Convert SMILES to kekulized form."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        Chem.Kekulize(mol)
        return Chem.MolToSmiles(mol, kekuleSmiles=True)
  

    def get_embeddings(smiles_list, model, tokenizer, device="cpu"):
        """Get embeddings for a list of SMILES strings."""
        # MIST was pretrained on Kekulize SMILES
        kekulized_smiles = [kekulize_smiles(s) for s in smiles_list]

        # Tokenize
        inputs = tokenizer(
            kekulized_smiles,
            padding=True,
            truncation=True,
            max_length=512,
            return_tensors="pt"
        )
  
        # Move to device
        inputs = {k: v.to(device) for k, v in inputs.items()}
  
        # Get embeddings
        model.eval()
        with torch.no_grad():
            outputs = model(**inputs)
            embeddings = outputs.last_hidden_state[:, 0, :].cpu().numpy()

        return embeddings


    def cosine_similarity(emb1, emb2):
        """Calculate cosine similarity between two embeddings."""
        dot_product = np.dot(emb1, emb2)
        norm1 = np.linalg.norm(emb1)
        norm2 = np.linalg.norm(emb2)
        return dot_product / (norm1 * norm2)


    def main():
        # Load pretrained model
        model_path = "mist-models/mist-28M-ti624ev1"
        print(f"Loading model from {model_path}...")
        tokenizer = AutoTokenizer.from_pretrained(model_path, trust_remote_code=True)
        model = AutoModel.from_pretrained(model_path, trust_remote_code=True)
  
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model = model.to(device)
        print(f"Using device: {device}")
  
        # Example SMILES
        smiles_list = [
            "CCO",           # Ethanol
            "CC(C)O",        # Isopropanol
            "CCCO",          # Propanol
            "c1ccccc1",      # Benzene
            "c1ccccc1O",     # Phenol
            "CC(=O)O",       # Acetic acid
        ]
  
        print(f"\nGetting embeddings for {len(smiles_list)} molecules...")
        embeddings = get_embeddings(smiles_list, model, tokenizer, device)
        print(f"Embedding shape: {embeddings.shape}")
  
       
        # Query similarity example
        query = "CCCCO"  # Butanol
        print(f"Query: {query}")
  
        query_emb = get_embeddings([query], model, tokenizer, device)[0]
  
        print("\nMost similar molecules:")
        query_sims = []
        for i, smiles in enumerate(smiles_list):
            sim = cosine_similarity(query_emb, embeddings[i])
            query_sims.append((smiles, sim))
  
        query_sims.sort(key=lambda x: x[1], reverse=True)
        for i, (smiles, sim) in enumerate(query_sims[:3], 1):
            print(f"  {i}. {smiles:15s} : {sim:.4f}")

    if __name__ == "__main__":
          main()
  ```
</details>

### Fine-tuning for Property Prediction
Fine-tune a pre-trained MIST encoder for a dataset/ property of your interest.
<details>
  <summary>Demo Code</summary>

  ```python
    """
    Minimal example for finetuning a pretrained model on a CSV using HuggingFace Trainer.

    CSV format:
        smiles,target
        CCO,1.23
        CC(C)O,2.45

    Usage:
        uv run finetune_minimal.py
    """

    # /// script
    # requires-python = ">=3.9"
    # dependencies = [
    #     "torch>=2.0.0",
    #     "transformers>=4.30.0",
    #     "datasets>=2.0.0",
    #     "smirk>=0.1.0",
    #     "accelerate>=0.26.0",
    #     "rdkit>=2022.0.0",
    # ]
    # ///

    import torch
    import torch.nn as nn
    from rdkit import Chem
    from smirk import SmirkTokenizerFast
    from datasets import load_dataset
    from transformers import (
        AutoModel,
        AutoTokenizer,
        Trainer,
        TrainingArguments,
        DataCollatorWithPadding,
    )
    from pathlib import Path


    def kekulize_smiles(smiles):
        """Convert SMILES to kekulized form."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        Chem.Kekulize(mol)
        return Chem.MolToSmiles(mol, kekuleSmiles=True)


    class RegressionModel(nn.Module):
        """Model with encoder + regression task head."""

        def __init__(self, encoder, hidden_size=768, dropout=0.1):
            super().__init__()
            self.encoder = encoder
            self.task_head = nn.Sequential(
                nn.Linear(hidden_size, hidden_size),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size, 1),
            )

        def forward(self, input_ids, attention_mask, labels=None):
            # Get encoder outputs
            encoder_output = self.encoder(
                input_ids=input_ids, attention_mask=attention_mask
            )
            # Use first token
            pooled = encoder_output.last_hidden_state[:, 0, :]

            # Regression prediction
            logits = self.task_head(pooled)

            loss = None
            if labels is not None:
                loss_fn = nn.MSELoss()
                loss = loss_fn(logits.squeeze(-1), labels)

            return {"loss": loss, "logits": logits} if loss is not None else {"logits": logits}


    def tokenize_function(examples, tokenizer):
        """Tokenize SMILES strings (kekulized)."""
        # MIST was pretrained on kekulized SMILES
        kekulized = [kekulize_smiles(s) for s in examples["smiles"]]
        return tokenizer(
            kekulized,
            padding="max_length",
            truncation=True,
            max_length=512,
        )


    def main():
        # 1. Load dataset from CSV
        dataset = load_dataset(
            "csv",
            data_files={"train": "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/qm9.csv"},
        )["train"]

        # Split into train/test
        dataset = dataset.train_test_split(test_size=0.2, seed=42)

        # 2. Load pretrained encoder and tokenizer
        model_path = "mist-models/mist-28M-ti624ev1"
        tokenizer = AutoTokenizer.from_pretrained(model_path, trust_remote_code=True)
        encoder = AutoModel.from_pretrained(model_path, trust_remote_code=True)

        # 3. Create regression model with task head
        model = RegressionModel(
            encoder=encoder,
            hidden_size=encoder.config.hidden_size,
            dropout=0.1,
        )

        # 4. Tokenize dataset
        dataset = dataset.map(
            lambda x: tokenize_function(x, tokenizer),
            batched=True,
            desc="Tokenizing",
        )

        # Rename target column to labels (Trainer expects this)
        dataset = dataset.rename_column("lumo", "labels")

        # Set format for PyTorch
        dataset.set_format(
            type="torch",
            columns=["input_ids", "attention_mask", "labels"],
        )

        # 5. Setup training arguments
        training_args = TrainingArguments(
            output_dir="./finetuned_model",
            num_train_epochs=10,
            per_device_train_batch_size=32,
            per_device_eval_batch_size=32,
            learning_rate=1e-5,
            warmup_ratio=0.1,
            logging_steps=10,
            eval_strategy="epoch",
            save_strategy="epoch",
            load_best_model_at_end=True,
            metric_for_best_model="eval_loss",
            save_total_limit=2,
            report_to="none",  # Disable wandb/tensorboard
        )

        # 6. Create data collator
        data_collator = DataCollatorWithPadding(tokenizer)

        # 7. Create Trainer
        trainer = Trainer(
            model=model,
            args=training_args,
            train_dataset=dataset["train"],
            eval_dataset=dataset["test"],
            processing_class=tokenizer,
            data_collator=data_collator,
        )

        # 8. Train!
        print("Starting training...")
        trainer.train()

        # 9. Save final model
        print("Saving model...")
        trainer.save_model("./finetuned_model")
        tokenizer.save_pretrained("./finetuned_model")
        print("Done!")

        # 10. Inference example with the finetuned model
        test_smiles = [
            "CCO",
            "CC(C)O",
            "CCCC",
            "c1ccccc1",
            "CC(=O)O"
        ]

        # Kekulize and tokenize
        kekulized_test = [kekulize_smiles(s) for s in test_smiles]
        inputs = tokenizer(
            kekulized_test,
            padding=True,
            truncation=True,
            max_length=512,
            return_tensors="pt"
        )

        # Move to device and run inference
        device = next(model.parameters()).device
        inputs = {k: v.to(device) for k, v in inputs.items()}

        model.eval()
        with torch.no_grad():
            outputs = model(**inputs)
            predictions = outputs["logits"].squeeze(-1).cpu()

        print("\nPredictions:")
        for smiles, pred in zip(test_smiles, predictions):
            print(f"  {smiles} → {pred.item():.4f} Hartree")


      if __name__ == "__main__":
          main()
  ```
</details>

### Use and Restrictions

Model weights are provided as-is for research purposes only, without guarantees of correctness, fitness for purpose, or warranties of any kind.
- Research use only
- No redistribution without permission
- No commercial use without licensing agreement

## Training Details

### Training Data

We use the the Enamine REAL Space dataset to pretrain MIST models.
At time of writing, Enamine REAL Space is the largest database of commercially available compounds.
The dataset was constructed using *forward synthetic analysis*: experimentally validated building blocks were converted into synthons annotated with reactivity features.
Enamine REAL Space was selected as the pretraining dataset since it was the largest database of molecular SMILES at the time of training, it is easily accessible for academic use and molecules relevant to downstream tasks, such as drug candidates, electrolytes, fragrances, live in synthetically accessible regions of chemical space.


### Training Procedure

<!-- This relates heavily to the Technical Specifications. Content here should link to that section when it is relevant to the training procedure. -->

### Inputs

The input to MIST models are SMILES strings for molecules.
This model was pretrained and on Kekulized SMILES strings.

### Model Architecture and Objective

- Encoder: [``RoBERTa-PreLayerNorm``](https://huggingface.co/docs/transformers/en/model_doc/roberta-prelayernorm) encoder with 8 layers, a hidden size of 512, intermediate size of 2048, 8 attention heads and maximum sequence length of 2048.
- Objective: MLM (Masked Language Modeling)
- Loss:  Cross-Entropy Loss
- Optimizer: ``deepspeed.ops.lamb.FusedLAMB``

### Compute Infrastructure

#### Hardware
This model was pre-trained on 2 NVIDIA A100-SXM4-80GB in ~12 hours and 15 minutes.

#### Software

This model was trained with PyTorchLightning using the DeepSpeed strategy for data distributed parallelism.
Model are exported in a Safetensors format.

## Citation

<!-- If there is a paper or blog post introducing the model, the APA and Bibtex information for that should go in this section. -->

If you use this model in your research, please cite:

```bibtex
@online{MIST,
  title = {Foundation Models for Discovery and Exploration in Chemical Space},
  author = {Wadell, Alexius and Bhutani, Anoushka and Azumah, Victor and Ellis-Mohr, Austin R. and Kelly, Celia and Zhao, Hancheng and Nayak, Anuj K. and Hegazy, Kareem and Brace, Alexander and Lin, Hongyi and Emani, Murali and Vishwanath, Venkatram and Gering, Kevin and Alkan, Melisa and Gibbs, Tom and Wells, Jack and Varshney, Lav R. and Ramsundar, Bharath and Duraisamy, Karthik and Mahoney, Michael W. and Ramanathan, Arvind and Viswanathan, Venkatasubramanian},
  date = {2025-10-20},
  eprint = {2510.18900},
  eprinttype = {arXiv},
  eprintclass = {physics},
  doi = {10.48550/arXiv.2510.18900},
  url = {http://arxiv.org/abs/2510.18900},  
}
```

## Model Card Authors
Anoushka Bhutani, Alexius Wadell

## Model Card Contact

For questions, issues, or licensing inquiries, please contact Venkat Viswanathan [venkvis@umich.edu](mailto:venkvis@umich.edu).