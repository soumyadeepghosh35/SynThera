#!/usr/bin/env python

'''
### Use generated molecules as targets (find pathways to everything)
'''

import sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import csv
import time

DORANET_PATH = Path("/users/sghosh6/DTRA_project/MACAW/doranet")
sys.path.insert(0, str(DORANET_PATH))

import doranet.modules.enzymatic as enzymatic
import doranet.modules.post_processing as post_processing

start_time = time.time()

# Starter molecule
user_starters = {'CCCCCCCCCCCCCCCCOCCCOP(=O)(CO[C@@H](CN1C=CC(=NC1=O)N)CO)O'}

user_helpers = {
    'O', 'O=O', '[H][H]', 'O=C=O', 'C=O', '[C-]#[O+]', 'Br', '[Br][Br]', 'CO',
    'C=C', 'O=S(O)O', 'N', 'O=S(=O)(O)O', 'O=NO', 'N#N', 'O=[N+]([O-])O', 'NO',
    'C#N', 'S', 'O=S=O', 'N#CO'
}

job_name = "BrincidofovirOnly_gen2"

print(f"Starter: {list(user_starters)[0]}")

# Generate network
forward_network = enzymatic.generate_network(
    job_name=job_name,
    starters=user_starters,
    gen=2,
    max_atoms={'C': 41, 'N': 5, 'O': 11, 'S': 0},
    direction="forward",
    ruleset="JN3604IMT"
)

# Extract generated molecules
#smiles_list = [mol.uid for mol in forward_network.mols]
#print(f"\nGenerated {len(smiles_list)} molecules")

# Ensure starters are included in the molecule list
smiles_list = list(user_starters) + [mol.uid for mol in forward_network.mols if mol.uid not in user_starters]
print(f"\nGenerated {len(smiles_list) - len(user_starters)} new molecules + {len(user_starters)} starters")

# Use ALL generated molecules (except starters and helpers) as targets
all_targets = set(smiles_list) - user_starters - user_helpers

print(f"Using {len(all_targets)} generated molecules as targets for pathway finding")

# Save molecules to CSV
output_path = Path(f"{job_name}_molecules.csv")
with open(output_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['SMILES', 'Is_Starter', 'MolFormula', 'MolWeight', 'NumHeavyAtoms'])
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            formula = rdMolDescriptors.CalcMolFormula(mol)
            mol_weight = round(rdMolDescriptors.CalcExactMolWt(mol), 4)
            num_heavy = mol.GetNumHeavyAtoms()
        else:
            formula, mol_weight, num_heavy = "N/A", 0, 0
        writer.writerow([smi, smi in user_starters, formula, mol_weight, num_heavy])

print(f"\nSaved molecules to {output_path}")

# Generate pathways to all generated molecules
if all_targets:
    post_processing.one_step(
        networks={forward_network},
        total_generations=2,
        starters=user_starters,
        helpers=user_helpers,
        target=all_targets,  # Use generated molecules as targets
        job_name=job_name,
    )
    print(f"\nPathway files generated with prefix: {job_name}")
else:
    print("No targets found for pathway generation")

print(f"Time: {time.time() - start_time:.2f} s")
