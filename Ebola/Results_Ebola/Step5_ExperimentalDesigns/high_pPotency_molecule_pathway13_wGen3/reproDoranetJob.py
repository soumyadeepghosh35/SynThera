#!/usr/bin/env python3
import ast, csv, os, shutil, sys, time
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

doranetPath = Path(r"/users/sghosh6/DTRA_project/MACAW/doranet")
sys.path.insert(0, str(doranetPath))

import doranet.modules.enzymatic as enzymatic
import doranet.modules.post_processing as postProcessing

jobName      = "high_pPotency_molecule_pathway13_wGen3"
starters     = {'Nc1nc2c(ncn2[C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)[nH]1'}
helpers      = {'[C-]#[O+]', 'O=S=O', 'O', 'O=S(O)O', '[H][H]', 'C=C', 'CO', 'N', 'S', 'Br', 'C#N', '[Br][Br]', 'O=C=O', 'N#N', 'O=[N+]([O-])O', 'O=NO', 'C=O', 'O=O', 'O=S(=O)(O)O', 'NO', 'N#CO'}
target       = {'Nc1nc2c(ncn2[C@@H]2OC(C(O)O)=C[C@H]2O)c(=O)[nH]1'}
maxAtoms     = {'C': 15, 'N': 8, 'O': 8, 'S': 1}
generations  = 3
ruleset      = "JN3604IMT"
searchDepth  = 3
maxNumRxns   = 3
minRxnAtomEconomy   = 0.0
filterExactNumRxns  = True
exactNumRxns        = 3
runRanking          = True
runVisualization    = False

def getSmilesProps(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return "N/A", 0.0, 0
    return rdMolDescriptors.CalcMolFormula(mol), round(rdMolDescriptors.CalcExactMolWt(mol), 4), mol.GetNumHeavyAtoms()

def splitPathwayBlocks(p):
    if not p.exists():
        return []
    lines = p.read_text(encoding="utf-8").splitlines()
    blocks, current = [], []
    for line in lines:
        if line.startswith("pathway number ") and current:
            blocks.append(current)
            current = [line]
        else:
            current.append(line)
    if current:
        blocks.append(current)
    return blocks

def countSteps(block):
    import ast
    prefix = "reaction SMILES stoichiometry "
    for line in block:
        if line.startswith(prefix):
            try:
                parsed = ast.literal_eval(line[len(prefix):].strip())
                return len(parsed) if isinstance(parsed, list) else None
            except Exception:
                return None
    return None

t0 = time.time()
network = enzymatic.generate_network(
    job_name=jobName, starters=starters, gen=generations, max_atoms=maxAtoms,
    direction="forward", targets=target, ruleset=ruleset,
)

smilesList = list(starters) + [m.uid for m in network.mols if m.uid not in starters]
with Path(f"{jobName}_molecules.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["SMILES", "isStarter", "molFormula", "molWeight", "numHeavyAtoms"])
    for s in smilesList:
        w.writerow([s, s in starters, *getSmilesProps(s)])

postProcessing.pretreat_networks(networks={network}, total_generations=generations, starters=starters, helpers=helpers, job_name=jobName)
postProcessing.pathway_finder(starters=starters, helpers=helpers, target=target, search_depth=searchDepth, max_num_rxns=maxNumRxns, min_rxn_atom_economy=minRxnAtomEconomy, job_name=jobName)

if filterExactNumRxns:
    pPath = Path(f"{jobName}_pathways.txt")
    ePath  = Path(f"{jobName}_pathways_exact{exactNumRxns}.txt")
    if pPath.exists():
        blocks = splitPathwayBlocks(pPath)
        exact  = [b for b in blocks if countSteps(b) == exactNumRxns]
        txt    = "\n\n".join("\n".join(b) for b in exact)
        ePath.write_text(txt + "\n" if txt else "", encoding="utf-8")
        pPath.write_text(txt + "\n" if txt else "", encoding="utf-8")
        print(f"Exact-step filter Gen{generations}: {len(exact)}/{len(blocks)} pathways retained")

if runRanking:
    postProcessing.pathway_ranking(starters=starters, helpers=helpers, target=target, job_name=jobName, num_process=1)
if runVisualization:
    postProcessing.pathway_visualization(starters=starters, helpers=helpers, job_name=jobName, num_process=1)

print(f"Done {jobName} in {time.time() - t0:.2f} s")
