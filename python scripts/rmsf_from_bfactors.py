#!/usr/bin/env python3
"""
batch_rmsf_from_cif.py

Batch compute per-residue Cα RMSF between experimental apo and holo structures,
attach pLDDT scores from AlphaFold predictions, and save a single combined CSV
across all pairs (with apo & holo IDs as columns). Figures are generated per pair.

Experimental structures may have arbitrary chain IDs; AlphaFold models use chain 'A'.
Chain mappings (experimental → chosen chain) are provided in CHAIN_MAP_CSV.

Output CSV columns: apo_id, holo_id, chain, residue, rmsf, plddt_apo, plddt_holo
PNG figures per pair are unchanged.

Dependencies:
    pip install biopython pandas numpy matplotlib

Usage:
    python batch_rmsf_from_cif.py
"""

import os
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, Superimposer
import matplotlib.pyplot as plt

# --- User parameters ---
CIF_DIR       = './structures_from_pdb/content/cif_files'
ALL_FOLDS     = './all_folds'
PAIRS_FILE    = './reviewed_apo_holo_pair_ligand_and_ion_data.tsv'
CHAIN_MAP_CSV = './apo_holo_pairs.csv'
MODEL_W_IDS   = './model_w_pdb_ids.csv'
OUTPUT_DIR    = './rmsf_outputs_apo_holo'
COMBINED_CSV  = os.path.join(OUTPUT_DIR, 'combined_apo_holo_rmsf.csv')

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load chain mapping
chain_df = pd.read_csv(CHAIN_MAP_CSV)
chain_df['holo_key'] = chain_df['holo_structure'].str[:4]
chain_map = chain_df.set_index('holo_key')['holo_aligned_chain'].to_dict()

# Load pairs and model IDs
pairs_df    = pd.read_csv(PAIRS_FILE, sep='\t', usecols=['Apo_PDB_ID','Holo_PDB_ID'])
model_id_df = pd.read_csv(MODEL_W_IDS)
parser      = MMCIFParser(QUIET=True)

# Helpers
def ca_atoms(model, chain='A'):
    atoms = {}
    if chain not in model:
        return atoms
    for res in model[chain]:
        if 'CA' in res:
            atoms[res.id[1]] = res['CA']
    return atoms

def residue_plddt(model, chain='A'):
    scores = {}
    if chain not in model:
        return scores
    for res in model[chain]:
        num = res.id[1]
        if 'CA' in res:
            scores[num] = float(res['CA'].get_bfactor())
    return scores

def compute_rmsf(res_nums, atoms1, atoms2):
    sup = Superimposer()
    fixed  = [atoms1[i] for i in res_nums]
    moving = [atoms2[i] for i in res_nums]
    sup.set_atoms(fixed, moving)
    sup.apply(moving)
    vals = []
    for a1, a2 in zip(fixed, moving):
        c1 = a1.get_coord()
        c2 = a2.get_coord()
        mid = 0.5*(c1 + c2)
        vals.append(np.sqrt(((c1-mid)**2 + (c2-mid)**2).sum()/2.0))
    return vals

# Collect results
all_results = []

for idx, row in pairs_df.iterrows():
    apo_id, holo_id = row['Apo_PDB_ID'], row['Holo_PDB_ID']
    key = holo_id[:4]
    chain_id = chain_map.get(key)
    if pd.isna(chain_id) or chain_id is None:
        print(f"[{idx}] No chain mapping for {holo_id}, skipping")
        continue
    chain_id = str(chain_id).split(',')[0].strip()

    # File paths
    apo_cif     = os.path.join(CIF_DIR,    f"{apo_id}.cif")
    holo_cif    = os.path.join(CIF_DIR,    f"{holo_id}.cif")
    apo_model   = model_id_df.loc[model_id_df['apo_pdb_id']==apo_id, 'alphafold_model_id'].iloc[0]
    holo_model  = model_id_df.loc[model_id_df['holo_pdb_id']==holo_id, 'alphafold_model_id'].iloc[0]
    af_apo_cif  = os.path.join(ALL_FOLDS,  f"{apo_id}_apo",  f"fold_{apo_id}_apo_model_{apo_model}.cif")
    af_holo_cif = os.path.join(ALL_FOLDS,  f"{holo_id}_holo", f"fold_{holo_id}_holo_model_{holo_model}.cif")

    # Check existence
    missing = [p for p in (apo_cif, holo_cif, af_apo_cif, af_holo_cif) if not os.path.isfile(p)]
    if missing:
        print(f"[{idx}] Skipping {apo_id}-{holo_id}: missing {missing}")
        continue

    # Load structures
    apo_exp   = parser.get_structure(f"{apo_id}_exp",  apo_cif)[0]
    holo_exp  = parser.get_structure(f"{holo_id}_exp", holo_cif)[0]
    af_apo    = parser.get_structure(f"{apo_id}_af",   af_apo_cif)[0]
    af_holo   = parser.get_structure(f"{holo_id}_af",  af_holo_cif)[0]

    # Extract data
    plddt_apo  = residue_plddt(af_apo,  chain='A')
    plddt_holo = residue_plddt(af_holo, chain='A')
    ca_apo_exp  = ca_atoms(apo_exp,  chain_id)
    ca_holo_exp = ca_atoms(holo_exp, chain_id)
    common = sorted(set(ca_apo_exp) & set(ca_holo_exp))
    if not common:
        print(f"[{idx}] No matching Cα for {apo_id}-{holo_id}")
        continue

    rmsf_vals      = compute_rmsf(common, ca_apo_exp, ca_holo_exp)
    for res, rmsf in zip(common, rmsf_vals):
        all_results.append({
            'apo_id':    apo_id,
            'holo_id':   holo_id,
            'chain':     chain_id,
            'residue':   res,
            'rmsf':   rmsf,
            'plddt_apo':  plddt_apo.get(res),
            'plddt_holo': plddt_holo.get(res)
        })

    # Generate unchanged figure\ n    
    df_tmp = pd.DataFrame({'residue':common, 'rmsf':rmsf_vals})
    plt.figure(figsize=(8,4))
    plt.plot(df_tmp['residue'], df_tmp['rmsf'], marker='o', linestyle='-')
    plt.xlabel('Residue Number')
    plt.ylabel('RMSF (Å)')
    plt.title(f"Exp apo–holo RMSF: {apo_id}-{holo_id} ({chain_id})")
    plt.grid(True)
    plt.tight_layout()
    png = os.path.join(OUTPUT_DIR, f"{apo_id}_{holo_id}_{chain_id}_rmsf.png")
    plt.savefig(png, dpi=300)
    plt.close()
    print(f"[{idx}] Done {apo_id}-{holo_id}")

# Save combined CSV
df_all = pd.DataFrame(all_results)
df_all.to_csv(COMBINED_CSV, index=False)
print("Combined CSV written to", COMBINED_CSV)
