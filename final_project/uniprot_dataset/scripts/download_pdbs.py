# download_pdbs.py
import json, requests, os

with open("../raw/has_pdb.json") as f:
    has_pdb = json.load(f)

out_dir = "../structures/from_pdb"
os.makedirs(out_dir, exist_ok=True)

for uid, pdb_ids in has_pdb.items():
    pdb_id = pdb_ids[0]  # 1st PDB
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    r = requests.get(url)
    if r.ok:
        with open(f"{out_dir}/{uid}_{pdb_id}.cif", "wb") as f:
            f.write(r.content)
        print(f"Downloaded {uid} -> {pdb_id}")
    else:
        print(f"Failed: {uid} {pdb_id}")
