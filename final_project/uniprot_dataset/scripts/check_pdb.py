import json, requests, time

raw_dir = "../raw"

with open(f"{raw_dir}/uniprot_proteins_with_seq.json") as f:
    all_seqs = json.load(f)

with open(f"{raw_dir}/filtered_proteins.txt") as f:
    filtered = [l.strip() for l in f if l.strip()]

seqs = {k: v for k, v in all_seqs.items() if k in filtered}
print(f"Total filtered proteins: {len(seqs)}")

has_pdb = {}   # uid -> [pdb_ids]
no_pdb = []

for uid in seqs:
    url = f"https://www.ebi.ac.uk/proteins/api/proteins/{uid}"
    r = requests.get(url, headers={"Accept": "application/json"})
    if r.ok:
        entry = r.json()
        pdbs = [x['id'] for x in entry.get('dbReferences', []) if x['type'] == 'PDB']
        if pdbs:
            has_pdb[uid] = pdbs
        else:
            no_pdb.append(uid)
    else:
        no_pdb.append(uid)
    time.sleep(0.2)  # avoid rate limit

print(f"Has PDB: {len(has_pdb)}")
print(f"Needs AF3 folding: {len(no_pdb)}")

with open(f"{raw_dir}/has_pdb.json", "w") as f:
    json.dump(has_pdb, f, indent=2)

with open(f"{raw_dir}/no_pdb.txt", "w") as f:
    f.write("\n".join(no_pdb))
