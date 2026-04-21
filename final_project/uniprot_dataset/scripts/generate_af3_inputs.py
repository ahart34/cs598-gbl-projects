import json, os

with open("../raw/uniprot_proteins_with_seq.json") as f:
    all_seqs = json.load(f)

with open("../raw/filtered_proteins.txt") as f:
    filtered = [l.strip() for l in f if l.strip()]

af3_input_dir = "/media/Data/shanezhu/alphafold3_install/alphafold3/input/uniprot"
os.makedirs(af3_input_dir, exist_ok=True)

generated = 0
skipped = 0

for uid in filtered:
    val = all_seqs.get(uid)
    if not val:
        print(f"No entry found for {uid}")
        skipped += 1
        continue

    seq = None
    for item in val:
        if isinstance(item, dict) and 'sequence' in item:
            seq = item['sequence']
            break

    if not seq:
        print(f"Skipping {uid}: no sequence found")
        skipped += 1
        continue

    inp = {
        "name": uid,
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": seq
                }
            }
        ],
        "modelSeeds": [42],
        "dialect": "alphafold3",
        "version": 1
    }
    with open(f"{af3_input_dir}/{uid}.json", "w") as f:
        json.dump(inp, f)
    generated += 1

print(f"Generated: {generated}, Skipped: {skipped}")
