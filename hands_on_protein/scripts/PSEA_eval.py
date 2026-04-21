import biotite.structure as struc
import biotite.structure.io.pdb as pdb
from pathlib import Path
from collections import Counter

def calc_psea_metrics(pdb_path):
    # calcuclate structure % using only c_a backbone
    # we use biotite as its based on psea algorithm
    pdb_path = str(Path(pdb_path).resolve())
    pdb_file = pdb.PDBFile.read(pdb_path)
    ca_bckbn_atoms = pdb_file.get_structure(model=1)

    sse = struc.annotate_sse(ca_bckbn_atoms)
    # a=helix b=strand c=coil
    valid_sse = [x for x in sse.tolist() if x in ("a", "b", "c")]

    counts = Counter(valid_sse)
    n = len(valid_sse)

    # coil + non_coil = 1, other are subtypes
    helix_percent = counts.get("a", 0) / n
    strand_percent = counts.get("b", 0) / n
    coil_percent = counts.get("c", 0) / n
    non_coil_percent = helix_percent + strand_percent

    return {
        "non_coil_percent": float(non_coil_percent),
        "coil_percent": float(coil_percent),
        "helix_percent": float(helix_percent),
        "strand_percent": float(strand_percent),
    }