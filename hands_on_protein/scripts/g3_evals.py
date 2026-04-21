import os
import shutil
import argparse
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from omegaconf import OmegaConf
from EvalRunner import EvalRunner
from PSEA_eval import calc_psea_metrics

class DictToObject:
    def __init__(self, dictionary):
        for key, value in dictionary.items():
            setattr(self, key, value)

def load_config(path):
    return OmegaConf.to_container(OmegaConf.load(path), resolve=True)

def find_inference_dirs(model_dir):
    # find directories in the inference folders which have .pdbs (i.e., not log files or csvs)
    inference_dirs = []
    
    for p in sorted(model_dir.iterdir()):
        if not p.is_dir():
            continue
        if p.name.startswith("."):
            continue
        if any(f.suffix == ".pdb" for f in p.iterdir() if f.is_file()):
            inference_dirs.append(p)
    
    return inference_dirs

def find_main_pdb(inference_dir):
    # filter out *.pdb_esm files and sample.pdb from a possible previous run
    pdbs = sorted(p for p in inference_dir.glob("*.pdb") if p.name != "sample.pdb")
    
    if not pdbs:
        raise FileNotFoundError(f"no *.pdb found in {inference_dir}")
    
    return pdbs[0]

def create_sample_input_dir(inference_dir, pdb_file):
    # create the dir to pass on for EvalRunner.py using the example naming of 'sample.pdb'
    # we are moving the .pdb we care for in inference dir to a cleaner dir
    input_dir = inference_dir / "g3_eval_input"
    input_dir.mkdir(exist_ok=True)
    sample_pdb = input_dir / "sample.pdb"
    shutil.copy2(pdb_file, sample_pdb)
    
    return input_dir, sample_pdb
    
def run_designability(EvalRunner_conf, input_dir, output_dir):
    # run self consistency metrics
    selfcon_dir = output_dir / "self_consistency"
    selfcon_dir.mkdir(parents=True, exist_ok=True)
    selfcon_results = EvalRunner_conf.calc_designability(str(selfcon_dir), str(input_dir))
    
    return selfcon_dir, selfcon_results

def run_all_metrics(EvalRunner_conf, input_dir, sample_pdb, output_dir, pdb_i):
    # run all metrics
    selfcon_dir, selfcon_results = run_designability(EvalRunner_conf, input_dir, output_dir)
    mdtraj_results = EvalRunner_conf.calc_mdtraj_metrics(str(sample_pdb))
    novelty = EvalRunner_conf.pdbTM(str(sample_pdb), pdb_i)
    
    df_row = {
        "pdbTM": novelty,
        "best_tm_score": float(selfcon_results["tm_score"].max()),
        "best_bb_rmsd": float(selfcon_results["bb_rmsd"].min()),
        **mdtraj_results,
    }
    
    summary_df = pd.DataFrame([df_row])
    summary_df.to_csv(output_dir / "all_metrics.csv", index=False)
    
    return selfcon_dir, selfcon_results, summary_df

def run_psea(pdb_file, output_dir):
    # run psea for structure %
    psea_metrics = calc_psea_metrics(str(pdb_file))
    summary_df = pd.DataFrame([psea_metrics])
    summary_df.to_csv(output_dir / "psea_metrics.csv", index=False)
    
    return summary_df

def update_summary_wpsea(model_dir, psea_df):
    # add psea results to g3_eval_all_metrics_summary.csv
    all_metrics_csv = model_dir / "g3_eval_all_metrics_summary.csv"
    
    if not all_metrics_csv.exists():
        return
    existing_df = pd.read_csv(all_metrics_csv)
    
    psea_cols = [
        "design_name",
        "pdb_file",
        "non_coil_percent",
        "coil_percent",
        "helix_percent",
        "strand_percent",
    ]
    
    # avoid pd getting confused if we had psea metrics already
    for col in ["non_coil_percent", "coil_percent", "helix_percent", "strand_percent"]:
        if col in existing_df.columns:
            existing_df = existing_df.drop(columns=col)
    
    merged_df = existing_df.merge(
        psea_df[psea_cols],
        on=["design_name", "pdb_file"],
        how="left",
    )
    merged_df.to_csv(all_metrics_csv, index=False)

def run_diversity(EvalRunner_conf, model_dir, inference_dirs):
    # run diversity with maxcluster
    cluster_dir = model_dir / "cluster"
    cluster_dir.mkdir(parents=True, exist_ok=True)
    pdb_paths = []
    
    prog_bar = tqdm(inference_dirs, desc="setting pdbs for clustering", unit="pdb", total=len(inference_dirs))
    for inference_dir in prog_bar:
        prog_bar.set_postfix_str(inference_dir.name)
        pdb_file = find_main_pdb(inference_dir)
        pdb_paths.append(str(pdb_file.resolve()))
    
    good_path_clustering = cluster_dir / "designable_paths.txt"
    with open(good_path_clustering, "w") as f:
        f.write("\n".join(pdb_paths))
    clusters = EvalRunner_conf.run_max_cluster(str(good_path_clustering), str(cluster_dir))
    
    summary_df = pd.DataFrame(
        [
            {
                "model_dir": str(model_dir),
                "num_designs": int(len(pdb_paths)),
                "num_clusters": int(clusters),
                "cluster_dir": str(cluster_dir),
            }
        ]
    )
    output_f = model_dir / "g3_eval_diversity_summary.csv"
    summary_df.to_csv(output_f, index=False)
    
    return summary_df, output_f

def main():
    parser = argparse.ArgumentParser(description="run g3's evals over one model inference folder")
    parser.add_argument(
        "--config",
        default="./configs/evaluation.yaml",
        help="path to evaluation config yaml",
    )
    parser.add_argument(
        "--pdb-path",
        required=True,
        help="path to one model inference folder (e.g. protein_inference/inference_ucond_200m_notri_vf)",
    )
    parser.add_argument(
        "--mode",
        choices=["designability", "all_metrics", "psea", "diversity"],
        default="all_metrics",
        help="options: designability, all metrics, psea, diversity",
    )
    args = parser.parse_args()
    
    config_dict = load_config(args.config)
    conf = DictToObject(config_dict)
    model_dir = Path(args.pdb_path).resolve()
    inference_dirs = find_inference_dirs(model_dir)
    
    EvalRunner_conf = None
    if args.mode in ["designability", "all_metrics", "diversity"]:
        EvalRunner_conf = EvalRunner(conf)
    if args.mode == "diversity":
        summary_df, output_f = run_diversity(EvalRunner_conf, model_dir, inference_dirs)
        print(f"diversity summary saved to: {output_f}", flush=True)
        print(f"cluster outputs saved to: {model_dir / 'cluster'}", flush=True)
        return

    df_rows = []

    prog_bar = tqdm(inference_dirs, desc="running pdbs", unit="pdb", total=len(inference_dirs))
    for i, inference_dir in enumerate(prog_bar, start=1):
        prog_bar.set_postfix_str(inference_dir.name)
        pdb_file = find_main_pdb(inference_dir)
        input_dir, sample_pdb = create_sample_input_dir(inference_dir, pdb_file)
        output_dir = inference_dir / "g3_eval"
        output_dir.mkdir(exist_ok=True)

        if args.mode == "designability":
            selfcon_dir, selfcon_results = run_designability(EvalRunner_conf, input_dir, output_dir)
            df_rows.append(
                {
                    "design_name": inference_dir.name,
                    "pdb_file": str(pdb_file),
                    "g3_eval_dir": str(output_dir),
                    "self_consistency_dir": str(selfcon_dir),
                    "best_tm_score": float(selfcon_results["tm_score"].max()),
                    "best_bb_rmsd": float(selfcon_results["bb_rmsd"].min()),
                }
            )
        elif args.mode == "all_metrics":
            selfcon_dir, selfcon_results, summary_df = run_all_metrics(
                EvalRunner_conf, input_dir, sample_pdb, output_dir, i
            )
            df_row = summary_df.iloc[0].to_dict()
            df_row.update(
                {
                    "design_name": inference_dir.name,
                    "pdb_file": str(pdb_file),
                    "g3_eval_dir": str(output_dir),
                    "self_consistency_dir": str(selfcon_dir),
                }
            )
            df_rows.append(df_row)
        elif args.mode == "psea":
            summary_df = run_psea(pdb_file, output_dir)
            df_row = summary_df.iloc[0].to_dict()
            df_row.update(
                {
                    "design_name": inference_dir.name,
                    "pdb_file": str(pdb_file),
                    "g3_eval_dir": str(output_dir),
                }
            )
            df_rows.append(df_row)

    summary_df = pd.DataFrame(df_rows)

    if args.mode == "designability":
        output_f = model_dir / "g3_eval_designability_summary.csv"
        summary_df.to_csv(output_f, index=False)
    elif args.mode == "all_metrics":
        output_f = model_dir / "g3_eval_all_metrics_summary.csv"
        summary_df.to_csv(output_f, index=False)
    elif args.mode == "psea":
        output_f = model_dir / "g3_eval_psea_summary.csv"
        summary_df.to_csv(output_f, index=False)
        update_summary_wpsea(model_dir, summary_df)

    print(f"summary saved to: {output_f}", flush=True)

if __name__ == "__main__":
    main()