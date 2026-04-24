import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

BASE_DIR = Path(__file__).parent / "evals_proteina"
OUTPUT_DIR = Path(__file__).parent / "plots"
OUTPUT_DIR.mkdir(exist_ok=True)

METRIC_COLS = [
    "pdbTM", "best_bb_rmsd", "diversity",
    "non_coil_percent", "coil_percent", "helix_percent", "strand_percent",
]

METRIC_COLS_PER_LENGTH = [m for m in METRIC_COLS if m != "diversity"]

METRIC_TRANSFORMS = {
    "pdbTM":        lambda x: 1 - x,
    "best_bb_rmsd": lambda x: 1 / x,
}

METRIC_LABELS = {
    "pdbTM":            "Novelty: 1 - pdbTM",
    "best_bb_rmsd":     "Designability: 1 / bbRMSD",
    "diversity":        "Diversity: # Clusters",
    "non_coil_percent": "Secondary Structure: Non-Coil",
    "coil_percent":     "Secondary Structure: Coil",
    "helix_percent":    "Secondary Structure: Helix",
    "strand_percent":   "Secondary Structure: Strand",
}

MODEL_LABELS = {
    "ucond_200m_notri_sc_noise045": "Stochastic: Noise 0.45",
    "ucond_200m_notri_sc_noise060": "Stochastic: Noise 0.60",
    "ucond_200m_notri_vf":          "Deterministic",
}

COLORS = plt.cm.Set2(np.linspace(0, 1, 8))

# Matplotlib defaults tuned for publication
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 14,
    "axes.titlesize": 14,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})


def load_data() -> dict[str, pd.DataFrame]:
    dfs = {}
    for folder in sorted(BASE_DIR.iterdir()):
        if not folder.is_dir():
            continue
        metrics_csv = next(
            (f for f in folder.glob("*.csv") if "metrics" in f.name), None
        )
        if metrics_csv is None:
            continue
        df = pd.read_csv(metrics_csv)
        df["seq_length"] = df["design_name"].str.extract(r"n_(\d+)_id_").astype(int)

        diversity_csv = folder / "g3_eval_diversity_summary.csv"
        if diversity_csv.exists():
            div = pd.read_csv(diversity_csv)
            df["diversity"] = int(div["num_clusters"].iloc[0])
        else:
            df["diversity"] = float("nan")

        label = MODEL_LABELS.get(folder.name, folder.name)
        dfs[label] = df
    return dfs


def plot_bar(dfs: dict[str, pd.DataFrame]) -> None:
    models = list(dfs.keys())
    n_models = len(models)
    ncols = 4
    nrows = (len(METRIC_COLS) + ncols - 1) // ncols
    width = 0.7 / n_models
    x = np.arange(n_models)

    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 3 * nrows))
    axes = axes.flatten()

    handles = []
    for idx, metric in enumerate(METRIC_COLS):
        ax = axes[idx]
        for i, (model, df) in enumerate(dfs.items()):
            values = METRIC_TRANSFORMS.get(metric, lambda x: x)(df[metric])
            mean = values.mean()
            std = values.std()
            yerr = [[min(std, mean)], [std]]  # don't allow value before 0 for error bar (wouldn't make sense)
            bar = ax.bar(i, mean, width * n_models * 0.9, yerr=yerr, label=model,
                         color=COLORS[i], capsize=3, alpha=0.88,
                         ecolor="black", linewidth=0.4)
            if idx == 0:
                handles.append(bar)

        ax.set_title(METRIC_LABELS[metric])
        ax.set_ylabel("Mean ± SD")
        ax.set_xticks([])
        ax.set_xlim(-0.6, n_models - 0.4)

    for idx in range(len(METRIC_COLS), len(axes)):
        axes[idx].set_visible(False)

    labels = list(dfs.keys())
    fig.legend(handles, labels, loc="lower center",
               ncol=n_models, bbox_to_anchor=(0.5, -0.04), fontsize=16)
    plt.tight_layout()
    out = OUTPUT_DIR / "bar_metrics_by_model.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


def plot_lines(dfs: dict[str, pd.DataFrame]) -> None:
    n_metrics = len(METRIC_COLS_PER_LENGTH)
    ncols = 4
    nrows = (n_metrics + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 3 * nrows))
    axes = axes.flatten()

    handles, labels = [], []
    for idx, metric in enumerate(METRIC_COLS_PER_LENGTH):
        ax = axes[idx]
        for i, (model, df) in enumerate(dfs.items()):
            transform = METRIC_TRANSFORMS.get(metric, lambda x: x)
            grouped = df.groupby("seq_length")[metric].apply(lambda g: transform(g).mean()).sort_index()
            line, = ax.plot(grouped.index, grouped.values, marker="o",
                            label=model, color=COLORS[i], linewidth=1.6, markersize=3.5)
            if idx == 0:
                handles.append(line)
                labels.append(model)
        ax.set_title(METRIC_LABELS[metric])
        ax.set_xlabel("Sequence Length")
        ax.set_ylabel("Mean Value")

    for idx in range(len(METRIC_COLS_PER_LENGTH), len(axes)):
        axes[idx].set_visible(False)

    fig.legend(handles, labels, loc="lower center",
               ncol=len(dfs), bbox_to_anchor=(0.5, -0.04), fontsize=16)
    plt.tight_layout()
    out = OUTPUT_DIR / "line_metrics_by_length.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


if __name__ == "__main__":
    dfs = load_data()
    for model, df in dfs.items():
        lengths = sorted(df["seq_length"].unique())

    plot_bar(dfs)
    plot_lines(dfs)

