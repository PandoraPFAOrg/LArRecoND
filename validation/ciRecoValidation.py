import argparse
import json
import sys
from pathlib import Path

import awkward as ak
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import uproot

matplotlib.use("Agg")

PLOT_BRANCHES = (
    ("n3DHits", "3D hits per PFP"),
    ("trackScore", "Track score"),
    ("energy", "Reconstructed energy"),
    ("purity", "MC match purity"),
    ("completeness", "MC match completeness"),
)


def describe(values):
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {
            "count": 0,
            "min": None,
            "max": None,
            "mean": None,
            "std": None,
            "median": None,
        }
    return {
        "count": int(values.size),
        "min": float(values.min()),
        "max": float(values.max()),
        "mean": float(values.mean()),
        "std": float(values.std()),
        "median": float(np.median(values)),
    }

def populate_stats(file_path):
    if file_path is None:
        return

    file_path = Path(file_path)
    if not file_path.is_file() or file_path.stat().st_size == 0:
        print(f"File is missing or empty: {file_path}")
        return

    with uproot.open(file_path) as root_file:
        if "LArRecoND" not in root_file:
            print(f"Missing LArRecoND tree in reference file: {file_path}")
            return

        tree = root_file["LArRecoND"]
        cluster_ids = tree["clusterId"].array(library="ak")
        pfps = ak.to_numpy(ak.flatten(cluster_ids, axis=1))
        branch_values = {
            branch: ak.to_numpy(ak.flatten(tree[branch].array(library="ak"), axis=None))
            for branch, _ in PLOT_BRANCHES
        }

        return pfps, branch_values


def validate_and_plot(input_path, dataset_id, dataset_name, output_dir, reference_path):
    input_path = Path(input_path)
    output_dir = Path(output_dir)

    if not input_path.is_file() or input_path.stat().st_size == 0:
        raise ValueError(f"Output file is missing or empty: {input_path}")

    with uproot.open(input_path) as root_file:
        if "LArRecoND" not in root_file:
            raise ValueError(f"Missing LArRecoND tree in {input_path}")

        tree = root_file["LArRecoND"]
        entries = int(tree.num_entries)
        if entries == 0:
            raise ValueError(f"LArRecoND tree has no entries: {input_path}")

        branches = set(tree.keys())
        required_branches = {"clusterId", *(branch for branch, _ in PLOT_BRANCHES)}
        missing_branches = sorted(required_branches - branches)
        if missing_branches:
            raise ValueError(
                f"Missing required branches: {', '.join(missing_branches)}"
            )

        cluster_ids = tree["clusterId"].array(library="ak")
        pfps_per_event = ak.to_numpy(ak.num(cluster_ids, axis=1))
        pfp_total = int(np.sum(pfps_per_event))
        if pfp_total == 0:
            raise ValueError(
                f"LArRecoND tree has entries but no reconstructed PFPs: {input_path}"
            )

        branches_to_read = [branch for branch, _ in PLOT_BRANCHES]
        branch_values = {
            branch: ak.to_numpy(ak.flatten(tree[branch].array(library="ak"), axis=None))
            for branch in branches_to_read
        }

    # Next, read in the reference file if provided and compare the distributions.
    ref_pfps = None
    ref_branch_values = {}
    ref_stats = populate_stats(reference_path)

    if ref_stats is not None:
        ref_pfps, ref_branch_values = ref_stats

    event_stats = describe(pfps_per_event)
    branch_stats = [
        {"branch": branch, **describe(branch_values[branch])}
        for branch in branches_to_read
    ]

    output_dir.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.flatten()

    if ref_pfps is not None:
        axes[0].hist(ref_pfps, bins="auto", color="gray", alpha=0.5, label="Main")

    axes[0].hist(pfps_per_event, bins="auto", color="#2878a5", edgecolor="white")
    axes[0].set_title("Reconstructed PFPs per event")
    axes[0].set_xlabel("PFP count")
    axes[0].set_ylabel("Events")

    if ref_pfps is not None:
        axes[0].legend()

    for axis, (branch, title) in zip(axes[1:], PLOT_BRANCHES):
        values = np.asarray(branch_values[branch], dtype=np.float64)
        values = values[np.isfinite(values)]

        if branch in ref_branch_values:
            ref_values = np.asarray(ref_branch_values[branch], dtype=np.float64)
            ref_values = ref_values[np.isfinite(ref_values)]
            axis.hist(
                ref_values, bins=50, color="gray", alpha=0.5, label="Reference"
            )

        if values.size:
            axis.hist(values, bins=50, color="#d06b3c", edgecolor="white")
        else:
            axis.text(
                0.5,
                0.5,
                "No values",
                ha="center",
                va="center",
                transform=axis.transAxes,
            )
        axis.set_title(title)
        axis.set_xlabel(branch)
        axis.set_ylabel("PFPs")

        if branch in ref_branch_values:
            axis.legend()

    figure.suptitle(f"{dataset_name}: {entries} events, {pfp_total:,} PFPs")
    figure.tight_layout()
    plot_path = output_dir / f"{dataset_id}.png"
    figure.savefig(plot_path, dpi=150)
    plt.close(figure)

    summary = {
        "dataset_id": dataset_id,
        "dataset_name": dataset_name,
        "entries": entries,
        "pfp_total": pfp_total,
        "pfps_per_event": [int(value) for value in pfps_per_event],
        "event_stats": event_stats,
        "branch_stats": branch_stats,
        "plot": plot_path.name,
    }
    summary_path = output_dir / f"{dataset_id}.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(
        description="Validate and plot a LArRecoND ROOT output file."
    )
    parser.add_argument("--input", required=True, help="LArRecoND.root path")
    parser.add_argument("--dataset-id", required=True, help="Stable dataset identifier")
    parser.add_argument(
        "--dataset-name", required=True, help="Human-readable dataset name"
    )
    parser.add_argument(
        "--output-dir", required=True, help="Directory for JSON summary and plot"
    )
    parser.add_argument(
        "--reference",
        required=False,
        help="Path to reference ROOT file for validation",
        default=None,
    )
    args = parser.parse_args()

    try:
        validate_and_plot(
            args.input,
            args.dataset_id,
            args.dataset_name,
            args.output_dir,
            args.reference,
        )
    except Exception as error:
        print(f"CI validation failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
