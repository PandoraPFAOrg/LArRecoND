import argparse
import json
import sys
from pathlib import Path

import awkward as ak
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import uproot

matplotlib.use("Agg")

PLOT_BRANCHES = (
    ("n3DHits", "3D hits per PFO"),
    ("trackScore", "Track score"),
    ("energy", "Reconstructed energy"),
    ("purity", "MC match purity"),
    ("completeness", "MC match completeness"),
)


N_BINS = 50
HEAVY_TAILED = ("n3DHits", "energy")  # robust range and log bins
LOG_SPAN = 100  # use log bins when the (robust) range spans more than this
MIN_RATIO_COUNTS = 5  # hide ratio points with fewer counts than this
CURRENT_COLOUR = "#d06b3c"
REFERENCE_COLOURS = ("#2878a5", "#4c9a5a", "#8e6bb0", "#b8a02a", "#6b6b6b")
REFERENCE_STYLES = ("-", "--", ":", "-.")


def finite(values):
    values = np.asarray(values, dtype=np.float64)
    return values[np.isfinite(values)]


def make_bins(values, ref_arrays=(), integer=False, heavy_tail=False):
    """Common bin edges for the current and all reference values.

    Bounded quantities use the full range. For heavy_tail quantities the range
    is robust to outliers (0.5-99.5 percentiles of the pooled values; values
    outside it are put in the edge bins by hist()), and strictly positive ones
    get log-spaced bins, as a linear axis would squash them into one bin.
    """
    arrays = [a for a in (values, *ref_arrays) if a is not None and a.size]
    if not arrays:
        return np.linspace(0.0, 1.0, N_BINS + 1)

    pooled = np.concatenate(arrays)
    low, high = pooled.min(), pooled.max()
    if integer and high - low < N_BINS:
        return np.arange(np.floor(low), np.floor(high) + 2) - 0.5

    if heavy_tail:
        low, high = np.percentile(pooled, [0.5, 99.5])
    if integer:
        low, high = np.floor(low), np.ceil(high)
    if high <= low:
        low, high = pooled.min() - 0.5, pooled.max() + 0.5
    if heavy_tail and low > 0 and high / low > LOG_SPAN:
        return np.geomspace(low, high, N_BINS + 1)
    if integer:
        step = max(1, int(np.ceil((high - low) / N_BINS)))
        return np.arange(low, high + step, step) - 0.5
    return np.linspace(low, high, N_BINS + 1)


def is_log_bins(bins):
    widths = np.diff(bins)
    return bins[0] > 0 and not np.allclose(widths, widths[0])


def hist(values, bins):
    """Histogram with outliers folded into the first and last bins."""
    return np.histogram(np.clip(values, bins[0], bins[-1]), bins=bins)[0]


def chi2_test(values, ref_values, bins):
    """Shape-only chi-squared between two histograms with differing totals.

    Uses chi2 = sum (sqrt(N_ref/N) n - sqrt(N/N_ref) n_ref)^2 / (n + n_ref) over
    bins that are populated in at least one histogram, with ndf = bins - 1.
    Returns (chi2, ndf), or (-1, -1) if the test is not possible.
    """
    if values is None or ref_values is None or values.size == 0 or ref_values.size == 0:
        return -1.0, -1

    n = hist(values, bins)
    n_ref = hist(ref_values, bins)
    used = (n + n_ref) > 0
    ndf = int(used.sum()) - 1
    if ndf < 1:
        return -1.0, -1

    scale = np.sqrt(n_ref.sum() / n.sum())
    diff = scale * n[used] - n_ref[used] / scale
    chi2 = float(np.sum(diff**2 / (n[used] + n_ref[used])))
    return chi2, ndf


def compare(values, ref_values, bins):
    """chi2 comparison of values against one reference (-1 if unavailable)."""
    chi2, ndf = chi2_test(values, ref_values, bins)
    return {
        "chi2": chi2,
        "ndf": ndf,
        "chi2_ndf": chi2 / ndf if ndf > 0 else -1.0,
        "ref_mean": float(ref_values.mean())
        if ref_values is not None and ref_values.size
        else None,
    }


def describe(values, refs=None, bins=None):
    """Summary statistics plus chi2 against each reference.

    refs maps reference name -> values. Per-reference results are under
    "comparisons"; the first reference's results are also copied to the top
    level (chi2, ndf, chi2_ndf, ref_mean), all -1/None if there are none.
    """
    values = finite(values)
    refs = {name: finite(ref) for name, ref in (refs or {}).items()}

    if bins is None:
        bins = make_bins(values, refs.values())
    comparisons = {name: compare(values, ref, bins) for name, ref in refs.items()}
    primary = next(iter(comparisons.values()), compare(values, None, bins))

    if values.size == 0:
        base = {
            "count": 0,
            "min": None,
            "max": None,
            "mean": None,
            "std": None,
            "median": None,
        }
    else:
        base = {
            "count": int(values.size),
            "min": float(values.min()),
            "max": float(values.max()),
            "mean": float(values.mean()),
            "std": float(values.std()),
            "median": float(np.median(values)),
        }
    return {**base, **primary, "comparisons": comparisons}


def read_tree(tree):
    """Return (PFOs per event, {branch: flattened values}) from a LArRecoND tree."""
    cluster_ids = tree["clusterId"].array(library="ak")
    pfps_per_event = ak.to_numpy(ak.num(cluster_ids, axis=1))
    branch_values = {
        branch: ak.to_numpy(ak.flatten(tree[branch].array(library="ak"), axis=None))
        for branch, _ in PLOT_BRANCHES
    }
    return pfps_per_event, branch_values


def parse_reference(spec):
    """Split "name=path" (or a bare path, named after its file stem)."""
    if "=" in spec:
        name, path = spec.split("=", 1)
        return name.strip(), path.strip()
    return Path(spec).stem, spec


def populate_stats(file_path):
    """Read a reference file; returns None (with a message) if unusable."""
    if file_path is None:
        return None

    file_path = Path(file_path).expanduser()

    if not file_path.is_file():
        print(f"File is missing: {file_path}", file=sys.stderr)
        return None

    if file_path.stat().st_size == 0:
        print(f"File is empty: {file_path}", file=sys.stderr)
        return None

    with uproot.open(file_path) as root_file:
        if "LArRecoND" not in root_file:
            print(f"Missing LArRecoND tree in reference file: {file_path}", file=sys.stderr)
            return None

        tree = root_file["LArRecoND"]
        missing = {"clusterId", *(b for b, _ in PLOT_BRANCHES)} - set(tree.keys())
        if missing:
            print(f"Reference file is missing branches: {', '.join(sorted(missing))}", file=sys.stderr)
            return None

        return read_tree(tree)


def load_references(reference_specs):
    """Load every usable reference as {name: (pfps_per_event, branch_values)}."""
    references = {}
    for spec in reference_specs or []:
        name, path = parse_reference(spec)
        base, suffix = name, 2
        while name in references:
            name, suffix = f"{base}_{suffix}", suffix + 1
        loaded = populate_stats(path)
        if loaded is not None:
            references[name] = loaded
    return references


def draw_panel(figure, spec, title, xlabel, ylabel, values, refs, bins, comparisons):
    """Draw one distribution; with references, add normalised overlays + ratio."""
    refs = {name: ref for name, ref in refs.items() if ref.size}
    has_ref = bool(refs)
    if has_ref:
        inner = spec.subgridspec(2, 1, height_ratios=(3, 1), hspace=0.05)
        axis = figure.add_subplot(inner[0])
        ratio_axis = figure.add_subplot(inner[1], sharex=axis)
        axis.tick_params(labelbottom=False)
    else:
        axis = figure.add_subplot(spec)
        ratio_axis = None

    axis.set_title(title)

    if values.size == 0:
        axis.text(0.5, 0.5, "No values", ha="center", va="center", transform=axis.transAxes)
        return

    # Compare shapes (unit area) when references are present, as the samples
    # generally contain different numbers of PFOs.
    n = hist(values, bins)
    clipped = np.mean((values < bins[0]) | (values > bins[-1]))
    log_x = is_log_bins(bins)
    cur_norm = n.sum() if has_ref else 1.0
    cur = n / cur_norm
    centres = 0.5 * (bins[1:] + bins[:-1])
    widths = np.diff(bins)

    if not has_ref:
        axis.stairs(cur, bins, fill=True, color=CURRENT_COLOUR, alpha=0.6)
        axis.stairs(cur, bins, color=CURRENT_COLOUR, linewidth=1.2)
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
    else:
        ratio_axis.axhline(1.0, color="black", linewidth=0.8)
        for index, (name, ref_values) in enumerate(refs.items()):
            colour = REFERENCE_COLOURS[index % len(REFERENCE_COLOURS)]
            style = REFERENCE_STYLES[index % len(REFERENCE_STYLES)]
            n_ref = hist(ref_values, bins)
            ref = n_ref / n_ref.sum()
            if index == 0:
                axis.stairs(ref, bins, fill=True, color=colour, alpha=0.25)
            axis.stairs(ref, bins, color=colour, linewidth=1.4, linestyle=style, label=name)

            # Ratio of current to this reference, offset a little per
            # reference so overlapping points stay readable.
            ok = (n_ref >= MIN_RATIO_COUNTS) & (n >= MIN_RATIO_COUNTS)
            ratio = np.full(len(cur), np.nan)
            err = np.full(len(cur), np.nan)
            ratio[ok] = cur[ok] / ref[ok]
            err[ok] = np.sqrt(n[ok]) / cur_norm / ref[ok]
            offset = (index - (len(refs) - 1) / 2) * 0.25 * widths
            ratio_axis.errorbar(
                centres + offset, ratio, yerr=err, fmt="o", color=colour,
                markersize=3, elinewidth=1,
            )

        # Poisson errors on the current file, propagated through the normalisation.
        axis.errorbar(
            centres, cur, yerr=np.sqrt(n) / cur_norm, fmt="none",
            ecolor=CURRENT_COLOUR, elinewidth=1, capsize=0,
        )
        axis.stairs(cur, bins, color=CURRENT_COLOUR, linewidth=1.8, label="Current")

        # Ticks stop short of the limits so labels don't collide where the
        # main and ratio axes meet.
        ratio_axis.set_ylim(0, 2)
        ratio_axis.set_yticks([0.5, 1.0, 1.5])
        ratio_axis.set_ylabel("Ratio")
        ratio_axis.set_xlabel(xlabel)
        axis.set_ylabel("Fraction (unit area)")
        axis.legend(loc="upper right", fontsize=8)

        lines = []
        for name in refs:
            comparison = comparisons[name]
            if comparison["ndf"] > 0:
                lines.append(
                    f"$\\chi^2$/ndf vs {name}: {comparison['chi2']:.1f}/{comparison['ndf']}"
                    f" = {comparison['chi2_ndf']:.2f}"
                )
            else:
                lines.append(f"$\\chi^2$ vs {name}: n/a")
        axis.text(
            0.03, 0.95, "\n".join(lines), transform=axis.transAxes, va="top", fontsize=8,
            bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8, "edgecolor": "0.7"},
        )

    axis.margins(x=0)
    axis.grid(alpha=0.25)
    axis.yaxis.set_major_locator(mticker.MaxNLocator(nbins=6, prune="lower"))
    if has_ref:
        # Leave headroom for the legend and chi2 box.
        top = max(cur.max(), *(hist(r, bins).max() / r.size for r in refs.values()))
        axis.set_ylim(0, top * (1.2 + 0.1 * len(refs)))

    bottom_axis = ratio_axis if has_ref else axis
    if log_x:
        bottom_axis.set_xscale("log")
        bottom_axis.xaxis.set_major_formatter(mticker.EngFormatter())
    else:
        bottom_axis.xaxis.set_major_locator(mticker.MaxNLocator(nbins=6))
        if max(abs(bins[0]), abs(bins[-1])) >= 1000:
            bottom_axis.xaxis.set_major_formatter(mticker.EngFormatter())
    if clipped > 0:
        bottom_axis.set_xlabel(f"{xlabel}  (outliers, {clipped:.1%}, in edge bins)")
    if ratio_axis is not None:
        ratio_axis.grid(alpha=0.25)


def validate_and_plot(input_path, dataset_id, dataset_name, output_dir, reference_specs):
    input_path = Path(input_path).expanduser()
    output_dir = Path(output_dir).expanduser()

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

        pfps_per_event, branch_values = read_tree(tree)

    pfp_total = int(np.sum(pfps_per_event))
    if pfp_total == 0:
        raise ValueError(
            f"LArRecoND tree has entries but no reconstructed PFOs: {input_path}"
        )

    # Read in any reference files so the distributions can be compared.
    references = load_references(reference_specs)

    # (key, title, xlabel, ylabel, current values, {ref name: values}, integer bins,
    #  heavy tailed)
    panels = [
        (
            "pfps_per_event", "Reconstructed PFOs per event", "PFO count", "Events",
            finite(pfps_per_event),
            {name: finite(ref[0]) for name, ref in references.items()}, True, False,
        )
    ]
    for branch, title in PLOT_BRANCHES:
        panels.append(
            (
                branch, title, branch, "PFOs", finite(branch_values[branch]),
                {name: finite(ref[1][branch]) for name, ref in references.items()},
                branch == "n3DHits",
                branch in HEAVY_TAILED,
            )
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    figure = plt.figure(figsize=(15, 9), constrained_layout=False)
    grid = figure.add_gridspec(2, 3, hspace=0.35, wspace=0.25)

    stats = {}
    for index, (key, title, xlabel, ylabel, values, refs, integer, heavy) in enumerate(panels):
        bins = make_bins(values, refs.values(), integer=integer, heavy_tail=heavy)
        stats[key] = describe(values, refs, bins)
        draw_panel(
            figure, grid[index // 3, index % 3], title, xlabel, ylabel, values,
            refs, bins, stats[key]["comparisons"],
        )

    suffix = " (vs " + ", ".join(references) + ")" if references else ""
    figure.suptitle(
        f"{dataset_name}: {entries} events, {pfp_total:,} PFOs{suffix}", fontsize=14
    )
    figure.subplots_adjust(top=0.9, bottom=0.07, left=0.06, right=0.98)
    plot_path = output_dir / f"{dataset_id}.png"
    figure.savefig(plot_path, dpi=150)
    plt.close(figure)

    event_stats = stats["pfps_per_event"]
    branch_stats = [{"branch": branch, **stats[branch]} for branch, _ in PLOT_BRANCHES]

    summary = {
        "dataset_id": dataset_id,
        "dataset_name": dataset_name,
        "entries": entries,
        "pfp_total": pfp_total,
        "references": list(references),
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
        action="append",
        metavar="[NAME=]PATH",
        help="Reference ROOT file to compare against, e.g. main=ref.root. "
        "May be given multiple times; unnamed files are named after their stem.",
        default=[],
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
