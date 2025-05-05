#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerPatch


class HandlerThinRectangle(HandlerPatch):
    """Custom legend handler for thin STOP rectangles."""
    def create_artists(self, legend, orig_handle, xdescent, ydescent,
                       width, height, fontsize, trans):
        thin_width = height * 0.16
        rect_x = xdescent + (width - thin_width) * 0.5
        rect_y = ydescent - height * 0.22
        rect = mpatches.Rectangle(
            [rect_x, rect_y], thin_width, height,
            facecolor=orig_handle.get_facecolor(),
            edgecolor=orig_handle.get_edgecolor(),
            transform=trans
        )
        return [rect]


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Plot coverage and isoform tracks with fixed intron length"
    )
    parser.add_argument("infiles",
                        help="Comma-separated TSV files with Coverage, Exons, and isoform columns")
    parser.add_argument("outfile", help="Path for output PDF")
    parser.add_argument("--intron_size", type=int, default=200,
                        help="Fixed width for intron segments")
    parser.add_argument("--reverse", default=None,
                        help="Comma-separated flags (true/false) to reverse each file")
    parser.add_argument("--aspect_ratio", type=float, default=None,
                        help="Width/height ratio (default width=12)")
    parser.add_argument("--align_files", action="store_true",
                        help="Align x-axes of all files to the longest file")
    parser.add_argument("--max_diff", type=int, default=None,
                        help="Max exon-length difference for intron adjustment "
                             "(requires --align_files)")
    parser.add_argument("--direction", choices=["forward", "reverse"],
                        default="forward", help="Alignment direction")
    return parser.parse_args()


def find_exon_runs_orig(df):
    """Identify intron runs and their original lengths."""
    runs = (df["Exons"] != df["Exons"].shift()).cumsum()
    ex_runs = []
    for rid, grp in df.groupby(runs, sort=False):
        if grp["Exons"].iat[0] == "INTRON":
            start, end = grp.index.min(), grp.index.max()
            ex_runs.append({"id": rid, "orig_len": end - start + 1})
    return ex_runs


def find_exon_blocks(df):
    """Locate exon blocks and return their start index and length."""
    runs = (df["Exons"] != df["Exons"].shift()).cumsum()
    exons = []
    for _, grp in df.groupby(runs, sort=False):
        if grp["Exons"].iat[0] == "EXON":
            exons.append((grp.index.min(), grp.shape[0]))
    return exons


def expand_exons(df0, ex_runs, iso_cols, intron_size):
    """Expand introns to a fixed size and retain exon segments."""
    runs = (df0["Exons"] != df0["Exons"].shift()).cumsum()
    rows = []
    for rid, grp in df0.groupby(runs, sort=False):
        label = grp["Exons"].iat[0]
        if label == "INTRON":
            orig = next(ex["orig_len"] for ex in ex_runs if ex["id"] == rid)
            for _ in range(intron_size):
                r = {"Coverage": 0, "Exons": "INTRON",
                     "ex_run_id": rid, "orig_len_ex": orig}
                for col in iso_cols:
                    r[col] = "INTRON"
                rows.append(r)
        else:
            for _, r0 in grp.iterrows():
                r = r0.to_dict()
                r["ex_run_id"] = None
                r["orig_len_ex"] = 0
                rows.append(r)
    return pd.DataFrame(rows)


def plot_isoform(ax, df, col, color_map, intron_size):
    """Plot isoform track (exons, introns, STOP) on the given axis."""
    runs = (df[col] != df[col].shift()).cumsum()
    segments = [
        (grp[col].iat[0], grp.index.min() + 1, grp.shape[0], grp)
        for _, grp in df.groupby(runs, sort=False)
    ]
    last_idx = len(segments) - 1

    for idx, (label, start, length, grp) in enumerate(segments):
        mid = start - 0.5 + length / 2
        if label == "INTRON":
            if idx in (0, last_idx):
                continue
            intron_rows = grp[grp["ex_run_id"].notna()]
            if intron_rows.empty:
                continue
            original_len = intron_rows["orig_len_ex"].iloc[0]
            ax.hlines(0.5, start - 0.5, start - 0.5 + length,
                      color="black", linewidth=2, zorder=2)
            ax.text(start - 0.5 + length / 2, 0.65,
                    f"{original_len:,}", ha="center", va="bottom", fontsize=12, zorder=3)

        elif label == "STOP":
            exp_s = max(1, start - 3)
            exp_e = min(len(df), start + length - 1 + 3)
            ax.broken_barh(
                [(exp_s - 0.5, exp_e - exp_s + 1)],
                (-0.1, 1.2),
                facecolors=color_map.get("STOP", "red"),
                edgecolors="black",
                zorder=5
            )

        else:
            ax.broken_barh(
                [(start - 0.5, length)],
                (0, 1),
                facecolors=color_map.get(label, "lightgray"),
                edgecolors="black",
                zorder=1
            )
            ax.text(mid, 1.02,
                    f"{length:,}", ha="center", va="bottom", fontsize=12, zorder=3)

    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def main():
    """Process input files and generate the combined coverage/isoform plot."""
    args = parse_args()
    paths = [p.strip() for p in args.infiles.split(",")]
    intron_size = args.intron_size

    if args.reverse:
        rev_flags = [f.lower() in ("1", "true", "yes", "t", "y")
                     for f in args.reverse.split(",")]
    else:
        rev_flags = [False] * len(paths)

    # Load data and identify isoform columns
    initial_data = []
    for path, rev in zip(paths, rev_flags):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Input file not found: {path}")
        df0 = pd.read_csv(path, sep="\t", encoding="utf-8")
        if not {"Coverage", "Exons"}.issubset(df0.columns):
            raise ValueError(f"{path} must contain Coverage and Exons columns")
        if rev:
            df0 = df0.iloc[::-1].reset_index(drop=True)
        df0["Exons"] = df0["Exons"].astype(str).str.strip().str.upper()
        iso_cols = [c for c in df0.columns if c not in ("Coverage", "Exons")]
        initial_data.append((path, df0, iso_cols))

    # Determine all annotation labels for color mapping
    all_labels = set()
    for _, df0, iso_cols in initial_data:
        all_labels.update(df0["Exons"].unique())
        for col in iso_cols:
            all_labels.update(df0[col].unique())

    # Build color map
    defaults = {"UTR": "lightgray", "TRANSIT": "lightgreen",
                "CDS": "lightblue", "STOP": "red"}
    cmap = plt.get_cmap("tab20")
    color_map = {k: v for k, v in defaults.items() if k in all_labels}
    idx = 0
    for label in sorted(all_labels):
        if label not in color_map and label not in ("INTRON", "STOP", "EXON"):
            color_map[label] = cmap(idx)
            idx += 1
    color_map.pop("INTRON", None)
    color_map.pop("EXON", None)

    # Expand introns and prepare for plotting
    processed = []
    for path, df0, iso_cols in initial_data:
        ex_runs = find_exon_runs_orig(df0)
        df_expanded = expand_exons(df0, ex_runs, iso_cols, intron_size)
        processed.append((path, df_expanded, iso_cols))

    # Align padding between first two files if requested
    if args.align_files and args.max_diff is not None and len(processed) >= 2:
        p1, df1, iso1 = processed[0]
        p2, df2, iso2 = processed[1]
        _, df1_orig, _ = next(item for item in initial_data if item[0] == p1)
        _, df2_orig, _ = next(item for item in initial_data if item[0] == p2)

        ex1 = find_exon_blocks(df1_orig)
        ex2 = find_exon_blocks(df2_orig)
        common = min(len(ex1), len(ex2))

        new_df1, new_df2 = df1.copy(), df2.copy()
        for i in range(common):
            diff = ex1[i][1] - ex2[i][1]
            if abs(diff) > args.max_diff:
                break
            pad = abs(diff)
            if pad == 0:
                continue
            target_df, iso_cols = (new_df2, iso2) if diff > 0 else (new_df1, iso1)
            blank = []
            for _ in range(pad):
                row = {"Coverage": 0, "Exons": "INTRON",
                       "ex_run_id": None, "orig_len_ex": 0}
                for col in iso_cols:
                    row[col] = "INTRON"
                blank.append(row)
            blank_df = pd.DataFrame(blank)

            if i == 0:
                if args.direction == "forward":
                    target_df = pd.concat([blank_df, target_df], ignore_index=True)
                else:
                    ex_blocks = find_exon_blocks(target_df)
                    if ex_blocks:
                        end = ex_blocks[0][0] + ex_blocks[0][1]
                        top, bot = target_df.iloc[:end], target_df.iloc[end:]
                        target_df = pd.concat([top, blank_df, bot], ignore_index=True)
                    else:
                        target_df = pd.concat([target_df, blank_df], ignore_index=True)
            else:
                ex_blocks = find_exon_blocks(target_df)
                if i < len(ex_blocks):
                    start = ex_blocks[i][0]
                    runs = (target_df["Exons"] != target_df["Exons"].shift()).cumsum()
                    groups = list(target_df.groupby(runs, sort=False))
                    grp_idx = next((j for j, (_, g) in enumerate(groups) if start in g.index), None)
                    if grp_idx is not None:
                        intron_grp = groups[grp_idx][1]
                        if intron_grp["Exons"].iat[0] == "INTRON":
                            s = intron_grp.index.min()
                            top, bot = target_df.iloc[:s], target_df.iloc[s:]
                            target_df = pd.concat([top, blank_df, bot], ignore_index=True)

            if diff > 0:
                new_df2 = target_df
            else:
                new_df1 = target_df

        processed[0] = (p1, new_df1, iso1)
        processed[1] = (p2, new_df2, iso2)

    maxN = max(df.shape[0] for _, df, _ in processed) if args.align_files else None

    # Compute layout for coverage and isoform tracks
    total_tracks = sum(1 + len(iso_cols) for _, _, iso_cols in processed)
    coverage_height, iso_height = 6, 1
    height_ratios = []
    for _, _, iso_cols in processed:
        height_ratios.append(coverage_height)
        height_ratios.extend([iso_height] * len(iso_cols))

    fig_height = 2 * total_tracks
    fig_width = args.aspect_ratio * fig_height if args.aspect_ratio else 12
    fig, axes = plt.subplots(
        nrows=total_tracks,
        figsize=(fig_width, fig_height),
        gridspec_kw={"height_ratios": height_ratios}
    )
    axes = axes.flatten()

    idx = 0
    for _, df, iso_cols in processed:
        # Plot coverage track
        cov = df["Coverage"].astype(int).values
        x = np.arange(0.5, len(cov) + 0.5)
        xf, yf = [0.5], [0]
        for xv, yv in zip(x, cov):
            xf.extend([xv, xv]); yf.extend([yv, yv])
        xf.append(x[-1] + 0.5); yf.append(0)
        ax = axes[idx]
        ax.fill(xf, yf, color="gray", zorder=0)
        ax.set_xlim(0.5, (maxN or len(cov)) + 0.5)
        ax.set_ylim(0, None)
        ax.yaxis.set_major_formatter(mtick.StrMethodFormatter('{x:,.0f}'))
        ax.set_ylabel("Read Coverage", fontsize=16)
        ax.tick_params(axis="y", labelsize=12)
        ax.set_xticks([])
        idx += 1

        # Plot isoform tracks
        for col in iso_cols:
            ax = axes[idx]
            plot_isoform(ax, df, col, color_map, intron_size)
            ax.set_xlim(0.5, (maxN or len(df)) + 0.5)
            ax.set_ylim(-0.1, 1.2)
            ax.set_ylabel(col, rotation=0, ha="right", va="center", fontsize=14)
            idx += 1

    fig.subplots_adjust(bottom=0.08)

    # Build and draw legend
    handles = [Line2D([0], [0], color="black", lw=2)]
    labels = ["INTRON"]
    seen = set()
    for _, df, iso_cols in processed:
        for col in iso_cols:
            for lbl in df[col].str.upper().tolist():
                if lbl not in ("INTRON", "STOP", "EXON") and lbl not in seen:
                    handles.append(Rectangle((0, 0), 1, 1,
                                             facecolor=color_map[lbl],
                                             edgecolor="black"))
                    labels.append(lbl)
                    seen.add(lbl)
    handler_map = {}
    if "STOP" in all_labels:
        stop_handle = Rectangle((0, 0), 0.05, 1,
                                 facecolor=color_map.get("STOP", "red"),
                                 edgecolor="black")
        handles.append(stop_handle)
        labels.append("STOP")
        handler_map = {stop_handle: HandlerThinRectangle()}

    fig.legend(handles, labels, loc="lower center", ncol=len(handles),
               frameon=False, handlelength=1, handleheight=1,
               fontsize=18, handler_map=handler_map)

    plt.savefig(args.outfile, bbox_inches="tight", pad_inches=0.2)


if __name__ == "__main__":
    main()
