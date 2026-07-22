#!/usr/bin/env python
"""Generate QA outputs for nf-junction_atlas subject results.

The figure styles intentionally mirror the manuscript visualization scripts.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import hsv_to_rgb, rgb_to_hsv
from skimage.segmentation import find_boundaries

matplotlib.use("Agg")


AXIAL_MM = [-26, -6, 8, 24, 49]
CORONAL_MM = [-60, -39, -27, -13, 14]
CLASS4_TRIPTYCH_MM = [-24, -12, 0]

FAMILY_BASE = {
    "A": np.array([0.95, 0.45, 0.15]),
    "C": np.array([0.15, 0.55, 0.95]),
    "P": np.array([0.15, 0.78, 0.35]),
    "Cb": np.array([0.67, 0.35, 0.95]),
}

CLASS4_ORDER = ["Asso", "Comm", "Proj", "Cereb"]
CLASS4_COLORS = {
    "Asso": np.array([0.30, 0.85, 0.20]),
    "Comm": np.array([0.95, 0.20, 0.20]),
    "Proj": np.array([0.18, 0.48, 0.95]),
    "Cereb": np.array([0.88, 0.38, 0.78]),
}


def first_match(root: Path, patterns: list[str]) -> Path:
    for pattern in patterns:
        matches = sorted(root.glob(pattern))
        if matches:
            return matches[0]
    raise FileNotFoundError(f"Could not find any of {patterns} under {root}")


def load_img(path: Path):
    img = nib.load(str(path))
    data = np.asanyarray(img.dataobj)
    return img, data


def normalize_t1(
    vol: np.ndarray, pmin: float = 1, pmax: float = 99.8, gamma: float = 0.9
) -> np.ndarray:
    v = vol.astype(float).copy()
    nz = v[v > 0]
    if nz.size == 0:
        return np.zeros_like(v)
    lo, hi = np.percentile(nz, [pmin, pmax])
    v = np.clip((v - lo) / (hi - lo + 1e-6), 0, 1)
    v = v**gamma
    v[vol <= 0] = 0
    return v


def index_from_world(
    affine: np.ndarray, axis: int, coord_mm: float, shape: tuple[int, ...]
) -> int:
    scale = float(affine[axis, axis])
    offset = float(affine[axis, 3])
    if abs(scale) < 1e-8:
        idx = int(round((np.linalg.inv(affine) @ np.eye(4)[axis])[axis]))
    else:
        idx = int(round((coord_mm - offset) / scale))
    return max(0, min(idx, shape[axis] - 1))


def world_coord(affine: np.ndarray, axis: int, index: int) -> float:
    return float(affine[axis, axis] * index + affine[axis, 3])


def get_slice(vol: np.ndarray, axis: int, idx: int) -> np.ndarray:
    if axis == 2:
        return vol[:, :, idx].T
    if axis == 1:
        return vol[:, idx, :].T
    if axis == 0:
        return vol[idx, :, :].T
    raise ValueError("axis must be 0, 1, or 2")


def get_bounds(mask: np.ndarray, axis: int, pad: int = 4):
    if axis == 2:
        proj = mask.any(axis=axis)
        coords = np.argwhere(proj)
        x0, y0 = coords.min(0)
        x1, y1 = coords.max(0) + 1
        return (
            max(x0 - pad, 0),
            min(x1 + pad, mask.shape[0]),
            max(y0 - pad, 0),
            min(y1 + pad, mask.shape[1]),
        )
    if axis == 1:
        proj = mask.any(axis=axis)
        coords = np.argwhere(proj)
        x0, z0 = coords.min(0)
        x1, z1 = coords.max(0) + 1
        return (
            max(x0 - pad, 0),
            min(x1 + pad, mask.shape[0]),
            max(z0 - pad, 0),
            min(z1 + pad, mask.shape[2]),
        )
    if axis == 0:
        proj = mask.any(axis=axis)
        coords = np.argwhere(proj)
        y0, z0 = coords.min(0)
        y1, z1 = coords.max(0) + 1
        return (
            max(y0 - pad, 0),
            min(y1 + pad, mask.shape[1]),
            max(z0 - pad, 0),
            min(z1 + pad, mask.shape[2]),
        )
    raise ValueError("axis must be 0, 1, or 2")


def crop_slice(slice2d: np.ndarray, bounds):
    a0, a1, b0, b1 = [int(v) for v in bounds]
    return slice2d[b0:b1, a0:a1]


def infer_type_from_row(row: dict[str, str]) -> str:
    candidate = row.get("TYPE", "")
    if candidate:
        comps = candidate.split("-")
        if all(comp in FAMILY_BASE for comp in comps):
            return candidate
    conn = row.get("CONNECTIONS", "")
    hemi = row.get("HEMI", "")
    tokens: list[str] = []
    cortical_terms = [
        "Frontal",
        "Parietal",
        "Temporal",
        "Occipital",
        "Insular",
        "Limbic",
    ]
    if any(term in conn for term in cortical_terms):
        tokens.append("C" if hemi == "C" else "A")
    if any(term in conn for term in ["Subcortical", "Brainstem"]):
        tokens.append("P")
    if "Cerebellum" in conn:
        tokens.append("Cb")
    if not tokens:
        tokens = ["A"]
    order = {"A": 0, "C": 1, "P": 2, "Cb": 3}
    return "-".join(sorted(set(tokens), key=lambda x: order[x]))


def load_label_records(labels_path: Path) -> list[dict[str, str]]:
    records = []
    with open(labels_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("LABEL"):
                row["TYPE_INFERRED"] = infer_type_from_row(row)
                records.append(row)
    return records


def build_semantic_palette(label_records: list[dict[str, str]]) -> np.ndarray:
    max_label = (
        max(int(r["LABEL"]) for r in label_records) if label_records else 0
    )
    palette = np.zeros((max_label + 1, 3), dtype=float)

    type_counts = {}
    for r in label_records:
        t = r["TYPE_INFERRED"]
        type_counts[t] = type_counts.get(t, 0) + 1

    type_seen = {}
    for row in label_records:
        t = row["TYPE_INFERRED"]
        k = type_seen.get(t, 0)
        type_seen[t] = k + 1
        gsize = type_counts[t]

        comps = t.split("-")
        base_rgb = np.mean(np.array([FAMILY_BASE[c] for c in comps]), axis=0)
        hsv = rgb_to_hsv(base_rgb)
        hemi = row.get("HEMI", "")
        if hemi in ("L", "LC"):
            hsv[0] = (hsv[0] - 0.025) % 1.0
            hsv[2] *= 0.90
        elif hemi in ("R", "RC"):
            hsv[0] = (hsv[0] + 0.025) % 1.0
            hsv[2] = min(hsv[2] * 1.03, 1.0)
        else:
            hsv[1] *= 0.75
            hsv[2] *= 0.95
        frac = k / max(gsize - 1, 1)
        hsv[1] = np.clip(hsv[1] * (0.85 + 0.35 * (frac - 0.5)), 0.45, 1.0)
        hsv[2] = np.clip(hsv[2] * (0.85 + 0.20 * (1 - frac)), 0.55, 1.0)
        palette[int(row["LABEL"])] = hsv_to_rgb(hsv)
    return palette


def make_fill_tile(
    t1_slice: np.ndarray,
    label_slice: np.ndarray,
    palette: np.ndarray,
    fill_alpha: float = 0.82,
) -> np.ndarray:
    tile = np.zeros((*t1_slice.shape, 4), dtype=float)
    brain = t1_slice > 0
    tile[brain, :3] = t1_slice[brain, None]
    tile[brain, 3] = 1.0
    m = label_slice > 0
    src = palette[label_slice[m]]
    dst = tile[m, :3]
    tile[m, :3] = src * fill_alpha + dst * (1 - fill_alpha)
    tile[m, 3] = 1.0
    return tile


def make_subpixel_outline(
    label_slice: np.ndarray, palette: np.ndarray, alpha: float = 0.72
) -> np.ndarray:
    b = find_boundaries(label_slice, mode="subpixel")
    rgba = np.zeros((*b.shape, 4), dtype=float)
    color = np.array([0.06, 0.06, 0.06])
    rgba[b, :3] = color
    rgba[b, 3] = alpha
    return rgba


def draw_semantic_tile(
    ax, x, y, w, h, t1_slice, label_slice, palette, zorder_base
):
    fill = make_fill_tile(t1_slice, label_slice, palette, fill_alpha=0.82)
    ax.imshow(
        fill,
        extent=(x, x + w, y, y + h),
        origin="lower",
        interpolation="bilinear",
        zorder=zorder_base,
    )
    outline = make_subpixel_outline(label_slice, palette, alpha=0.72)
    ax.imshow(
        outline,
        extent=(x, x + w, y, y + h),
        origin="lower",
        interpolation="nearest",
        zorder=zorder_base + 1,
    )


def make_semantic_choice_montage(
    atlas: np.ndarray,
    t1n: np.ndarray,
    atlas_img,
    palette: np.ndarray,
    out_png: Path,
):
    axial_slices = [
        index_from_world(atlas_img.affine, 2, mm, atlas.shape)
        for mm in AXIAL_MM
    ]
    coronal_slices = [
        index_from_world(atlas_img.affine, 1, mm, atlas.shape)
        for mm in CORONAL_MM
    ]

    fig = plt.figure(figsize=(14, 8), facecolor="black")
    ax = fig.add_axes([0, 0, 1, 1], facecolor="black")
    ax.axis("off")
    ax.set_aspect("equal")

    axial_bounds = get_bounds(t1n > 0, 2, pad=4)
    coronal_bounds = get_bounds(t1n > 0, 1, pad=4)

    top_tiles = []
    for idx in axial_slices:
        t1sl = crop_slice(get_slice(t1n, 2, idx), axial_bounds)
        labsl = crop_slice(get_slice(atlas, 2, idx), axial_bounds)
        top_tiles.append((idx, t1sl, labsl))

    bottom_tiles = []
    for idx in coronal_slices:
        t1sl = crop_slice(get_slice(t1n, 1, idx), coronal_bounds)
        labsl = crop_slice(get_slice(atlas, 1, idx), coronal_bounds)
        bottom_tiles.append((idx, t1sl, labsl))

    top_h = 1.0
    bot_h = 0.95
    gap_y = 0.18
    top_w = [sl.shape[1] / sl.shape[0] * top_h for _, sl, _ in top_tiles]
    bot_w = [sl.shape[1] / sl.shape[0] * bot_h for _, sl, _ in bottom_tiles]

    x_top = [0.0]
    for w in top_w[:-1]:
        x_top.append(x_top[-1] + w * 0.76)

    x_bot = [0.35]
    for w in bot_w[:-1]:
        x_bot.append(x_bot[-1] + w * 0.72)

    total_w = max(x_top[-1] + top_w[-1], x_bot[-1] + bot_w[-1]) + 0.10
    total_h = top_h + bot_h + gap_y
    ax.set_xlim(-0.05, total_w)
    ax.set_ylim(-0.05, total_h)

    for i, ((idx, t1sl, labsl), x, w) in enumerate(
        zip(top_tiles, x_top, top_w)
    ):
        y = bot_h + gap_y + (0.03 if i % 2 else 0.0)
        draw_semantic_tile(
            ax, x, y, w, top_h, t1sl, labsl, palette, zorder_base=10 + i
        )
        ax.text(
            x + 0.01 * w,
            y + 0.02,
            f"z={world_coord(atlas_img.affine, 2, idx):.0f}",
            color="white",
            fontsize=12,
            weight="bold",
            ha="left",
            va="bottom",
            zorder=100,
        )

    for i, ((idx, t1sl, labsl), x, w) in enumerate(
        zip(bottom_tiles, x_bot, bot_w)
    ):
        y = 0.02 if i % 2 else 0.0
        draw_semantic_tile(
            ax, x, y, w, bot_h, t1sl, labsl, palette, zorder_base=20 + i
        )
        ax.text(
            x + 0.01 * w,
            y + 0.02,
            f"y={world_coord(atlas_img.affine, 1, idx):.0f}",
            color="white",
            fontsize=12,
            weight="bold",
            ha="left",
            va="bottom",
            zorder=100,
        )

    fig.savefig(
        out_png,
        dpi=500,
        facecolor="black",
        bbox_inches="tight",
        pad_inches=0.02,
    )
    return fig


def build_class4_combo(
    class_masks: dict[str, np.ndarray],
) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    stack = np.stack(
        [class_masks[k].astype(np.uint8) for k in CLASS4_ORDER], axis=-1
    )
    combo = np.zeros(stack.shape[:-1], dtype=np.uint8)
    for i in range(4):
        combo |= stack[..., i] << i

    palette: dict[int, np.ndarray] = {}
    for bits in range(16):
        active = [CLASS4_ORDER[i] for i in range(4) if (bits >> i) & 1]
        if not active:
            palette[bits] = np.array([0, 0, 0], dtype=float)
        else:
            rgb = np.array(
                [CLASS4_COLORS[k] for k in active], dtype=float
            ).mean(0)
            mx = rgb.max()
            if mx > 0:
                rgb = np.clip(rgb / mx * min(mx * 1.08, 1.0), 0, 1)
            palette[bits] = rgb

    palette.update(
        {
            (1 << 0) | (1 << 1): np.array([0.90, 0.82, 0.22]),
            (1 << 0) | (1 << 2): np.array([0.20, 0.88, 0.82]),
            (1 << 1) | (1 << 2): np.array([0.90, 0.38, 0.95]),
            (1 << 1) | (1 << 2) | (1 << 0): np.array([0.95, 0.95, 0.98]),
            (1 << 2) | (1 << 3): np.array([0.72, 0.50, 0.97]),
        }
    )
    return combo, palette


def draw_class4_exact_panel(
    ax,
    t1n,
    class_masks,
    combo,
    combo_palette,
    axis,
    idx,
    fill_alpha=0.28,
    boundary_lw=1.0,
):
    bounds = get_bounds(t1n > 0, axis, pad=4)
    t1sl = crop_slice(get_slice(t1n, axis, idx), bounds)
    csl = crop_slice(get_slice(combo, axis, idx), bounds)

    ax.imshow(t1sl, cmap="gray", origin="lower", vmin=0, vmax=1)
    rgba = np.zeros((*csl.shape, 4), dtype=float)
    nz = csl > 0
    rgba[nz, :3] = np.array([combo_palette[int(v)] for v in csl[nz]])
    rgba[nz, 3] = fill_alpha
    ax.imshow(rgba, origin="lower", interpolation="nearest")

    for key in CLASS4_ORDER:
        msl = crop_slice(
            get_slice(class_masks[key].astype(np.uint8), axis, idx), bounds
        )
        if msl.any():
            ax.contour(
                msl.astype(float),
                levels=[0.5],
                colors=[CLASS4_COLORS[key]],
                linewidths=boundary_lw,
                origin="lower",
            )


def make_class4_triptych(
    t1n, label_img, hierarchy_dir: Path, out_png: Path
) -> tuple[plt.Figure, dict[str, int]]:
    class_masks = {}
    class4_counts = {}
    for key in CLASS4_ORDER:
        data = (
            np.asanyarray(
                nib.load(
                    str(hierarchy_dir / "class-4" / f"{key}.nii.gz")
                ).dataobj
            )
            > 0
        )
        class_masks[key] = data
        class4_counts[key] = int(np.count_nonzero(data))
    combo, combo_palette = build_class4_combo(class_masks)
    triptych_idxs = [
        index_from_world(label_img.affine, 1, mm, label_img.shape[:3])
        for mm in CLASS4_TRIPTYCH_MM
    ]

    fig, axs = plt.subplots(1, 3, figsize=(10.2, 3.6), facecolor="black")
    for ax, idx in zip(axs, triptych_idxs):
        ax.set_facecolor("black")
        draw_class4_exact_panel(
            ax,
            t1n,
            class_masks,
            combo,
            combo_palette,
            axis=1,
            idx=idx,
            fill_alpha=0.28,
            boundary_lw=1.0,
        )
        ax.axis("off")
    plt.subplots_adjust(
        wspace=0.02, left=0.01, right=0.99, top=0.99, bottom=0.01
    )
    fig.savefig(
        out_png,
        dpi=500,
        facecolor="black",
        bbox_inches="tight",
        pad_inches=0.01,
    )
    return fig, class4_counts


def write_label_counts(
    output_dir: Path,
    label_counts: dict[int, int],
    label_table: list[dict[str, str]],
):
    csv_path = output_dir / "junction_label_counts.csv"
    headers = [
        "LABEL",
        "VOXELS",
        "PRESENT",
        "CONNECTIONS",
        "TYPE",
        "HEMI",
        "CLASSIFICATION_NAME",
        "PAIR",
        "ATLAS_SIZE",
        "SUBJECT_TO_ATLAS_SIZE_RATIO",
    ]

    table_dict = {int(r["LABEL"]): r for r in label_table if r.get("LABEL")}

    rows = []
    for label in range(1, max(max(table_dict.keys(), default=281), 281) + 1):
        rec = table_dict.get(label, {})
        voxels = int(label_counts.get(label, 0))
        atlas_size_str = rec.get("SIZE", "")
        atlas_size = int(atlas_size_str) if atlas_size_str.isdigit() else None

        ratio = ""
        if atlas_size:
            ratio = f"{voxels / float(atlas_size):.3f}"

        rows.append(
            {
                "LABEL": label,
                "VOXELS": voxels,
                "PRESENT": bool(voxels),
                "CONNECTIONS": rec.get("CONNECTIONS", ""),
                "TYPE": rec.get("TYPE", ""),
                "HEMI": rec.get("HEMI", ""),
                "CLASSIFICATION_NAME": rec.get("CLASSIFICATION_NAME", ""),
                "PAIR": rec.get("PAIR", ""),
                "ATLAS_SIZE": atlas_size if atlas_size is not None else "",
                "SUBJECT_TO_ATLAS_SIZE_RATIO": ratio,
            }
        )

    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)

    return csv_path


def make_stats_page(label_counts, label_records, class4_counts, out_png: Path):
    table_dict = {int(r["LABEL"]): r for r in label_records if r.get("LABEL")}

    present = {k: v for k, v in label_counts.items() if v > 0}
    top = sorted(present.items(), key=lambda kv: kv[1], reverse=True)[:20]
    top_labels = []
    top_counts = []
    for label, count in top:
        name = table_dict.get(label, {}).get("CLASSIFICATION_NAME", "")
        top_labels.append(f"{label}: {name}" if name else str(label))
        top_counts.append(count)

    fig = plt.figure(figsize=(11, 8.5), dpi=180)
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.0, 1.15],
        width_ratios=[1.0, 1.0],
        hspace=0.34,
        wspace=0.34,
    )
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[1, :])

    ax0.axis("off")
    total_nonzero = sum(present.values())
    lines = [
        "Summary",
        f"Nonzero junction voxels: {total_nonzero:}",
        f"Present labels: {len(present):} / 281",
        f"Absent labels: {281 - len(present):}",
        (
            f"Top label: {top[0][0]} ({top[0][1]:,} voxels)"
            if top
            else "Top label: n/a"
        ),
    ]
    ax0.text(0.02, 0.95, lines[0], fontsize=17, fontweight="bold", va="top")
    ax0.text(
        0.02,
        0.76,
        "\n".join(lines[1:]),
        fontsize=13,
        va="top",
        linespacing=1.55,
    )

    vals = [class4_counts.get(n, 0) for n in CLASS4_ORDER]
    ax1.bar(CLASS4_ORDER, vals, color=[CLASS4_COLORS[n] for n in CLASS4_ORDER])
    ax1.set_title("Class-4 nonzero voxels", fontsize=14, fontweight="bold")
    ax1.set_ylabel("voxels")
    ax1.grid(axis="y", alpha=0.25)

    y = np.arange(len(top_labels))[::-1]
    ax2.barh(y, top_counts[::-1], color="#5b8fc9")
    ax2.set_yticks(y)
    ax2.set_yticklabels(top_labels[::-1], fontsize=8)
    ax2.set_xlabel("voxels")
    ax2.set_title(
        "Top 20 junction labels by voxel count", fontsize=14, fontweight="bold"
    )
    ax2.grid(axis="x", alpha=0.25)
    fig.suptitle(
        "Label existence and volume summary", fontsize=17, fontweight="bold"
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_png, bbox_inches="tight")
    return fig


def save_summary_text(
    path: Path,
    subject_id: str,
    labels,
    label_counts,
    class4_counts,
    hierarchy_counts,
):
    nonzero = labels[labels > 0]
    lines = [
        f"Subject: {subject_id}",
        f"Label shape: {labels.shape}",
        f"Label range: {int(labels.min())}..{int(labels.max())}",
        f"Nonzero voxels: {int(nonzero.size)}",
        f"Unique nonzero labels: {len(label_counts)}",
        f"Labels outside 1..281: {int(np.sum((nonzero < 1) | (nonzero > 281)))}",
        "",
        "Class-4 nonzero voxels:",
    ]
    for name in CLASS4_ORDER:
        lines.append(f"- {name}: {int(class4_counts.get(name, 0))}")
    lines.append("")
    lines.append("Hierarchy mask counts:")
    for name, count in hierarchy_counts.items():
        lines.append(f"- {name}: {count}")
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_dir", required=True, type=Path)
    parser.add_argument("--out_dir", type=Path, default=Path("QA"))
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Overwrite existing output directories",
    )
    args = parser.parse_args()

    in_dir = args.in_dir.resolve()
    out_dir = args.out_dir.resolve()

    if not in_dir.exists() or not in_dir.is_dir():
        print(f"Error: The input directory '{in_dir}' does not exist.")
        return 1

    out_dir.mkdir(parents=True, exist_ok=True)

    labeler_script = Path(__file__).parent / "junction_labeler.py"
    master_csv = Path(__file__).parent / "data" / "junction_labels.csv"

    if not master_csv.exists():
        print(f"Error: {master_csv} not found.")
        return 1

    for subject_dir in in_dir.iterdir():
        if not subject_dir.is_dir() or subject_dir.name in [
            "QA",
            "qa",
            "pipeline_info",
        ]:
            continue

        subject_id = subject_dir.name
        print(f"\nProcessing {subject_id}...")

        subj_out_dir = out_dir / subject_id
        if subj_out_dir.exists() and any(subj_out_dir.iterdir()):
            if not args.force:
                print(
                    f"Error: Output directory '{subj_out_dir}' already exists. Use -f/--force to overwrite."
                )
                return 1

        try:
            label_path = first_match(
                subject_dir,
                ["GENERATE_JUNCTION_SIGNATURES/*__junction_labels.nii.gz"],
            )
            t1_path = first_match(
                subject_dir, ["TRANSFORM_IMAGE_T1_MNI/*__t1_mni__warped.nii.gz"]
            )
        except FileNotFoundError as e:
            print(f"Skipping {subject_id}: {e}")
            continue

        hierarchy_dir = subject_dir / "WM_JUNCTION"

        print(f"Running labeler for {subject_id}...")
        labeler_cmd = [
            "python",
            str(labeler_script),
            str(label_path),
            "--out_dir",
            str(hierarchy_dir),
        ]
        if args.force:
            labeler_cmd.append("--force")

        subprocess.run(labeler_cmd, check=True)

        subj_out_dir = out_dir / subject_id
        fig_dir = subj_out_dir / "figures"
        table_dir = subject_dir / "LABEL_TABLE"
        subj_out_dir.mkdir(parents=True, exist_ok=True)
        fig_dir.mkdir(parents=True, exist_ok=True)
        table_dir.mkdir(parents=True, exist_ok=True)

        copied_csv = table_dir / "junction_labels.csv"
        shutil.copy2(master_csv, copied_csv)

        label_img, labels = load_img(label_path)
        _, t1 = load_img(t1_path)
        labels = labels.astype(np.int32, copy=False)
        t1n = normalize_t1(t1)

        nonzero = labels[labels > 0]
        unique, counts = np.unique(nonzero, return_counts=True)
        label_counts = {int(k): int(v) for k, v in zip(unique, counts)}

        label_records = load_label_records(copied_csv)
        csv_path = write_label_counts(table_dir, label_counts, label_records)

        semantic_png = fig_dir / f"{subject_id}_semantic_choice.png"
        class4_png = (
            fig_dir
            / f"{subject_id}_class4_triptych_yneg24_neg12_0_clean_transparent.png"
        )
        stats_png = fig_dir / f"{subject_id}_label_count_summary.png"

        palette = build_semantic_palette(label_records)
        semantic_fig = make_semantic_choice_montage(
            labels, t1n, label_img, palette, semantic_png
        )
        class4_fig, class4_counts = make_class4_triptych(
            t1n, label_img, hierarchy_dir, class4_png
        )
        stats_fig = make_stats_page(
            label_counts, label_records, class4_counts, stats_png
        )

        hierarchy_counts = {}
        for sub in ["class-4", "class-9", "class-31", "class-lobes"]:
            hierarchy_counts[sub] = len(
                list((hierarchy_dir / sub).glob("*.nii.gz"))
            )

        pdf_path = subj_out_dir / f"{subject_id}_junction_atlas_QA.pdf"
        with PdfPages(pdf_path) as pdf:
            pdf.savefig(semantic_fig, bbox_inches="tight", facecolor="black")
            pdf.savefig(class4_fig, bbox_inches="tight", facecolor="black")
            pdf.savefig(stats_fig, bbox_inches="tight")
        for fig in [semantic_fig, class4_fig, stats_fig]:
            plt.close(fig)

        summary_path = (
            subj_out_dir / f"{subject_id}_junction_atlas_QA_summary.txt"
        )
        save_summary_text(
            summary_path,
            subject_id,
            labels,
            label_counts,
            class4_counts,
            hierarchy_counts,
        )

        print(f"Completed {subject_id}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
