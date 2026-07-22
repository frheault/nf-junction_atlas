#!/usr/bin/env python
"""
One-command subject label classifier for the JUNCTION atlas.

This script parcellates white matter labels into various hierarchical classes
based on connectivity signatures described in the Junction Atlas manuscript.

Hierarchical levels:
- class-4: Broad pathway systems (Association, Projection, Commissural, Cereb).
- class-9: Simplified inter-system junction classes.
- class-31: Full intra-system overlap and bottleneck classes.
- class-lobes: Anatomical lobe parcellation (15 regions).

Usage:
  python junction_labeler.py <labels.nii.gz> [--out_dir <output_path>]
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Final, Optional

import nibabel as nib
import numpy as np
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import track

# =====================================================================
#                     LOAD MAPPING DATA
# =====================================================================

json_path = Path(__file__).parent / "data" / "junction_labels_mapping.json"
try:
    with open(json_path, "r") as f:
        mapping_data = json.load(f)
        LABEL_TO_FULL = {
            int(k): v for k, v in mapping_data["LABEL_TO_FULL"].items()
        }
        LABEL_TO_SIMPLE = {
            int(k): v for k, v in mapping_data["LABEL_TO_SIMPLE"].items()
        }
        LOBE_TO_LABELS = mapping_data["LOBE_TO_LABELS"]
        TRANSLATE_FULL = mapping_data["TRANSLATE_FULL"]
        TRANSLATE_SIMPLE = mapping_data["TRANSLATE_SIMPLE"]
except Exception as e:
    print(f"Error loading {json_path}: {e}")
    sys.exit(1)

# =====================================================================
#                             LOGGING CONFIG
# =====================================================================

console = Console()
logging.basicConfig(
    level="INFO",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, console=console)],
)
logger = logging.getLogger("junction_labeler")


# =====================================================================
#                             IMPLEMENTATION
# =====================================================================


def ensure_dir(path: Path) -> None:
    """Create directory if it does not exist."""
    path.mkdir(parents=True, exist_ok=True)


def translate_full_filename(key: str) -> str:
    """Convert a full token key like 'A-A-P-Cb' to 'Asso-Asso-Proj-Cereb'."""
    if not key or key == "nan":
        return "NA"
    parts = [t for t in key.split("-") if t]
    return "-".join(TRANSLATE_FULL.get(t, t) for t in parts)


def translate_simple_filename(letters: str) -> str:
    """Convert a simple key like 'ACP' to 'Asso-Comm-Proj'."""
    if not letters:
        return "NA"
    return "-".join(TRANSLATE_SIMPLE.get(ch, ch) for ch in letters)


def write_mask_like(
    mask_bool: np.ndarray, src_img: nib.Nifti1Image, out_path: Path
) -> None:
    """Write a mask, preserving dtype/header/affine of the source labels."""
    src_dtype = src_img.get_data_dtype()
    hdr = src_img.header.copy()
    hdr.set_data_dtype(src_dtype)
    out_img = nib.Nifti1Image(
        mask_bool.astype(src_dtype), src_img.affine, header=hdr
    )
    nib.save(out_img, str(out_path))


def process_labels(
    labels_path: Path,
    out_dir_root: Optional[Path],
    include_empty: bool,
    force: bool,
) -> None:
    """Process label image and save hierarchical masks."""
    if not labels_path.exists():
        logger.error(
            "[bold red]Error:[/bold red] The input labels file was not found: %s",
            labels_path,
        )
        sys.exit(1)

    out_root = out_dir_root or (labels_path.parent / "WM_JUNCTION")

    if out_root.exists() and any(out_root.iterdir()):
        if not force:
            logger.error(
                "[bold red]Error:[/bold red] Output directory '%s' already exists and is not empty. "
                "Use -f/--force to overwrite.",
                out_root,
            )
            sys.exit(1)

    class4_dir = out_root / "class-4"
    class9_dir = out_root / "class-9"
    class31_dir = out_root / "class-31"
    lobes_dir = out_root / "class-lobes"

    for d in (class4_dir, class9_dir, class31_dir, lobes_dir):
        ensure_dir(d)

    logger.info("Loading labels from [bold cyan]%s[/bold cyan]", labels_path)
    src_img = nib.load(str(labels_path))
    lab = np.asanyarray(src_img.dataobj)

    # ---------------- FULL (class-31 folder) ----------------
    logger.info("Processing class-31 (Full classes)...")
    full_to_labels: dict[str, list[int]] = {}
    for lbl, key in LABEL_TO_FULL.items():
        if key and key != "nan":
            full_to_labels.setdefault(key, []).append(lbl)

    for key, lbls in track(full_to_labels.items(), description="Writing masks"):
        mask = np.isin(lab, lbls)
        if not include_empty and not np.any(mask):
            continue
        name = f"{translate_full_filename(key)}.nii.gz"
        write_mask_like(mask, src_img, class31_dir / name)

    # ---------------- SIMPLIFIED (class-9 folder) ----------------
    logger.info("Processing class-9 (Simplified classes)...")
    simple_to_labels: dict[str, list[int]] = {}
    for lbl, letters in LABEL_TO_SIMPLE.items():
        if letters:
            simple_to_labels.setdefault(letters, []).append(lbl)

    for letters, lbls in track(
        simple_to_labels.items(), description="Writing masks"
    ):
        mask = np.isin(lab, lbls)
        if not include_empty and not np.any(mask):
            continue
        name = f"{translate_simple_filename(letters)}.nii.gz"
        write_mask_like(mask, src_img, class9_dir / name)

    # ---------------- BROAD (class-4 folder) ----------------
    logger.info("Processing class-4 (Broad systems)...")
    broad: dict[str, list[int]] = {"A": [], "P": [], "C": [], "B": []}
    for lbl, letters in LABEL_TO_SIMPLE.items():
        if not letters:
            continue
        for code in broad:
            if code in letters:
                broad[code].append(lbl)

    broad_names = {"A": "Asso", "P": "Proj", "C": "Comm", "B": "Cereb"}
    for code, lbls in broad.items():
        mask = np.isin(lab, lbls)
        if not include_empty and not np.any(mask):
            continue
        name = f"{broad_names[code]}.nii.gz"
        write_mask_like(mask, src_img, class4_dir / name)

    # ---------------- LOBES (class-lobes folder) ----------------
    logger.info("Processing class-lobes (Anatomical regions)...")
    for region, lbls in track(
        LOBE_TO_LABELS.items(), description="Writing masks"
    ):
        mask = np.isin(lab, lbls)
        if not include_empty and not np.any(mask):
            continue
        write_mask_like(mask, src_img, lobes_dir / f"{region}.nii.gz")

    logger.info(
        "[bold green]Success![/bold green] Results saved to: %s", out_root
    )


def main() -> None:
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Junction Atlas Subject Labeler"
    )
    parser.add_argument(
        "labels", type=Path, help="Subject labels NIfTI (.nii or .nii.gz)"
    )
    parser.add_argument(
        "--out_dir",
        type=Path,
        default=None,
        help="Output root (default: <labels_parent>/WM_JUNCTION)",
    )
    parser.add_argument(
        "--include_empty",
        action="store_true",
        help="Also write classes with no voxels in this subject",
    )

    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Overwrite existing output directories",
    )

    args = parser.parse_args()

    try:
        process_labels(
            args.labels, args.out_dir, args.include_empty, args.force
        )
    except Exception as e:
        logger.exception("An error occurred during processing: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
