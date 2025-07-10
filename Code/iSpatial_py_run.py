"""Example runner for iSpatial_py

Adjust `SPRNA_PATH` and `SCRNA_PATH` to point to your own datasets.
"""

from pathlib import Path
import argparse

import scanpy as sc
import anndata as ad
import numpy as np
from scipy import sparse

from iSpatial_py import infer  # or integrate_iSpatial

# ---------------------------------------------------------------------------
# 1. Configure paths to your data
# ---------------------------------------------------------------------------

# Default paths (can be overridden via CLI)

DEFAULT_SPRNA = Path("data/whole_brain_merlin/Forager/Forager_sizenorm.h5ad")
DEFAULT_SCRNA = Path("data/scRNAseq/sf_combined_cpm_log1p.h5ad")
DEFAULT_OUT   = Path("output_ispatial/enhanced_merfish_full_Forager_k30_scRNAstabilized_weight0.5.h5ad")
DEFAULT_PLOT_DIR = Path("output_ispatial/enhanced_plots_Forager_k30_scRNAstabilized_weight0.5")


def main() -> None:
    # -----------------------------------------------------------------------
    # Parse command-line arguments
    # -----------------------------------------------------------------------
    parser = argparse.ArgumentParser(description="Run iSpatial inference on MERFISH + scRNA data")
    parser.add_argument("--spatial", type=Path, default=DEFAULT_SPRNA, help="Path to spatial (MERFISH) AnnData .h5ad file")
    parser.add_argument("--scrna", type=Path,   default=DEFAULT_SCRNA, help="Path to scRNA-seq AnnData .h5ad file")
    parser.add_argument("--out",   type=Path,   default=DEFAULT_OUT,   help="Destination .h5ad for enhanced output")
    parser.add_argument("--plotdir", type=Path, default=DEFAULT_PLOT_DIR, help="Directory to save gene spatial plots (PNG)")
    args = parser.parse_args()

    SPRNA_PATH = args.spatial
    SCRNA_PATH = args.scrna
    OUT_PATH   = args.out
    PLOT_DIR   = args.plotdir

    # -----------------------------------------------------------------------
    # 1. Load AnnData objects
    # -----------------------------------------------------------------------
    if not SPRNA_PATH.exists():
        raise FileNotFoundError(f"Spatial file not found: {SPRNA_PATH}")
    if not SCRNA_PATH.exists():
        raise FileNotFoundError(f"scRNA file not found: {SCRNA_PATH}")

    spRNA = sc.read_h5ad(SPRNA_PATH)
    scRNA = sc.read_h5ad(SCRNA_PATH)

    # -----------------------------------------------------------------------
    # 3. Run iSpatial inference
    # -----------------------------------------------------------------------
    enhanced = infer(
        spRNA,
        scRNA,
        k_neighbor=30,           # choose K
        weighted_KNN=True,       # correlation-weighted assignment
        infered_layer="enhanced",
        correct_spRNA=False,
        correct_scRNA=True,
        correct_weight_NN=0.5
    )

    # -----------------------------------------------------------------------
    # 3b. Save enhanced AnnData (MERFISH cells only)
    # -----------------------------------------------------------------------
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    enhanced.write_h5ad(OUT_PATH)

    # -----------------------------------------------------------------------
    # 4. Quick visualisation: plot an inferred gene
    # -----------------------------------------------------------------------
    gene_to_plot = "DopR2"  # change to any gene present in scRNA dataset

    # Ensure plot directory exists
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    ax = sc.pl.spatial(
        enhanced,
        color=gene_to_plot,
        layer="enhanced",
        spot_size=25,
        cmap="viridis",
        show=False,
    )
    # sc.pl.spatial returns an Axes object list when show=False; grab the current figure
    import matplotlib.pyplot as plt

    fig = plt.gcf()
    plot_path = PLOT_DIR / f"{gene_to_plot}_enhanced.png"
    fig.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved spatial plot to {plot_path}")

    # Also add the enhanced expression as a layer for backward compatibility
    enhanced.layers["enhanced"] = enhanced.X.copy()


if __name__ == "__main__":
    main()