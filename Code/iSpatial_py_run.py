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

SPRNA_PATH = Path("/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/data/SvsF_col116/f11/f11_transformed.h5ad")
SCRNA_PATH = Path("/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/data/scRNAseq/scRNA_forager_cpm_log1p.h5ad")
OUT_PATH = Path("/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/code/ispatial/output_ispatial/f11/f11_ispatial.h5ad")



def main() -> None:
    # -----------------------------------------------------------------------
    # Parse command-line arguments
    # -----------------------------------------------------------------------
    parser = argparse.ArgumentParser(description="Run iSpatial inference on MERFISH + scRNA data")
    parser.add_argument("--spatial", type=Path, default=SPRNA_PATH, help="Path to spatial (MERFISH) AnnData .h5ad file")
    parser.add_argument("--scrna", type=Path,   default=SCRNA_PATH, help="Path to scRNA-seq AnnData .h5ad file")
    parser.add_argument("--out",   type=Path,   default=OUT_PATH,   help="Destination .h5ad for enhanced output")
    args = parser.parse_args()

    spatial_path = args.spatial
    scrna_path = args.scrna
    output_path = args.out

    # -----------------------------------------------------------------------
    # 1. Load AnnData objects
    # -----------------------------------------------------------------------
    if not spatial_path.exists():
        raise FileNotFoundError(f"Spatial file not found: {spatial_path}")
    if not scrna_path.exists():
        raise FileNotFoundError(f"scRNA file not found: {scrna_path}")

    spRNA = sc.read_h5ad(spatial_path)
    scRNA = sc.read_h5ad(scrna_path)

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
    output_path.parent.mkdir(parents=True, exist_ok=True)
    enhanced.write_h5ad(output_path)
    
    print(f"iSpatial inference completed. Enhanced data saved to {output_path}")
    print(f"Enhanced expression available in layer: 'enhanced'")


if __name__ == "__main__":
    main()