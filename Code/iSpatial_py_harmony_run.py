"""Example runner for iSpatial_py - Modified for Harmony Coordinate Averaging

This script processes both soldier and forager scRNA and spatial data,
integrates them in Harmony space, and averages Harmony coordinates
instead of gene expression values.
"""

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from iSpatial_py_harmony import run_harmony, stabilize_expr

# ---------------------------------------------------------------------------
# 1. Configure paths to your data
# ---------------------------------------------------------------------------

# Default paths (can be overridden via CLI)
DEFAULT_SPRNA = Path("/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/data/SvsF_col116/combined_s39_f11_transformed_norm_log1p.h5ad")
DEFAULT_SCRNA = Path("/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/data/scRNAseq/sf_combined_cpm_log1p.h5ad")
DEFAULT_OUT   = Path("/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/code/ispatial/output_ispatial/ispatial_harmony_coords_combined_s39_f11.h5ad")
DEFAULT_PLOT_DIR = Path("/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/code/ispatial/output_ispatial/analysis_plots_harmony_coords")


def infer_harmony_coordinates(
    spRNA,
    scRNA,
    dims=tuple(range(30)),
    k_neighbor=30,
    correct_scRNA=False,
    correct_neighbor=5,
    correct_weight_NN=0.5,
):
    """
    Infer spatial cell coordinates by averaging nearest neighbor scRNA coordinates
    in Harmony space instead of averaging gene expression.
    
    For original coordinates: Apply Harmony correction only within spRNA (using source as batch)
    For updated coordinates: Apply Harmony correction between spRNA and scRNA (using tech as batch)
    
    Parameters
    ----------
    spRNA : ad.AnnData
        Spatial RNA-seq data (both soldier and forager)
    scRNA : ad.AnnData
        Single-cell RNA-seq data (both soldier and forager)
    dims : tuple
        Principal components to use for integration
    k_neighbor : int
        Number of neighbors to use for coordinate averaging
    correct_scRNA : bool
        Whether to apply expression stabilization (smoothing) to scRNA data
    correct_neighbor : int
        Number of neighbors for scRNA smoothing
    correct_weight_NN : float
        Weight for neighbor smoothing in scRNA data
        
    Returns
    -------
    ad.AnnData
        Enhanced spatial data with both original and updated Harmony coordinates
    """
    
    # Step 0: Prepare data copies and add tech labels
    spRNA = spRNA.copy()
    scRNA = scRNA.copy()

    spRNA.obs["tech"] = "spatial"
    scRNA.obs["tech"] = "scRNA"

    # Remove barcode genes if present
    if any(g.startswith("barcode_") for g in spRNA.var_names):
        keep_mask = [not g.startswith("barcode_") for g in spRNA.var_names]
        spRNA = spRNA[:, keep_mask].copy()

    # Case-insensitive gene matching
    sp_genes_lower = {gene.lower(): gene for gene in spRNA.var_names}
    sc_genes_lower = {gene.lower(): gene for gene in scRNA.var_names}
    
    # Find intersecting genes (case-insensitive)
    intersect_lower = set(sp_genes_lower.keys()) & set(sc_genes_lower.keys())
    
    # Create mapping and gene lists
    sp_to_intersect = {}
    sc_to_intersect = {}
    intersect_genes_sp = []
    intersect_genes_sc = []
    
    for gene_lower in intersect_lower:
        sp_gene = sp_genes_lower[gene_lower]
        sc_gene = sc_genes_lower[gene_lower]
        sp_to_intersect[sp_gene] = gene_lower
        sc_to_intersect[sc_gene] = gene_lower
        intersect_genes_sp.append(sp_gene)
        intersect_genes_sc.append(sc_gene)
    
    if len(intersect_genes_sp) < 10:
        raise ValueError(f"Too few intersected genes between scRNA and spRNA: {len(intersect_genes_sp)}")
    
    print(f"Found {len(intersect_genes_sp)} intersecting genes (case-insensitive)")

    # For integration, work with intersected genes only
    spRNA_subset = spRNA[:, intersect_genes_sp].copy()
    scRNA_subset = scRNA[:, intersect_genes_sc].copy()
    
    # Rename genes to have consistent naming for concatenation
    canonical_names = [sp_to_intersect[g] for g in intersect_genes_sp]
    spRNA_subset.var_names = canonical_names
    scRNA_subset.var_names = canonical_names
    
    # ======================================================================
    # STEP 0.5: Apply scRNA smoothing if requested
    # ======================================================================
    if correct_scRNA:
        print(f"🧬 Applying scRNA expression smoothing...")
        print(f"   Smoothing parameters: neighbors={correct_neighbor}, weight={correct_weight_NN}")
        scRNA_subset = stabilize_expr(
            scRNA_subset,
            neighbor=correct_neighbor,
            npcs=len(dims),
            weight_NN=correct_weight_NN,
        )
        print(f"   ✅ scRNA smoothing completed")
    else:
        print("⚠️  Skipping scRNA smoothing (correct_scRNA=False)")
    
    # ======================================================================
    # STEP 1: Generate "ORIGINAL" Harmony coordinates (spRNA only, by source)
    # ======================================================================
    print("🔬 Computing ORIGINAL Harmony coordinates (spRNA only, correcting by source)...")
    
    # Check if source information exists
    source_col = None
    for col in spRNA.obs.columns:
        if 'source' in col.lower():
            source_col = col
            break
    
    if source_col is None:
        # Try to infer source from cell names
        spRNA.obs['source'] = 'Unknown'
        for i, cell_name in enumerate(spRNA.obs_names):
            if 's39' in cell_name.lower() or 'soldier' in cell_name.lower():
                spRNA.obs.iloc[i, spRNA.obs.columns.get_loc('source')] = 'Soldier'
            elif 'f11' in cell_name.lower() or 'forager' in cell_name.lower():
                spRNA.obs.iloc[i, spRNA.obs.columns.get_loc('source')] = 'Forager'
        source_col = 'source'
    
    # Apply Harmony to spRNA only, using source as batch key
    spRNA_only = spRNA_subset.copy()
    
    # Normalize and process
    sc.pp.normalize_total(spRNA_only, target_sum=1e6)
    sc.pp.log1p(spRNA_only)
    sc.pp.highly_variable_genes(spRNA_only, n_top_genes=len(canonical_names), subset=False)
    
    # Run Harmony on spatial data only
    n_pcs_needed = max(dims) + 1
    spRNA_harmony = run_harmony(spRNA_only, canonical_names, npc=n_pcs_needed, batch_key=source_col)
    X_harmony_original = spRNA_harmony.obsm["X_harmony"]
    
    print(f"   ✓ Original Harmony: {X_harmony_original.shape} (corrected by {source_col})")
    
    # ======================================================================
    # STEP 2: Generate "UPDATED" Harmony coordinates (spRNA + scRNA, by tech)
    # ======================================================================
    print("🔬 Computing UPDATED Harmony coordinates (spRNA + scRNA, correcting by tech)...")
    
    # Concatenate datasets for tech-based correction
    adata = sc.concat([spRNA_subset, scRNA_subset], join="inner", 
                      label="batch", keys=["spatial", "scRNA"], index_unique=None)

    # Normalize and log transform
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)

    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=len(canonical_names), subset=False)

    # Run Harmony integration (tech-based)
    adata = run_harmony(adata, canonical_names, npc=n_pcs_needed, batch_key="tech")

    # Get Harmony coordinates
    X_harmony = adata.obsm["X_harmony"][:, dims]

    # Get masks for spatial and scRNA cells
    spatial_mask = adata.obs["tech"] == "spatial"
    sc_mask = ~spatial_mask

    idx_spatial = np.where(spatial_mask)[0]
    idx_sc = np.where(sc_mask)[0]

    # Find k nearest scRNA neighbors for each spatial cell
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=k_neighbor, algorithm="auto").fit(X_harmony[idx_sc])
    distances, indices = nn.kneighbors(X_harmony[idx_spatial])

    # Average Harmony coordinates instead of gene expression
    print("Averaging Harmony coordinates from nearest neighbors...")
    # Get full Harmony coordinates (all dimensions)
    X_harmony_full = adata.obsm["X_harmony"]
    updated_harmony_coords = X_harmony_full.copy()
    
    for i_spatial, neigh_idx in enumerate(indices):
        # Get the original spatial cell index
        spatial_cell_idx = idx_spatial[i_spatial]
        
        # Get neighbor coordinates in full Harmony space
        neigh_coords_full = X_harmony_full[idx_sc[neigh_idx]]
        
        # Average neighbor coordinates
        avg_coords = neigh_coords_full.mean(axis=0)
        
        # Update the spatial cell's Harmony coordinates
        updated_harmony_coords[spatial_cell_idx] = avg_coords

    # Store updated coordinates
    adata.obsm["X_harmony_updated"] = updated_harmony_coords
    
    print(f"   ✓ Updated Harmony: {updated_harmony_coords.shape} (averaged from {k_neighbor} scRNA neighbors)")
    
    # ======================================================================
    # STEP 3: Create output AnnData with both coordinate sets
    # ======================================================================
    print("📦 Creating enhanced AnnData with both coordinate sets...")
    
    # Create output AnnData with spatial cells and both coordinate sets
    enhanced_adata = spRNA.copy()
    
    # Add ORIGINAL Harmony coordinates (spatial-only correction by source)
    enhanced_adata.obsm["X_harmony_original"] = X_harmony_original
    
    # Add UPDATED Harmony coordinates (averaged from scRNA neighbors)
    enhanced_adata.obsm["X_harmony_updated"] = updated_harmony_coords[idx_spatial]
    
    # Also store the scRNA Harmony coordinates for reference
    scRNA_harmony = X_harmony_full[idx_sc]
    enhanced_adata.uns["scRNA_harmony_coords"] = scRNA_harmony
    
    # Store metadata about the integration
    enhanced_adata.uns["integration_info"] = {
        "k_neighbors": k_neighbor,
        "dims_used": list(dims),  # Convert tuple to list for HDF5 compatibility
        "n_intersect_genes": len(intersect_genes_sp),
        "intersect_genes": intersect_genes_sp,
        "original_correction": f"spRNA-only by {source_col}",
        "updated_correction": "spRNA+scRNA by tech, then neighbor averaging",
        "scRNA_smoothing": {
            "enabled": correct_scRNA,
            "neighbors": correct_neighbor if correct_scRNA else None,
            "weight": correct_weight_NN if correct_scRNA else None
        }
    }
    
    print(f"✅ Enhanced spatial data created:")
    print(f"   • Original coordinates: spRNA-only Harmony correction by {source_col}")
    print(f"   • Updated coordinates: Averaged from {k_neighbor} scRNA neighbors in tech-corrected space")
    print(f"   • scRNA smoothing: {'✅ Applied' if correct_scRNA else '❌ Disabled'}")
    if correct_scRNA:
        print(f"     - Neighbors: {correct_neighbor}, Weight: {correct_weight_NN}")
    print(f"   • Spatial cells: {len(idx_spatial)}")
    print(f"   • scRNA reference cells: {len(idx_sc)}")
    
    return enhanced_adata


def main() -> None:
    # -----------------------------------------------------------------------
    # Parse command-line arguments
    # -----------------------------------------------------------------------
    parser = argparse.ArgumentParser(description="Run Harmony coordinate averaging on MERFISH + scRNA data")
    parser.add_argument("--spatial", type=Path, default=DEFAULT_SPRNA, help="Path to spatial (MERFISH) AnnData .h5ad file")
    parser.add_argument("--scrna", type=Path,   default=DEFAULT_SCRNA, help="Path to scRNA-seq AnnData .h5ad file")
    parser.add_argument("--out",   type=Path,   default=DEFAULT_OUT,   help="Destination .h5ad for enhanced output")
    parser.add_argument("--plotdir", type=Path, default=DEFAULT_PLOT_DIR, help="Directory to save analysis plots")
    parser.add_argument("--k", type=int, default=30, help="Number of nearest neighbors for coordinate averaging")
    parser.add_argument("--smooth_scrna", action="store_true", default=True, help="Apply scRNA expression smoothing before integration")
    parser.add_argument("--smooth_neighbors", type=int, default=5, help="Number of neighbors for scRNA smoothing")
    parser.add_argument("--smooth_weight", type=float, default=0.2, help="Weight for neighbor smoothing in scRNA data")
    args = parser.parse_args()

    SPRNA_PATH = args.spatial
    SCRNA_PATH = args.scrna
    OUT_PATH   = args.out
    PLOT_DIR   = args.plotdir
    K_NEIGHBORS = args.k

    # -----------------------------------------------------------------------
    # 1. Load AnnData objects
    # -----------------------------------------------------------------------
    if not SPRNA_PATH.exists():
        raise FileNotFoundError(f"Spatial file not found: {SPRNA_PATH}")
    if not SCRNA_PATH.exists():
        raise FileNotFoundError(f"scRNA file not found: {SCRNA_PATH}")

    print(f"Loading spatial data from: {SPRNA_PATH}")
    spRNA = sc.read_h5ad(SPRNA_PATH)
    print(f"Spatial data shape: {spRNA.shape}")
    
    print(f"Loading scRNA data from: {SCRNA_PATH}")
    scRNA = sc.read_h5ad(SCRNA_PATH)
    print(f"scRNA data shape: {scRNA.shape}")

    # -----------------------------------------------------------------------
    # 2. Run Harmony coordinate averaging
    # -----------------------------------------------------------------------
    print(f"Running Harmony coordinate averaging with k={K_NEIGHBORS} neighbors...")
    print(f"scRNA smoothing settings:")
    print(f"  • Enable smoothing: {args.smooth_scrna}")
    print(f"  • Smoothing neighbors: {args.smooth_neighbors}")
    print(f"  • Smoothing weight: {args.smooth_weight}")
    
    enhanced = infer_harmony_coordinates(
        spRNA,
        scRNA,
        k_neighbor=K_NEIGHBORS,
        correct_scRNA=args.smooth_scrna,
        correct_neighbor=args.smooth_neighbors,
        correct_weight_NN=args.smooth_weight
    )

    # -----------------------------------------------------------------------
    # 3. Save results
    # -----------------------------------------------------------------------
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    print(f"Saving enhanced data to: {OUT_PATH}")
    enhanced.write_h5ad(OUT_PATH)
    
    print("✅ Harmony coordinate averaging completed successfully!")
    print(f"📁 Results saved to: {OUT_PATH}")
    print(f"📊 Use the analysis script to cluster and visualize results")


if __name__ == "__main__":
    main()