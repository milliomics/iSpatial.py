"""
Harmony Coordinates Analysis Script
==================================

This script reads the output from iSpatial_py_harmony_run.py and performs
clustering and visualization based on the averaged Harmony coordinates.

Configure your analysis parameters in the CONFIGURATION section below and
run the script.

Features:
- Load enhanced spatial data with updated Harmony coordinates
- Perform clustering using both original and updated coordinates
- Generate comparison plots
- Create spatial visualizations
- Export cluster assignments and plots
"""

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score
import warnings
warnings.filterwarnings('ignore')

# Configure scanpy
sc.settings.verbosity = 3
sc.settings.figdir = './figures/'

# ---------------------------------------------------------------------------
# CONFIGURATION: Edit these parameters for your analysis
# ---------------------------------------------------------------------------

# Input/Output paths
INPUT_PATH = Path("/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/code/ispatial/output_ispatial/ispatial_harmony_coords_combined_s39_f11_20250722.h5ad")
OUTPUT_DIR = Path("/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/code/ispatial/output_ispatial/analysis_plots_harmony_coords_20250722")

# Clustering parameters
RESOLUTION_ORIGINAL = 4.0    # Resolution for clustering original Harmony coordinates (spRNA-only)
RESOLUTION_UPDATED = 3.28   # Resolution for clustering updated Harmony coordinates (scRNA-averaged)
RANDOM_STATE = 42            # Random state for reproducible results

# Plot settings
PLOT_DPI = 600               # DPI for saved images (higher = better quality, larger file size)
PLOT_SIZE = 4                # Point size for spatial plots
FIGURE_FACECOLOR = 'white'   # Background color for saved images

# Figure size settings
FIGURE_WIDTH_PER_SOURCE = 10  # Width per source (for multi-source plots)
FIGURE_HEIGHT = 12           # Height for multi-source plots
FIGURE_WIDTH_FALLBACK = 26   # Width for fallback plots (when only one source)
FIGURE_HEIGHT_FALLBACK = 8   # Height for fallback plots


def cluster_harmony_coordinates(adata, use_rep='X_harmony_updated', resolution=1.0):
    """
    Cluster cells based on Harmony coordinates using Leiden clustering.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with Harmony coordinates
    use_rep : str
        Which representation to use for clustering (e.g., 'X_harmony_original' or 'X_harmony_updated')
    resolution : float
        Resolution for Leiden clustering
        
    Returns
    -------
    dict
        Dictionary with clustering results and metrics
    """
    print(f"🔬 Clustering cells using coordinates from: {use_rep}")
    print(f"   Leiden clustering with resolution: {resolution}")
    
    if use_rep not in adata.obsm:
        raise ValueError(f"Representation '{use_rep}' not found in adata.obsm. Available: {list(adata.obsm.keys())}")
    
    coords = adata.obsm[use_rep]
    print(f"   Coordinate matrix shape: {coords.shape}")
    
    # Create a temporary copy for clustering
    adata_temp = adata.copy()
    adata_temp.obsm['X_harmony'] = coords
    
    # Compute neighbors and run Leiden clustering (using scanpy defaults)
    print(f"   Computing neighbors and running Leiden clustering...")
    sc.pp.neighbors(adata_temp, use_rep='X_harmony')
    sc.tl.leiden(adata_temp, resolution=resolution, key_added=f'leiden_{resolution:.1f}')
    
    leiden_clusters = adata_temp.obs[f'leiden_{resolution:.1f}'].astype(int).values
    leiden_sil_score = silhouette_score(coords, leiden_clusters)
    n_leiden_clusters = len(np.unique(leiden_clusters))
    
    print(f"   Leiden: resolution={resolution}, n_clusters={n_leiden_clusters}, silhouette={leiden_sil_score:.3f}")
    
    # Store results
    results = {
        'clusters': leiden_clusters,
        'silhouette_score': leiden_sil_score,
        'n_clusters': n_leiden_clusters,
        'resolution': resolution,
        'coordinates_used': use_rep
    }
    
    # Add clustering results to adata as categorical
    adata.obs[f'leiden_{use_rep}_res{resolution:.1f}'] = pd.Categorical(leiden_clusters)
    
    return results


def compare_clustering_methods(adata, clustering_results_orig, clustering_results_updated):
    """
    Compare clustering results between original and updated Harmony coordinates.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with clustering results
    clustering_results_orig : dict
        Clustering results for original coordinates
    clustering_results_updated : dict
        Clustering results for updated coordinates
        
    Returns
    -------
    dict
        Comparison metrics
    """
    print("📊 Comparing clustering results...")
    
    comparison = {}
    
    # Compare Leiden results
    orig_clusters = clustering_results_orig['clusters']
    updated_clusters = clustering_results_updated['clusters']
    
    # Calculate Adjusted Rand Index
    ari_leiden = adjusted_rand_score(orig_clusters, updated_clusters)
    comparison['leiden_ari'] = ari_leiden
    
    # Calculate silhouette score improvement
    orig_sil = clustering_results_orig['silhouette_score']
    updated_sil = clustering_results_updated['silhouette_score']
    sil_improvement = updated_sil - orig_sil
    comparison['silhouette_improvement'] = sil_improvement
    
    print(f"   Leiden ARI (original vs updated): {ari_leiden:.3f}")
    print(f"   Silhouette improvement: {sil_improvement:.3f} ({orig_sil:.3f} → {updated_sil:.3f})")
    print(f"   Cluster counts: {clustering_results_orig['n_clusters']} → {clustering_results_updated['n_clusters']}")
    
    return comparison


def create_visualization_plots(adata, clustering_results_orig, clustering_results_updated, comparison_metrics, output_dir):
    """
    Create comprehensive visualization plots.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with all results
    clustering_results_orig : dict
        Original coordinate clustering results
    clustering_results_updated : dict
        Updated coordinate clustering results
    comparison_metrics : dict
        Comparison metrics between methods
    output_dir : Path
        Directory to save plots
    """
    print("🎨 Creating visualization plots...")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. UMAP visualization of both coordinate sets
    print("   Creating UMAP embeddings...")
    
    # UMAP for ORIGINAL coordinates (X_harmony_original - spRNA-only correction)
    print("     Computing UMAP for X_harmony_original...")
    adata_orig = adata.copy()
    adata_orig.obsm['X_harmony'] = adata.obsm['X_harmony_original']
    sc.pp.neighbors(adata_orig, use_rep='X_harmony')
    sc.tl.umap(adata_orig, random_state=RANDOM_STATE)
    adata.obsm['X_umap_original'] = adata_orig.obsm['X_umap']
    
    # UMAP for UPDATED coordinates (X_harmony_updated - scRNA neighbor-averaged)
    print("     Computing UMAP for X_harmony_updated...")
    adata_updated = adata.copy()
    adata_updated.obsm['X_harmony'] = adata.obsm['X_harmony_updated']
    sc.pp.neighbors(adata_updated, use_rep='X_harmony')
    sc.tl.umap(adata_updated, random_state=RANDOM_STATE)
    adata.obsm['X_umap_updated'] = adata_updated.obsm['X_umap']
    
    # 2. Spatial-based clustering plots
    # Check if we have source information (soldier/forager)
    has_source_info = 'source' in adata.obs.columns or any(col for col in adata.obs.columns if 'soldier' in col.lower() or 'forager' in col.lower())
    
    if has_source_info:
        # Try to identify source column
        source_col = None
        for col in adata.obs.columns:
            if 'source' in col.lower() or 'type' in col.lower() or 'caste' in col.lower():
                source_col = col
                break
        
        if source_col is None:
            # Create source info based on cell names or other patterns
            adata.obs['source'] = 'Unknown'
            for i, cell_name in enumerate(adata.obs_names):
                if 's39' in cell_name.lower() or 'soldier' in cell_name.lower():
                    adata.obs.iloc[i, adata.obs.columns.get_loc('source')] = 'Soldier'
                elif 'f11' in cell_name.lower() or 'forager' in cell_name.lower():
                    adata.obs.iloc[i, adata.obs.columns.get_loc('source')] = 'Forager'
            source_col = 'source'
        
        # Create separate plots for each source
        sources = adata.obs[source_col].unique()
        n_sources = len([s for s in sources if s not in ['Unknown', 'nan', None]])
        
        if n_sources > 1:
            fig, axes = plt.subplots(2, n_sources, figsize=(FIGURE_WIDTH_PER_SOURCE * n_sources, FIGURE_HEIGHT))
            
            # Plot clustering results for each source
            res_orig = clustering_results_orig['resolution']
            res_updated = clustering_results_updated['resolution']
            cluster_col_orig = f'leiden_X_harmony_original_res{res_orig:.1f}'
            cluster_col_updated = f'leiden_X_harmony_updated_res{res_updated:.1f}'
            
            col_idx = 0
            for source in sources:
                if source in ['Unknown', 'nan', None]:
                    continue
                    
                source_mask = adata.obs[source_col] == source
                adata_source = adata[source_mask].copy()
                
                if len(adata_source) == 0:
                    continue
                
                # Original coordinates
                sc.pl.embedding(adata_source, basis='spatial', color=cluster_col_orig,
                               ax=axes[0, col_idx], show=False, frameon=False, 
                               legend_loc='right margin', size=PLOT_SIZE, palette='tab20')
                axes[0, col_idx].set_title(f'{source} - X_harmony_original\n(spRNA-only) Leiden (n={clustering_results_orig["n_clusters"]})')
                
                # Updated coordinates
                sc.pl.embedding(adata_source, basis='spatial', color=cluster_col_updated,
                               ax=axes[1, col_idx], show=False, frameon=False,
                               legend_loc='right margin', size=PLOT_SIZE, palette='tab20')
                axes[1, col_idx].set_title(f'{source} - X_harmony_updated\n(scRNA-averaged) Leiden (n={clustering_results_updated["n_clusters"]})')
                
                col_idx += 1
            
        else:
            # Fallback to combined plot if only one source
            fig, axes = plt.subplots(1, 2, figsize=(FIGURE_WIDTH_FALLBACK, FIGURE_HEIGHT_FALLBACK))
            
            res_orig = clustering_results_orig['resolution']
            res_updated = clustering_results_updated['resolution']
            cluster_col_orig = f'leiden_X_harmony_original_res{res_orig:.1f}'
            cluster_col_updated = f'leiden_X_harmony_updated_res{res_updated:.1f}'
            
            sc.pl.embedding(adata, basis='spatial', color=cluster_col_orig,
                           ax=axes[0], show=False, frameon=False, legend_loc='right margin', 
                           size=PLOT_SIZE, palette='tab20')
            axes[0].set_title(f'X_harmony_original - Spatial View\n(spRNA-only) Leiden (n={clustering_results_orig["n_clusters"]})')
            
            sc.pl.embedding(adata, basis='spatial', color=cluster_col_updated,
                           ax=axes[1], show=False, frameon=False, legend_loc='right margin', 
                           size=PLOT_SIZE, palette='tab20')
            axes[1].set_title(f'X_harmony_updated - Spatial View\n(scRNA-averaged) Leiden (n={clustering_results_updated["n_clusters"]})')
    else:
        # No source information available, create combined spatial plot
        fig, axes = plt.subplots(1, 2, figsize=(FIGURE_WIDTH_FALLBACK, FIGURE_HEIGHT_FALLBACK))
        
        res_orig = clustering_results_orig['resolution']
        res_updated = clustering_results_updated['resolution']
        cluster_col_orig = f'leiden_X_harmony_original_res{res_orig:.1f}'
        cluster_col_updated = f'leiden_X_harmony_updated_res{res_updated:.1f}'
        
        sc.pl.embedding(adata, basis='spatial', color=cluster_col_orig,
                       ax=axes[0], show=False, frameon=False, legend_loc='right margin', 
                       size=PLOT_SIZE, palette='tab20')
        axes[0].set_title(f'X_harmony_original - Spatial View\n(spRNA-only) Leiden (n={clustering_results_orig["n_clusters"]})')
        
        sc.pl.embedding(adata, basis='spatial', color=cluster_col_updated,
                       ax=axes[1], show=False, frameon=False, legend_loc='right margin', 
                       size=PLOT_SIZE, palette='tab20')
        axes[1].set_title(f'X_harmony_updated - Spatial View\n(scRNA-averaged) Leiden (n={clustering_results_updated["n_clusters"]})')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'spatial_clustering_by_source.png', dpi=PLOT_DPI, bbox_inches='tight', facecolor=FIGURE_FACECOLOR)
    plt.close()
    
    # Only generate the main spatial clustering plot
    
    print(f"   ✅ All plots saved to: {output_dir}")


def save_results_summary(adata, clustering_results_orig, clustering_results_updated, 
                        comparison_metrics, output_dir):
    """
    Save a comprehensive summary of the analysis results.
    """
    print("💾 Saving analysis summary...")
    
    output_dir = Path(output_dir)
    
    # Save cluster assignments
    cluster_assignments = adata.obs.copy()
    cluster_assignments.to_csv(output_dir / 'cluster_assignments.csv')
    
    # Create summary report
    summary = {
        'Analysis Summary': {
            'Total cells analyzed': adata.n_obs,
            'Total genes': adata.n_vars,
            'Integration info': adata.uns.get('integration_info', 'Not available')
        },
        'Original Harmony Clustering': {
            'Leiden resolution': clustering_results_orig['resolution'],
            'Leiden clusters found': clustering_results_orig['n_clusters'],
            'Leiden silhouette': clustering_results_orig['silhouette_score']
        },
        'Updated Harmony Clustering': {
            'Leiden resolution': clustering_results_updated['resolution'],
            'Leiden clusters found': clustering_results_updated['n_clusters'],
            'Leiden silhouette': clustering_results_updated['silhouette_score']
        },
        'Comparison Metrics': comparison_metrics
    }
    
    # Save as text file
    with open(output_dir / 'analysis_summary.txt', 'w') as f:
        for section, contents in summary.items():
            f.write(f"{section}:\n")
            f.write("=" * len(section) + ":\n")
            for key, value in contents.items():
                f.write(f"  {key}: {value}\n")
            f.write("\n")
    
    # Save updated AnnData with all results
    adata.write_h5ad(output_dir / 'enhanced_with_clustering.h5ad')
    
    print(f"   ✅ Summary saved to: {output_dir}")


def main():
    print("🔬 Starting Harmony Coordinates Analysis")
    print("=" * 50)
    
    # Print configuration
    print("📋 Configuration:")
    print(f"  • Input file: {INPUT_PATH.name}")
    print(f"  • Output directory: {OUTPUT_DIR.name}")
    print(f"  • Original coordinates resolution: {RESOLUTION_ORIGINAL}")
    print(f"  • Updated coordinates resolution: {RESOLUTION_UPDATED}")
    print(f"  • Plot DPI: {PLOT_DPI}")
    print(f"  • Plot point size: {PLOT_SIZE}")
    print(f"  • Figure size (multi-source): {FIGURE_WIDTH_PER_SOURCE}×{FIGURE_HEIGHT} per source")
    print(f"  • Figure size (fallback): {FIGURE_WIDTH_FALLBACK}×{FIGURE_HEIGHT_FALLBACK}")
    print(f"  • Random state: {RANDOM_STATE}")
    print("=" * 50)
    
    # Check input file exists
    if not INPUT_PATH.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_PATH}")
    

    
    # Load data
    print("📂 Loading enhanced spatial data...")
    adata = sc.read_h5ad(INPUT_PATH)
    print(f"   Data shape: {adata.shape}")
    print(f"   Available coordinate representations: {list(adata.obsm.keys())}")
    
    # Verify required data is present
    required_keys = ['X_harmony_original', 'X_harmony_updated']
    missing_keys = [key for key in required_keys if key not in adata.obsm]
    if missing_keys:
        raise ValueError(f"Missing required coordinate representations: {missing_keys}")
    
    print(f"   ✅ Found X_harmony_original: {adata.obsm['X_harmony_original'].shape}")
    print(f"   ✅ Found X_harmony_updated: {adata.obsm['X_harmony_updated'].shape}")
    
    # Check integration info if available
    if 'integration_info' in adata.uns:
        info = adata.uns['integration_info']
        print(f"   📋 Integration info:")
        print(f"      Original correction: {info.get('original_correction', 'N/A')}")
        print(f"      Updated correction: {info.get('updated_correction', 'N/A')}")
        print(f"      K neighbors used: {info.get('k_neighbors', 'N/A')}")
        
        # Show scRNA smoothing info if available
        if 'scRNA_smoothing' in info:
            smooth_info = info['scRNA_smoothing']
            enabled = smooth_info.get('enabled', False)
            print(f"      scRNA smoothing: {'✅ Applied' if enabled else '❌ Disabled'}")
            if enabled:
                neighbors = smooth_info.get('neighbors', 'N/A')
                weight = smooth_info.get('weight', 'N/A')
                print(f"         - Neighbors: {neighbors}, Weight: {weight}")
    
    # Cluster using original Harmony coordinates (spRNA-only correction by source)
    print("\n🔍 Clustering using ORIGINAL Harmony coordinates (X_harmony_original)...")
    clustering_results_orig = cluster_harmony_coordinates(
        adata, 
        use_rep='X_harmony_original',
        resolution=RESOLUTION_ORIGINAL
    )
    
    # Cluster using updated Harmony coordinates (scRNA neighbor-averaged)
    print("\n🔍 Clustering using UPDATED Harmony coordinates (X_harmony_updated)...")
    clustering_results_updated = cluster_harmony_coordinates(
        adata, 
        use_rep='X_harmony_updated',
        resolution=RESOLUTION_UPDATED
    )
    
    # Compare clustering results
    print("\n📊 Comparing clustering results...")
    comparison_metrics = compare_clustering_methods(
        adata, clustering_results_orig, clustering_results_updated
    )
    
    # Create visualizations
    print("\n🎨 Creating visualization plots...")
    create_visualization_plots(
        adata, clustering_results_orig, clustering_results_updated, comparison_metrics, OUTPUT_DIR
    )
    
    # Save comprehensive results
    print("\n💾 Saving analysis results...")
    save_results_summary(
        adata, clustering_results_orig, clustering_results_updated, 
        comparison_metrics, OUTPUT_DIR
    )
    
    print("\n🎉 Analysis completed successfully!")
    print(f"📁 All results saved to: {OUTPUT_DIR}")
    print(f"📊 Check the plots and summary files for detailed results")


if __name__ == "__main__":
    main() 