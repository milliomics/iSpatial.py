#!/usr/bin/env python3
"""
Run Spatial Analysis for Bee Brain Data
=======================================

This script runs the complete spatial analysis pipeline for combining
Soldier and Forager bee brain datasets with Harmony batch correction,
Leiden clustering, and comprehensive visualization.

Usage:
    python run_spatial_analysis.py

Author: Enhanced iSpatial Analysis Pipeline
Date: 2025
"""

import sys
import os
from pathlib import Path

# Add the current directory to path to import our tools
current_dir = Path(__file__).parent
sys.path.insert(0, str(current_dir))

# Import our spatial analysis tools
from spatial_analysis_tools import (
    combine_soldier_forager_with_clustering,
    plot_all_in_one_comprehensive,
    analyze_combined_data,
    add_leiden_clustering
)

# ---------------------------------------------------------------------------
# CONFIGURATION: Edit these parameters for your analysis
# ---------------------------------------------------------------------------

# Plot settings (matching harmony_coords_analysis.py)
PLOT_DPI = 600               # DPI for saved images (higher = better quality, larger file size)
PLOT_SIZE = 6                # Point size for spatial plots
FIGURE_FACECOLOR = 'white'   # Background color for saved images
COLOR_PALETTE = 'tab20'      # Color palette for clusters

# Figure size settings
FIGURE_WIDTH_FORAGER_SOLDIER = 22  # Width for forager vs soldier spatial plot
FIGURE_HEIGHT_FORAGER_SOLDIER = 8  # Height for forager vs soldier spatial plot

def plot_forager_soldier_spatial_simple(adata, cluster_key='clusters', save_path=None):
    """
    Create a simple side-by-side spatial plot of Forager and Soldier colored by clusters.
    Uses configuration parameters from the top of the file for consistent styling.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Combined dataset with both Forager and Soldier data
    cluster_key : str, optional
        Key for clustering results (default: 'clusters')
    save_path : str, optional
        Directory path to save plot. If None, don't save.
    """
    import matplotlib.pyplot as plt
    import scanpy as sc
    from pathlib import Path
    
    # Create save directory if specified
    if save_path:
        save_dir = Path(save_path)
        save_dir.mkdir(parents=True, exist_ok=True)
    
    # Set figure parameters
    sc.set_figure_params(facecolor=FIGURE_FACECOLOR)
    
    # Create side-by-side subplots
    fig, axes = plt.subplots(1, 2, figsize=(FIGURE_WIDTH_FORAGER_SOLDIER, FIGURE_HEIGHT_FORAGER_SOLDIER))
    
    # Separate by caste
    adata_forager = adata[adata.obs["source"] == "Forager"]
    adata_soldier = adata[adata.obs["source"] == "Soldier"]
    
    # Plot Forager spatial
    if adata_forager.n_obs > 0:
        sc.pl.embedding(adata_forager, basis="spatial", color=cluster_key,
                       title=f'Forager Spatial Clusters ({adata_forager.n_obs:,} cells)',
                       ax=axes[0], show=False, frameon=False, size=PLOT_SIZE, palette=COLOR_PALETTE)
    else:
        axes[0].text(0.5, 0.5, 'No Forager cells', ha='center', va='center',
                    transform=axes[0].transAxes, fontsize=16)
        axes[0].set_title("Forager Spatial Clusters")
    
    # Plot Soldier spatial
    if adata_soldier.n_obs > 0:
        sc.pl.embedding(adata_soldier, basis="spatial", color=cluster_key,
                       title=f'Soldier Spatial Clusters ({adata_soldier.n_obs:,} cells)',
                       ax=axes[1], show=False, frameon=False, size=PLOT_SIZE, palette=COLOR_PALETTE)
    else:
        axes[1].text(0.5, 0.5, 'No Soldier cells', ha='center', va='center',
                    transform=axes[1].transAxes, fontsize=16)
        axes[1].set_title("Soldier Spatial Clusters")
    
    # Add main title
    n_clusters = len(adata.obs[cluster_key].unique()) if cluster_key in adata.obs.columns else 0
    fig.suptitle(f'Spatial Clustering: Forager vs Soldier ({n_clusters} clusters)', 
                fontsize=16, fontweight='bold', y=0.95)
    
    plt.tight_layout()
    
    # Save if path specified
    if save_path:
        filename = f"forager_soldier_spatial_clusters.png"
        plt.savefig(save_dir / filename, dpi=PLOT_DPI, bbox_inches='tight', facecolor=FIGURE_FACECOLOR)
        print(f"  ✓ Saved Forager vs Soldier spatial plot: {filename}")
        
        # Get file size for reporting
        file_path = save_dir / filename
        if file_path.exists():
            file_size_mb = file_path.stat().st_size / (1024*1024)
            print(f"    File size: {file_size_mb:.1f} MB")
    
    plt.show()
    
    # Print summary
    print(f"\n📋 Forager vs Soldier Spatial Summary:")
    if adata_forager.n_obs > 0:
        print(f"  • Forager: {adata_forager.n_obs:,} cells")
    if adata_soldier.n_obs > 0:
        print(f"  • Soldier: {adata_soldier.n_obs:,} cells")
    if cluster_key in adata.obs.columns:
        print(f"  • Clusters displayed: {n_clusters}")
    print(f"  • High resolution: {PLOT_DPI} DPI for publication quality")
    print(f"  • Plot configuration: {PLOT_SIZE}px points, {COLOR_PALETTE} palette")

def main():
    """
    Main analysis function - modify the file paths below for your data.
    """
    print("🧠 Starting Bee Brain Spatial Analysis Pipeline...")
    print("="*60)
    
    # =================================================================
    # MODIFY THESE PATHS FOR YOUR DATA
    # =================================================================
    
    soldier_path = "/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/code/ispatial/output_ispatial/s31/s31_ispatial.h5ad"
    
    forager_path = "/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/code/ispatial/output_ispatial/f11/f11_ispatial.h5ad"
    
    output_path = "/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/code/ispatial/output_ispatial/f11_s31/f11_s31_col116.h5ad"
    
    # IMAGE SAVE PATH - MODIFY THIS TO YOUR DESIRED DIRECTORY
    save_images_path = "/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/code/ispatial/output_ispatial/f11_s31/analysis_plots"
    
    # Analysis parameters
    analysis_params = {
        'apply_harmony': True,
        'n_pcs': 50,
        'harmony_theta': 2.0,
        'run_clustering': True,
        'leiden_resolution': 2.2
    }
    
    # =================================================================
    # RUN ANALYSIS
    # =================================================================
    
    print(f"📁 Input files:")
    print(f"  • Soldier: {Path(soldier_path).name}")
    print(f"  • Forager: {Path(forager_path).name}")
    print(f"📁 Output file: {Path(output_path).name}")
    print(f"🖼️  Images will be saved to: {save_images_path}")
    print(f"\n⚙️  Analysis parameters:")
    for key, value in analysis_params.items():
        print(f"  • {key}: {value}")
    print(f"\n🎨 Plot configuration:")
    print(f"  • DPI: {PLOT_DPI}")
    print(f"  • Point size: {PLOT_SIZE}")
    print(f"  • Color palette: {COLOR_PALETTE}")
    print(f"  • Figure size: {FIGURE_WIDTH_FORAGER_SOLDIER}×{FIGURE_HEIGHT_FORAGER_SOLDIER}")
    print(f"  • Background: {FIGURE_FACECOLOR}")
    print("\n")
    
    # Create image save directory
    Path(save_images_path).mkdir(parents=True, exist_ok=True)
    
    # Check if input files exist
    if not Path(soldier_path).exists():
        print(f"❌ ERROR: Soldier file not found: {soldier_path}")
        return
    
    if not Path(forager_path).exists():
        print(f"❌ ERROR: Forager file not found: {forager_path}")
        return
    
    try:
        # Run the complete analysis pipeline
        print("🚀 Running complete analysis pipeline...")
        combined_data = combine_soldier_forager_with_clustering(
            soldier_path=soldier_path,
            forager_path=forager_path,
            output_path=output_path,
            **analysis_params
        )
        
        print("\n✅ Analysis pipeline completed successfully!")
        
        # =================================================================
        # GENERATE COMPREHENSIVE ANALYSIS SUMMARY
        # =================================================================
        
        print("\n📊 Generating analysis summary...")
        analyze_combined_data(combined_data, cluster_key='clusters', show_plots=False)
        
        # =================================================================
        # CREATE SINGLE COMPREHENSIVE VISUALIZATION
        # =================================================================
        
        print("\n🎨 Creating comprehensive visualization...")
        print("  🔬 Generating single comprehensive plot with ALL results...")
        print("  📐 Creating 4×4 grid with 16 visualization panels...")
        
        # Create one comprehensive plot that includes everything
        plot_all_in_one_comprehensive(combined_data, cluster_key='clusters', save_path=save_images_path)
        
        # =================================================================
        # CREATE SIMPLE FORAGER VS SOLDIER SPATIAL PLOT
        # =================================================================
        
        print("\n🗺️  Creating Forager vs Soldier spatial comparison...")
        print("  📐 Generating side-by-side spatial plots colored by clusters...")
        
        # Create the simple forager vs soldier spatial plot
        plot_forager_soldier_spatial_simple(combined_data, cluster_key='clusters', save_path=save_images_path)
        
        # =================================================================
        # FINAL SUMMARY
        # =================================================================
        
        print(f"\n✅ Analysis completed successfully!")
        print(f"💾 Data saved to: {output_path}")
        print(f"🖼️  Analysis plots saved to: {save_images_path}")
        
        # List saved images
        saved_images = list(Path(save_images_path).glob("*.png"))
        if saved_images:
            print(f"📸 Saved analysis images:")
            for img in sorted(saved_images):
                print(f"  • {img.name}")
                print(f"    Size: {img.stat().st_size / (1024*1024):.1f} MB")
        
        print(f"\n📋 What you got:")
        print(f"  🔬 Comprehensive Analysis:")
        print(f"    • 16 analysis panels in one comprehensive image")
        print(f"    • UMAP views (caste, batch, clusters)")
        print(f"    • Spatial views (combined, by caste, by cluster)")
        print(f"    • Statistical summaries and cluster breakdowns")
        print(f"    • Technical details and methods")
        print(f"  🗺️  Spatial Comparison:")
        print(f"    • Simple side-by-side Forager vs Soldier spatial plot")
        print(f"    • Clusters colored consistently across both castes")
        print(f"    • Clean, focused visualization for presentations")
        print(f"  📐 All images: Publication-ready {PLOT_DPI} DPI PNG format")
        
        print("\n" + "="*60)
        print("🎉 COMPREHENSIVE ANALYSIS COMPLETE!")
        print("="*60)
        
        return combined_data
        
    except Exception as e:
        print(f"\n❌ ERROR during analysis: {str(e)}")
        import traceback
        traceback.print_exc()
        return None




if __name__ == "__main__":
    """
    Run the main analysis when script is executed directly.
    """
    print("Bee Brain Spatial Transcriptomics Analysis")
    print("Enhanced iSpatial Pipeline with Harmony + Clustering")
    print("="*60)
    
    # Check if we're in the right environment
    try:
        import scanpy as sc
        import harmonypy as hm
        print("✅ Required packages available")
    except ImportError as e:
        print(f"❌ Missing required packages: {e}")
        print("Please install: pip install scanpy harmonypy matplotlib")
        sys.exit(1)
    
    # Run the main analysis
    combined_data = main()
    
    if combined_data is not None:
        print("\n📋 Data object available as 'combined_data' for further analysis")
        print("🔍 Access your results:")
        print("  • combined_data.obs['clusters'] - cluster assignments")
        print("  • combined_data.obs['source'] - caste information") 
        print("  • combined_data.obsm['X_umap'] - UMAP coordinates")
        print("  • combined_data.obsm['X_pca_harmony'] or ['X_harmony'] - batch-corrected embeddings")
        print("  • combined_data.obsm['spatial'] - spatial coordinates")
        print("\n🖼️  Generated Images:")
        print("  • complete_bee_brain_analysis_clusters.png - Comprehensive 16-panel analysis")
        print("  • forager_soldier_spatial_clusters.png - Simple Forager vs Soldier spatial comparison")
    else:
        print("\n❌ Analysis failed. Please check error messages above.") 