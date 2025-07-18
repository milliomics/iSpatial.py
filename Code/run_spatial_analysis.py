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

def main():
    """
    Main analysis function - modify the file paths below for your data.
    """
    print("🧠 Starting Bee Brain Spatial Analysis Pipeline...")
    print("="*60)
    
    # =================================================================
    # MODIFY THESE PATHS FOR YOUR DATA
    # =================================================================
    
    soldier_path = "/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/code/ispatial/output_ispatial/s39_col116_enhanced_merfish_full_Soldier_k30_scRNAstabilized_weight0.5_20250714.h5ad"
    
    forager_path = "/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/code/ispatial/output_ispatial/f11_col116_enhanced_merfish_full_Forager_k30_scRNAstabilized_weight0.5_20250714.h5ad"
    
    output_path = "/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/code/ispatial/output_ispatial/f11_s39_col_116_combined_soldier_forager_k30_scRNAstabilized_weight0.5_20250714_harmony.h5ad"
    
    # IMAGE SAVE PATH - MODIFY THIS TO YOUR DESIRED DIRECTORY
    save_images_path = "/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/code/ispatial/output_ispatial/analysis_plots"
    
    # Analysis parameters
    analysis_params = {
        'apply_harmony': True,
        'n_pcs': 50,
        'harmony_theta': 3.0,
        'run_clustering': True,
        'leiden_resolution': 0.2
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
        # FINAL SUMMARY
        # =================================================================
        
        print(f"\n✅ Analysis completed successfully!")
        print(f"💾 Data saved to: {output_path}")
        print(f"🖼️  Comprehensive plot saved to: {save_images_path}")
        
        # List saved images
        saved_images = list(Path(save_images_path).glob("*.png"))
        if saved_images:
            print(f"📸 Saved comprehensive image:")
            for img in sorted(saved_images):
                print(f"  • {img.name}")
                print(f"    Size: {img.stat().st_size / (1024*1024):.1f} MB")
        
        print(f"\n📋 What you got:")
        print(f"  • 16 analysis panels in one image")
        print(f"  • UMAP views (caste, batch, clusters)")
        print(f"  • Spatial views (combined, by caste, by cluster)")
        print(f"  • Statistical summaries and cluster breakdowns")
        print(f"  • Technical details and methods")
        print(f"  • Publication-ready 300 DPI PNG format")
        
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
    else:
        print("\n❌ Analysis failed. Please check error messages above.") 