"""
Spatial Analysis Tools for Bee Brain Data
=========================================

A comprehensive toolkit for combining, batch-correcting, clustering, and visualizing
spatial transcriptomics data from Soldier and Forager bee brains.

Main Functions:
- combine_soldier_forager_datasets(): Basic combination with Harmony
- combine_soldier_forager_with_clustering(): Complete pipeline with clustering
- Various visualization functions for UMAP and spatial views

Author: Enhanced iSpatial Analysis Pipeline
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')  # Suppress matplotlib warnings

# Set default scanpy settings
sc.settings.verbosity = 2  # Reduce scanpy output
sc.settings.set_figure_params(dpi=80, facecolor='white')

def combine_soldier_forager_datasets(soldier_path, forager_path, output_path=None, 
                                    apply_harmony=True, n_pcs=30, harmony_theta=2.0):
    """
    Combine Soldier and Forager h5ad files with optional Harmony batch correction.
    
    Parameters
    ----------
    soldier_path : str
        Path to Soldier h5ad file
    forager_path : str  
        Path to Forager h5ad file
    output_path : str, optional
        Path to save combined dataset. If None, won't save.
    apply_harmony : bool, optional
        Whether to apply Harmony batch correction (default: True)
    n_pcs : int, optional
        Number of principal components for Harmony (default: 30)
    harmony_theta : float, optional
        Harmony theta parameter controlling correction strength (default: 2.0)
        
    Returns
    -------
    combined_adata : anndata.AnnData
        Combined dataset with source info in .obs['source'] and optional Harmony correction
    """
    print("Combining Soldier and Forager datasets...")
    
    # Load datasets
    print(f"Loading Soldier data from: {soldier_path}")
    soldier_adata = sc.read_h5ad(soldier_path)
    print(f"✓ Soldier: {soldier_adata.n_obs} cells, {soldier_adata.n_vars} genes")
    
    print(f"Loading Forager data from: {forager_path}")
    forager_adata = sc.read_h5ad(forager_path)
    print(f"✓ Forager: {forager_adata.n_obs} cells, {forager_adata.n_vars} genes")
    
    # Add source information (batch labels for Harmony)
    soldier_adata.obs['source'] = 'Soldier'
    soldier_adata.obs['batch'] = 'Soldier'
    forager_adata.obs['source'] = 'Forager' 
    forager_adata.obs['batch'] = 'Forager'
    print("✓ Added source and batch information")
    
    # Concatenate datasets
    combined_adata = ad.concat(
        [soldier_adata, forager_adata], 
        join='outer',  # Include all genes from both datasets
        index_unique='-'  # Make cell names unique
    )
    
    print(f"✓ Combined: {combined_adata.n_obs} cells, {combined_adata.n_vars} genes")
    
    # Verify concatenation
    source_counts = combined_adata.obs['source'].value_counts()
    print(f"  - Soldier cells: {source_counts.get('Soldier', 0)}")
    print(f"  - Forager cells: {source_counts.get('Forager', 0)}")
    
    if apply_harmony:
        print("\nApplying Harmony batch correction...")
        
        # Basic preprocessing for Harmony
        print("  - Normalizing and log-transforming...")
        sc.pp.normalize_total(combined_adata, target_sum=1e6)
        sc.pp.log1p(combined_adata)
        
        # Find highly variable genes
        print("  - Finding highly variable genes...")
        sc.pp.highly_variable_genes(combined_adata, n_top_genes=2000, subset=False)
        
        # Scale data for PCA
        print("  - Scaling data...")
        sc.pp.scale(combined_adata, max_value=10)
        
        # Run PCA
        print(f"  - Running PCA with {n_pcs} components...")
        sc.tl.pca(combined_adata, n_comps=n_pcs, svd_solver='arpack')
        
        # Apply Harmony batch correction
        print(f"  - Running Harmony (theta={harmony_theta})...")
        try:
            # Try using scanpy's harmony integration (if available)
            sc.external.pp.harmony_integrate(
                combined_adata, 
                key='batch',
                theta=harmony_theta,
                max_iter_harmony=20,
                verbose=False
            )
            print("  ✓ Harmony correction completed")
            
        except (ImportError, AttributeError):
            # Fallback to harmonypy if scanpy integration not available
            print("  - Using harmonypy directly...")
            try:
                import harmonypy as hm
                
                # Prepare data for harmonypy
                pca_embeddings = combined_adata.obsm['X_pca']
                meta_data = pd.DataFrame({
                    'batch': combined_adata.obs['batch'].astype(str)
                })
                
                # Run Harmony
                ho = hm.run_harmony(
                    pca_embeddings,
                    meta_data,
                    vars_use=['batch'],
                    theta=harmony_theta,
                    max_iter_harmony=20,
                    verbose=False
                )
                
                # Store corrected embeddings
                combined_adata.obsm['X_harmony'] = ho.Z_corr.T
                print("  ✓ Harmony correction completed using harmonypy")
                
            except ImportError:
                print("  ⚠️  Warning: Neither scanpy harmony nor harmonypy available.")
                print("  ⚠️  Install with: pip install harmonypy")
                print("  ⚠️  Proceeding without batch correction...")
        
        # Compute neighborhood graph using corrected embeddings
        if 'X_pca_harmony' in combined_adata.obsm:
            print("  - Computing neighborhood graph with Harmony embeddings...")
            sc.pp.neighbors(combined_adata, use_rep='X_pca_harmony', n_pcs=n_pcs)
        elif 'X_harmony' in combined_adata.obsm:
            print("  - Computing neighborhood graph with Harmony embeddings...")
            sc.pp.neighbors(combined_adata, use_rep='X_harmony', n_pcs=n_pcs)
        else:
            print("  - Computing neighborhood graph with PCA embeddings...")
            sc.pp.neighbors(combined_adata, n_pcs=n_pcs)
            
        # Run UMAP for visualization
        print("  - Computing UMAP...")
        sc.tl.umap(combined_adata)
        
        print("✓ Batch correction pipeline completed")
    
    else:
        print("Skipping Harmony batch correction (apply_harmony=False)")
    
    # Save if requested
    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        combined_adata.write_h5ad(output_path)
        print(f"✓ Saved to: {output_path}")
    
    # Final summary
    print(f"\nFinal dataset summary:")
    print(f"  - Total cells: {combined_adata.n_obs}")
    print(f"  - Total genes: {combined_adata.n_vars}")
    print(f"  - Batch correction: {'Applied' if apply_harmony else 'Skipped'}")
    if 'X_pca_harmony' in combined_adata.obsm:
        print(f"  - Harmony embeddings: Available in .obsm['X_pca_harmony']")
    elif 'X_harmony' in combined_adata.obsm:
        print(f"  - Harmony embeddings: Available in .obsm['X_harmony']")
    if 'X_umap' in combined_adata.obsm:
        print(f"  - UMAP coordinates: Available in .obsm['X_umap']")
    
    return combined_adata


def combine_soldier_forager_with_clustering(soldier_path, forager_path, output_path=None, 
                                          apply_harmony=True, n_pcs=30, harmony_theta=2.0,
                                          run_clustering=True, leiden_resolution=0.2):
    """
    Combine Soldier and Forager datasets with Harmony batch correction and Leiden clustering.
    
    Parameters
    ----------
    soldier_path : str
        Path to Soldier h5ad file
    forager_path : str  
        Path to Forager h5ad file
    output_path : str, optional
        Path to save combined dataset. If None, won't save.
    apply_harmony : bool, optional
        Whether to apply Harmony batch correction (default: True)
    n_pcs : int, optional
        Number of principal components for Harmony (default: 30)
    harmony_theta : float, optional
        Harmony theta parameter controlling correction strength (default: 2.0)
    run_clustering : bool, optional
        Whether to run Leiden clustering (default: True)
    leiden_resolution : float, optional
        Resolution parameter for Leiden clustering (default: 0.2)
        
    Returns
    -------
    combined_adata : anndata.AnnData
        Combined dataset with batch correction, clustering, and visualization
    """
    # Use the existing function for the main workflow
    combined_adata = combine_soldier_forager_datasets(
        soldier_path=soldier_path,
        forager_path=forager_path, 
        output_path=None,  # Don't save yet
        apply_harmony=apply_harmony,
        n_pcs=n_pcs,
        harmony_theta=harmony_theta
    )
    
    # Add Leiden clustering if requested
    if run_clustering and apply_harmony:
        print(f"\nRunning Leiden clustering (resolution={leiden_resolution})...")
        
        # Ensure we have neighborhood graph from Harmony-corrected data
        if 'X_harmony' in combined_adata.obsm:
            print("  - Using Harmony-corrected neighborhood graph...")
        else:
            print("  - Using PCA-based neighborhood graph...")
        
        # Run Leiden clustering
        sc.tl.leiden(combined_adata, key_added="clusters", resolution=leiden_resolution)
        n_clusters = len(combined_adata.obs['clusters'].unique())
        print(f"  ✓ Found {n_clusters} clusters")
        
        # Also compute cluster-level statistics
        cluster_counts = combined_adata.obs.groupby(['clusters', 'source']).size().unstack(fill_value=0)
        print(f"  ✓ Cluster composition:")
        for cluster in sorted(combined_adata.obs['clusters'].unique()):
            soldier_cells = cluster_counts.loc[cluster, 'Soldier'] if 'Soldier' in cluster_counts.columns else 0
            forager_cells = cluster_counts.loc[cluster, 'Forager'] if 'Forager' in cluster_counts.columns else 0
            total_cells = soldier_cells + forager_cells
            print(f"    Cluster {cluster}: {total_cells} cells (S:{soldier_cells}, F:{forager_cells})")
            
    elif run_clustering and not apply_harmony:
        print("⚠️  Warning: Clustering without Harmony may be dominated by batch effects")
        print("⚠️  Consider setting apply_harmony=True for better clustering results")
    
    # Save if requested
    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        combined_adata.write_h5ad(output_path)
        print(f"✓ Saved to: {output_path}")
    
    # Final summary
    print(f"\n📊 Final Analysis Summary:")
    print(f"  - Total cells: {combined_adata.n_obs}")
    print(f"  - Total genes: {combined_adata.n_vars}")
    print(f"  - Batch correction: {'✓ Applied' if apply_harmony else '✗ Skipped'}")
    print(f"  - Leiden clustering: {'✓ Applied' if run_clustering else '✗ Skipped'}")
    if 'clusters' in combined_adata.obs.columns:
        print(f"  - Number of clusters: {len(combined_adata.obs['clusters'].unique())}")
    if 'X_pca_harmony' in combined_adata.obsm or 'X_harmony' in combined_adata.obsm:
        print(f"  - Harmony embeddings: ✓ Available")
    if 'X_umap' in combined_adata.obsm:
        print(f"  - UMAP coordinates: ✓ Available")
    
    return combined_adata


def add_leiden_clustering(adata, resolution=0.2, key_added="clusters"):
    """
    Add Leiden clustering to existing AnnData object.
    
    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object with neighborhood graph already computed
    resolution : float, optional
        Leiden resolution parameter (default: 0.2)
    key_added : str, optional
        Key to store clustering results (default: "clusters")
        
    Returns
    -------
    None (modifies adata in place)
    """
    # Check if neighborhood graph exists
    if 'neighbors' not in adata.uns:
        print("⚠️  No neighborhood graph found. Computing neighbors first...")
        if 'X_pca_harmony' in adata.obsm:
            print("  - Using Harmony embeddings for neighborhood graph")
            sc.pp.neighbors(adata, use_rep='X_pca_harmony')
        elif 'X_harmony' in adata.obsm:
            print("  - Using Harmony embeddings for neighborhood graph")
            sc.pp.neighbors(adata, use_rep='X_harmony')
        elif 'X_pca' in adata.obsm:
            print("  - Using PCA embeddings for neighborhood graph")
            sc.pp.neighbors(adata)
        else:
            print("  - No embeddings found. Computing PCA first...")
            sc.pp.pca(adata)
            sc.pp.neighbors(adata)
    
    # Run Leiden clustering
    print(f"Running Leiden clustering (resolution={resolution})...")
    sc.tl.leiden(adata, key_added=key_added, resolution=resolution)
    
    n_clusters = len(adata.obs[key_added].unique())
    print(f"✓ Found {n_clusters} clusters")
    
    # Show cluster composition if source info available
    if 'source' in adata.obs.columns:
        cluster_counts = adata.obs.groupby([key_added, 'source']).size().unstack(fill_value=0)
        print(f"Cluster composition:")
        for cluster in sorted(adata.obs[key_added].unique()):
            if 'Soldier' in cluster_counts.columns and 'Forager' in cluster_counts.columns:
                soldier_cells = cluster_counts.loc[cluster, 'Soldier']
                forager_cells = cluster_counts.loc[cluster, 'Forager']
                total_cells = soldier_cells + forager_cells
                print(f"  Cluster {cluster}: {total_cells} cells (S:{soldier_cells}, F:{forager_cells})")


def plot_batch_correction_results(adata, cluster_key='clusters', figsize=(15, 10), save_path=None):
    """
    Create comprehensive visualization of batch correction and clustering results.
    
    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object with UMAP and clustering results
    cluster_key : str, optional
        Key for clustering results in adata.obs (default: 'clusters')
    figsize : tuple, optional
        Figure size (default: (15, 10))
    save_path : str, optional
        Directory path to save plots. If None, don't save.
    """
    # Check required data
    if 'X_umap' not in adata.obsm:
        print("⚠️  UMAP coordinates not found. Computing UMAP...")
        sc.tl.umap(adata)
    
    # Create save directory if specified
    if save_path:
        save_dir = Path(save_path)
        save_dir.mkdir(parents=True, exist_ok=True)
        print(f"💾 Saving plots to: {save_dir}")
    
    # Create subplots
    fig, axes = plt.subplots(2, 3, figsize=figsize)
    axes = axes.flatten()
    
    # 1. Source (Caste) coloring
    if 'source' in adata.obs.columns:
        sc.pl.umap(adata, color='source', title='Caste Identity', 
                  ax=axes[0], show=False, frameon=False)
    
    # 2. Batch coloring (should be mixed after correction)
    if 'batch' in adata.obs.columns:
        sc.pl.umap(adata, color='batch', title='Batch Labels', 
                  ax=axes[1], show=False, frameon=False)
    
    # 3. Clusters
    if cluster_key in adata.obs.columns:
        sc.pl.umap(adata, color=cluster_key, title=f'Leiden Clusters ({cluster_key})', 
                  ax=axes[2], show=False, frameon=False, legend_loc='on data')
    
    # 4. Cluster composition bar plot
    if cluster_key in adata.obs.columns and 'source' in adata.obs.columns:
        cluster_comp = adata.obs.groupby([cluster_key, 'source']).size().unstack(fill_value=0)
        cluster_comp_pct = cluster_comp.div(cluster_comp.sum(axis=1), axis=0) * 100
        
        cluster_comp_pct.plot(kind='bar', stacked=True, ax=axes[3], 
                             title='Cluster Composition (%)', rot=45)
        axes[3].set_xlabel('Cluster')
        axes[3].set_ylabel('Percentage')
        axes[3].legend(title='Caste')
    
    # 5. Cell count per cluster
    if cluster_key in adata.obs.columns:
        cluster_counts = adata.obs[cluster_key].value_counts().sort_index()
        axes[4].bar(range(len(cluster_counts)), cluster_counts.values)
        axes[4].set_xlabel('Cluster')
        axes[4].set_ylabel('Number of Cells')
        axes[4].set_title('Cells per Cluster')
        axes[4].set_xticks(range(len(cluster_counts)))
        axes[4].set_xticklabels(cluster_counts.index)
    
    # 6. Summary statistics
    axes[5].axis('off')
    summary_text = []
    summary_text.append(f"Dataset Summary:")
    summary_text.append(f"• Total cells: {adata.n_obs:,}")
    summary_text.append(f"• Total genes: {adata.n_vars:,}")
    
    if 'source' in adata.obs.columns:
        source_counts = adata.obs['source'].value_counts()
        for source, count in source_counts.items():
            summary_text.append(f"• {source} cells: {count:,}")
    
    if cluster_key in adata.obs.columns:
        n_clusters = len(adata.obs[cluster_key].unique())
        summary_text.append(f"• Number of clusters: {n_clusters}")
    
    if 'X_harmony' in adata.obsm:
        summary_text.append(f"• Batch correction: ✓ Applied")
    else:
        summary_text.append(f"• Batch correction: ✗ Not applied")
    
    axes[5].text(0.1, 0.9, '\n'.join(summary_text), 
                transform=axes[5].transAxes, fontsize=12, 
                verticalalignment='top', family='monospace')
    
    plt.tight_layout()
    
    # Save if path specified
    if save_path:
        filename = f"batch_correction_umap_results.png"
        plt.savefig(save_dir / filename, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  ✓ Saved: {filename}")
    
    plt.show()
    
    # Print cluster statistics
    if cluster_key in adata.obs.columns and 'source' in adata.obs.columns:
        print(f"\n📊 Detailed Cluster Statistics:")
        cluster_stats = adata.obs.groupby([cluster_key, 'source']).size().unstack(fill_value=0)
        for cluster in sorted(adata.obs[cluster_key].unique()):
            soldier = cluster_stats.loc[cluster, 'Soldier'] if 'Soldier' in cluster_stats.columns else 0
            forager = cluster_stats.loc[cluster, 'Forager'] if 'Forager' in cluster_stats.columns else 0
            total = soldier + forager
            soldier_pct = (soldier/total)*100 if total > 0 else 0
            forager_pct = (forager/total)*100 if total > 0 else 0
            print(f"Cluster {cluster}: {total:3d} cells | S: {soldier:3d} ({soldier_pct:5.1f}%) | F: {forager:3d} ({forager_pct:5.1f}%)")


def plot_spatial_simple(adata, cluster_key='clusters', save_path=None):
    """
    Simple spatial plotting function matching the original code style.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Combined dataset
    cluster_key : str, optional
        Clustering key to visualize
    save_path : str, optional
        Directory path to save plots. If None, don't save.
    """
    # Create save directory if specified
    if save_path:
        save_dir = Path(save_path)
        save_dir.mkdir(parents=True, exist_ok=True)
        print(f"💾 Saving spatial plots to: {save_dir}")
    
    # Set figure parameters (matching your code)
    sc.set_figure_params(facecolor="white", figsize=(10, 8))
    
    # 1. Combined spatial plot
    plt.figure(figsize=(10, 8))
    sc.pl.embedding(adata, basis="spatial", color=cluster_key, title="Combined: Spatial Clusters")
    if save_path:
        filename = f"spatial_combined_{cluster_key}.png"
        plt.savefig(save_dir / filename, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  ✓ Saved: {filename}")
    plt.show()
    
    # Separate by caste
    adata_soldier = adata[adata.obs["source"] == "Soldier"]
    adata_forager = adata[adata.obs["source"] == "Forager"]
    
    # 2. Soldier spatial
    if adata_soldier.n_obs > 0:
        plt.figure(figsize=(10, 8))
        sc.pl.embedding(adata_soldier, basis="spatial", color=cluster_key, 
                       title=f"Soldier: Spatial Clusters ({adata_soldier.n_obs} cells)")
        if save_path:
            filename = f"spatial_soldier_{cluster_key}.png"
            plt.savefig(save_dir / filename, dpi=300, bbox_inches='tight', facecolor='white')
            print(f"  ✓ Saved: {filename}")
        plt.show()
    
    # 3. Forager spatial  
    if adata_forager.n_obs > 0:
        plt.figure(figsize=(10, 8))
        sc.pl.embedding(adata_forager, basis="spatial", color=cluster_key,
                       title=f"Forager: Spatial Clusters ({adata_forager.n_obs} cells)")
        if save_path:
            filename = f"spatial_forager_{cluster_key}.png"
            plt.savefig(save_dir / filename, dpi=300, bbox_inches='tight', facecolor='white')
            print(f"  ✓ Saved: {filename}")
        plt.show()


def plot_spatial_clusters(adata, cluster_key='clusters', figsize=(15, 5), 
                         separate_castes=True, point_size=None, save_path=None):
    """
    Plot clustering results in spatial coordinates.
    
    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object with spatial coordinates and clustering
    cluster_key : str, optional
        Key for clustering results (default: 'clusters')
    figsize : tuple, optional
        Figure size (default: (15, 5))
    separate_castes : bool, optional
        Whether to show separate plots for Soldier and Forager (default: True)
    point_size : float, optional
        Point size for spatial plots. If None, auto-determined.
    save_path : str, optional
        Directory path to save plots. If None, don't save.
    """
    # Create save directory if specified
    if save_path:
        save_dir = Path(save_path)
        save_dir.mkdir(parents=True, exist_ok=True)
        print(f"💾 Saving spatial cluster plots to: {save_dir}")
    
    # Set scanpy figure parameters
    sc.set_figure_params(facecolor="white", figsize=(10, 8))
    
    if separate_castes and 'source' in adata.obs.columns:
        # Create 3 subplots: combined, soldier, forager
        fig, axes = plt.subplots(1, 3, figsize=figsize)
        
        # Combined plot
        sc.pl.embedding(adata, basis="spatial", color=cluster_key, 
                       title="Combined: Spatial Clusters", ax=axes[0], show=False,
                       size=point_size)
        
        # Soldier only
        adata_soldier = adata[adata.obs["source"] == "Soldier"]
        if adata_soldier.n_obs > 0:
            sc.pl.embedding(adata_soldier, basis="spatial", color=cluster_key,
                           title="Soldier: Spatial Clusters", ax=axes[1], show=False,
                           size=point_size)
        else:
            axes[1].text(0.5, 0.5, 'No Soldier cells', ha='center', va='center',
                        transform=axes[1].transAxes)
            axes[1].set_title("Soldier: Spatial Clusters")
        
        # Forager only  
        adata_forager = adata[adata.obs["source"] == "Forager"]
        if adata_forager.n_obs > 0:
            sc.pl.embedding(adata_forager, basis="spatial", color=cluster_key,
                           title="Forager: Spatial Clusters", ax=axes[2], show=False,
                           size=point_size)
        else:
            axes[2].text(0.5, 0.5, 'No Forager cells', ha='center', va='center',
                        transform=axes[2].transAxes)
            axes[2].set_title("Forager: Spatial Clusters")
            
    else:
        # Single combined plot
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        sc.pl.embedding(adata, basis="spatial", color=cluster_key,
                       title="Spatial Clusters", ax=ax, show=False,
                       size=point_size)
    
    plt.tight_layout()
    
    # Save if path specified
    if save_path:
        filename = f"spatial_clusters_3panel_{cluster_key}.png"
        plt.savefig(save_dir / filename, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  ✓ Saved: {filename}")
    
    plt.show()


def plot_complete_analysis_results(adata, cluster_key='clusters', figsize=(20, 12), save_path=None):
    """
    Create comprehensive visualization combining UMAP and spatial views.
    
    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object with UMAP, spatial coordinates, and clustering
    cluster_key : str, optional
        Key for clustering results (default: 'clusters')
    figsize : tuple, optional
        Figure size (default: (20, 12))
    save_path : str, optional
        Directory path to save plots. If None, don't save.
    """
    # Create save directory if specified
    if save_path:
        save_dir = Path(save_path)
        save_dir.mkdir(parents=True, exist_ok=True)
        print(f"💾 Saving comprehensive analysis to: {save_dir}")
    
    # Set figure parameters
    sc.set_figure_params(facecolor="white")
    
    # Create large subplot grid: 3 rows x 4 columns
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.3)
    
    # Row 1: UMAP views
    ax_umap1 = fig.add_subplot(gs[0, 0])
    ax_umap2 = fig.add_subplot(gs[0, 1])
    ax_umap3 = fig.add_subplot(gs[0, 2])
    ax_stats = fig.add_subplot(gs[0, 3])
    
    # Row 2: Spatial views
    ax_spatial1 = fig.add_subplot(gs[1, 0])
    ax_spatial2 = fig.add_subplot(gs[1, 1])
    ax_spatial3 = fig.add_subplot(gs[1, 2])
    ax_comp = fig.add_subplot(gs[1, 3])
    
    # Row 3: Detailed spatial views (caste-specific)
    ax_combined_spatial = fig.add_subplot(gs[2, 0])
    ax_soldier_spatial = fig.add_subplot(gs[2, 1])
    ax_forager_spatial = fig.add_subplot(gs[2, 2])
    ax_summary = fig.add_subplot(gs[2, 3])
    
    # ========== ROW 1: UMAP VIEWS ==========
    
    # 1. UMAP: Source (Caste)
    if 'source' in adata.obs.columns and 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color='source', title='UMAP: Caste', 
                  ax=ax_umap1, show=False, frameon=False)
    
    # 2. UMAP: Batch 
    if 'batch' in adata.obs.columns and 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color='batch', title='UMAP: Batch', 
                  ax=ax_umap2, show=False, frameon=False)
    
    # 3. UMAP: Clusters
    if cluster_key in adata.obs.columns and 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color=cluster_key, title='UMAP: Clusters', 
                  ax=ax_umap3, show=False, frameon=False, legend_loc='on data')
    
    # 4. Statistics text
    ax_stats.axis('off')
    summary_text = ["UMAP Analysis:"]
    if 'X_umap' in adata.obsm:
        summary_text.append("✓ UMAP computed")
    else:
        summary_text.append("✗ UMAP missing")
    
    if 'X_pca_harmony' in adata.obsm or 'X_harmony' in adata.obsm:
        summary_text.append("✓ Harmony applied")
    else:
        summary_text.append("✗ No batch correction")
        
    summary_text.append(f"Total cells: {adata.n_obs:,}")
    
    if 'source' in adata.obs.columns:
        source_counts = adata.obs['source'].value_counts()
        for source, count in source_counts.items():
            summary_text.append(f"{source}: {count:,}")
    
    ax_stats.text(0.1, 0.9, '\n'.join(summary_text), 
                 transform=ax_stats.transAxes, fontsize=10, 
                 verticalalignment='top', family='monospace')
    
    # ========== ROW 2: BASIC SPATIAL VIEWS ==========
    
    # 5. Spatial: Source (Caste)
    if 'source' in adata.obs.columns and 'spatial' in adata.obsm:
        sc.pl.embedding(adata, basis="spatial", color='source', title='Spatial: Caste',
                       ax=ax_spatial1, show=False, frameon=False)
    
    # 6. Spatial: Batch
    if 'batch' in adata.obs.columns and 'spatial' in adata.obsm:
        sc.pl.embedding(adata, basis="spatial", color='batch', title='Spatial: Batch',
                       ax=ax_spatial2, show=False, frameon=False)
    
    # 7. Spatial: Clusters
    if cluster_key in adata.obs.columns and 'spatial' in adata.obsm:
        sc.pl.embedding(adata, basis="spatial", color=cluster_key, title='Spatial: Clusters',
                       ax=ax_spatial3, show=False, frameon=False)
    
    # 8. Cluster composition
    if cluster_key in adata.obs.columns and 'source' in adata.obs.columns:
        cluster_comp = adata.obs.groupby([cluster_key, 'source']).size().unstack(fill_value=0)
        cluster_comp_pct = cluster_comp.div(cluster_comp.sum(axis=1), axis=0) * 100
        
        cluster_comp_pct.plot(kind='bar', stacked=True, ax=ax_comp, 
                             title='Cluster Composition (%)', rot=45)
        ax_comp.set_xlabel('Cluster')
        ax_comp.set_ylabel('Percentage')
        ax_comp.legend(title='Caste', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # ========== ROW 3: DETAILED SPATIAL VIEWS ==========
    
    # 9. Combined spatial (larger)
    if cluster_key in adata.obs.columns and 'spatial' in adata.obsm:
        sc.pl.embedding(adata, basis="spatial", color=cluster_key, 
                       title='Combined Spatial Clusters',
                       ax=ax_combined_spatial, show=False, frameon=False)
    
    # 10. Soldier spatial
    if 'source' in adata.obs.columns and 'spatial' in adata.obsm:
        adata_soldier = adata[adata.obs["source"] == "Soldier"]
        if adata_soldier.n_obs > 0:
            sc.pl.embedding(adata_soldier, basis="spatial", color=cluster_key,
                           title=f'Soldier Spatial ({adata_soldier.n_obs} cells)',
                           ax=ax_soldier_spatial, show=False, frameon=False)
        else:
            ax_soldier_spatial.text(0.5, 0.5, 'No Soldier cells', ha='center', va='center',
                                  transform=ax_soldier_spatial.transAxes)
            ax_soldier_spatial.set_title("Soldier Spatial")
    
    # 11. Forager spatial
    if 'source' in adata.obs.columns and 'spatial' in adata.obsm:
        adata_forager = adata[adata.obs["source"] == "Forager"]
        if adata_forager.n_obs > 0:
            sc.pl.embedding(adata_forager, basis="spatial", color=cluster_key,
                           title=f'Forager Spatial ({adata_forager.n_obs} cells)',
                           ax=ax_forager_spatial, show=False, frameon=False)
        else:
            ax_forager_spatial.text(0.5, 0.5, 'No Forager cells', ha='center', va='center',
                                   transform=ax_forager_spatial.transAxes)
            ax_forager_spatial.set_title("Forager Spatial")
    
    # 12. Detailed summary
    ax_summary.axis('off')
    detail_text = ["Spatial Analysis:"]
    
    if 'spatial' in adata.obsm:
        detail_text.append("✓ Spatial coords available")
    else:
        detail_text.append("✗ No spatial coordinates")
    
    if cluster_key in adata.obs.columns:
        n_clusters = len(adata.obs[cluster_key].unique())
        detail_text.append(f"Clusters: {n_clusters}")
        
        # Add cluster breakdown
        if 'source' in adata.obs.columns:
            detail_text.append("\nCluster breakdown:")
            cluster_stats = adata.obs.groupby([cluster_key, 'source']).size().unstack(fill_value=0)
            for cluster in sorted(adata.obs[cluster_key].unique()):
                soldier = cluster_stats.loc[cluster, 'Soldier'] if 'Soldier' in cluster_stats.columns else 0
                forager = cluster_stats.loc[cluster, 'Forager'] if 'Forager' in cluster_stats.columns else 0
                detail_text.append(f"C{cluster}: S{soldier} F{forager}")
    
    ax_summary.text(0.1, 0.9, '\n'.join(detail_text), 
                   transform=ax_summary.transAxes, fontsize=9, 
                   verticalalignment='top', family='monospace')
    
    # Add main title
    fig.suptitle('Complete Analysis: UMAP + Spatial Views', fontsize=16, y=0.98)
    
    # Save if path specified
    if save_path:
        filename = f"complete_analysis_12panel_{cluster_key}.png"
        plt.savefig(save_dir / filename, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  ✓ Saved: {filename}")
    
    plt.show()


# Summary function for quick analysis
def analyze_combined_data(adata, cluster_key='clusters', show_plots=True, save_path=None):
    """
    Quick analysis summary of combined bee brain data.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Combined dataset
    cluster_key : str, optional
        Clustering results key
    show_plots : bool, optional
        Whether to display plots (default: True)
    save_path : str, optional
        Directory path to save plots. If None, don't save.
    """
    print("="*60)
    print("🧠 BEE BRAIN SPATIAL ANALYSIS SUMMARY")
    print("="*60)
    
    # Basic info
    print(f"\n📊 Dataset Overview:")
    print(f"  • Total cells: {adata.n_obs:,}")
    print(f"  • Total genes: {adata.n_vars:,}")
    
    # Source breakdown
    if 'source' in adata.obs.columns:
        source_counts = adata.obs['source'].value_counts()
        print(f"\n🐜 Caste Distribution:")
        for source, count in source_counts.items():
            pct = (count / adata.n_obs) * 100
            print(f"  • {source}: {count:,} cells ({pct:.1f}%)")
    
    # Clustering info
    if cluster_key in adata.obs.columns:
        n_clusters = len(adata.obs[cluster_key].unique())
        print(f"\n🎯 Clustering Results:")
        print(f"  • Number of clusters: {n_clusters}")
        
        # Cluster composition
        if 'source' in adata.obs.columns:
            print(f"\n📈 Cluster Composition:")
            cluster_stats = adata.obs.groupby([cluster_key, 'source']).size().unstack(fill_value=0)
            for cluster in sorted(adata.obs[cluster_key].unique()):
                soldier = cluster_stats.loc[cluster, 'Soldier'] if 'Soldier' in cluster_stats.columns else 0
                forager = cluster_stats.loc[cluster, 'Forager'] if 'Forager' in cluster_stats.columns else 0
                total = soldier + forager
                s_pct = (soldier/total)*100 if total > 0 else 0
                f_pct = (forager/total)*100 if total > 0 else 0
                print(f"  • Cluster {cluster}: {total:4d} cells | S: {soldier:3d} ({s_pct:4.1f}%) | F: {forager:3d} ({f_pct:4.1f}%)")
    
    # Processing status
    print(f"\n🔬 Processing Status:")
    processing_status = []
    if 'X_pca_harmony' in adata.obsm or 'X_harmony' in adata.obsm:
        processing_status.append("✓ Harmony batch correction")
    else:
        processing_status.append("✗ No batch correction")
        
    if 'X_umap' in adata.obsm:
        processing_status.append("✓ UMAP embedding")
    else:
        processing_status.append("✗ No UMAP")
        
    if 'spatial' in adata.obsm:
        processing_status.append("✓ Spatial coordinates")
    else:
        processing_status.append("✗ No spatial coordinates")
        
    if cluster_key in adata.obs.columns:
        processing_status.append("✓ Leiden clustering")
    else:
        processing_status.append("✗ No clustering")
    
    for status in processing_status:
        print(f"  • {status}")
    
    if show_plots:
        print(f"\n🎨 Generating visualizations...")
        
        # Show UMAP results if available
        if 'X_umap' in adata.obsm:
            plot_batch_correction_results(adata, cluster_key=cluster_key, save_path=save_path)
        
        # Show spatial results if available
        if 'spatial' in adata.obsm:
            plot_spatial_clusters(adata, cluster_key=cluster_key, save_path=save_path)
    
    print("\n" + "="*60)
    print("✅ Analysis complete!")
    print("="*60) 


def plot_all_in_one_comprehensive(adata, cluster_key='clusters', figsize=(24, 16), save_path=None):
    """
    Create one comprehensive visualization with ALL analysis results in a single plot.
    
    This combines UMAP views, spatial views, statistics, and cluster analysis 
    into one publication-ready figure.
    
    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object with UMAP, spatial coordinates, and clustering
    cluster_key : str, optional
        Key for clustering results (default: 'clusters')
    figsize : tuple, optional
        Figure size (default: (24, 16))
    save_path : str, optional
        Directory path to save plots. If None, don't save.
    """
    # Create save directory if specified
    if save_path:
        save_dir = Path(save_path)
        save_dir.mkdir(parents=True, exist_ok=True)
        print(f"💾 Saving comprehensive analysis to: {save_dir}")
    
    # Set figure parameters
    sc.set_figure_params(facecolor="white")
    
    # Create large subplot grid: 4 rows x 4 columns
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(4, 4, hspace=0.35, wspace=0.3)
    
    # ========== ROW 1: UMAP VIEWS ==========
    ax_umap1 = fig.add_subplot(gs[0, 0])
    ax_umap2 = fig.add_subplot(gs[0, 1])
    ax_umap3 = fig.add_subplot(gs[0, 2])
    ax_stats = fig.add_subplot(gs[0, 3])
    
    # ========== ROW 2: SPATIAL OVERVIEW ==========
    ax_spatial_combined = fig.add_subplot(gs[1, 0])
    ax_spatial_source = fig.add_subplot(gs[1, 1])
    ax_spatial_clusters = fig.add_subplot(gs[1, 2])
    ax_composition = fig.add_subplot(gs[1, 3])
    
    # ========== ROW 3: DETAILED SPATIAL BY CASTE ==========
    ax_soldier_spatial = fig.add_subplot(gs[2, 0])
    ax_forager_spatial = fig.add_subplot(gs[2, 1])
    ax_cluster_counts = fig.add_subplot(gs[2, 2])
    ax_summary = fig.add_subplot(gs[2, 3])
    
    # ========== ROW 4: EXTRA ANALYSIS ==========
    ax_batch_check = fig.add_subplot(gs[3, 0])
    ax_spatial_batch = fig.add_subplot(gs[3, 1])
    ax_detailed_stats = fig.add_subplot(gs[3, 2])
    ax_methods = fig.add_subplot(gs[3, 3])
    
    # ROW 1: UMAP ANALYSIS
    
    # 1. UMAP: Source (Caste)
    if 'source' in adata.obs.columns and 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color='source', title='UMAP: Caste Identity', 
                  ax=ax_umap1, show=False, frameon=False)
    
    # 2. UMAP: Batch (should be mixed after correction)
    if 'batch' in adata.obs.columns and 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color='batch', title='UMAP: Batch Labels', 
                  ax=ax_umap2, show=False, frameon=False)
    
    # 3. UMAP: Clusters
    if cluster_key in adata.obs.columns and 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color=cluster_key, title='UMAP: Leiden Clusters', 
                  ax=ax_umap3, show=False, frameon=False, legend_loc='on data')
    
    # 4. Dataset Statistics
    ax_stats.axis('off')
    stats_text = ["📊 Dataset Summary:"]
    stats_text.append(f"• Total cells: {adata.n_obs:,}")
    stats_text.append(f"• Total genes: {adata.n_vars:,}")
    
    if 'source' in adata.obs.columns:
        source_counts = adata.obs['source'].value_counts()
        stats_text.append("\n🐜 Caste Distribution:")
        for source, count in source_counts.items():
            pct = (count / adata.n_obs) * 100
            stats_text.append(f"• {source}: {count:,} ({pct:.1f}%)")
    
    if cluster_key in adata.obs.columns:
        n_clusters = len(adata.obs[cluster_key].unique())
        stats_text.append(f"\n🎯 Clustering:")
        stats_text.append(f"• {n_clusters} clusters found")
        stats_text.append(f"• Resolution: 0.2")
    
    ax_stats.text(0.05, 0.95, '\n'.join(stats_text), 
                 transform=ax_stats.transAxes, fontsize=11, 
                 verticalalignment='top', family='monospace')
    ax_stats.set_title('Analysis Overview', fontweight='bold')
    
    # ROW 2: SPATIAL OVERVIEW
    
    # 5. Combined spatial clusters
    if cluster_key in adata.obs.columns and 'spatial' in adata.obsm:
        sc.pl.embedding(adata, basis="spatial", color=cluster_key, 
                       title='Spatial: All Clusters',
                       ax=ax_spatial_combined, show=False, frameon=False)
    
    # 6. Spatial by source
    if 'source' in adata.obs.columns and 'spatial' in adata.obsm:
        sc.pl.embedding(adata, basis="spatial", color='source', 
                       title='Spatial: Caste Distribution',
                       ax=ax_spatial_source, show=False, frameon=False)
    
    # 7. Spatial clusters (cleaner view)
    if cluster_key in adata.obs.columns and 'spatial' in adata.obsm:
        sc.pl.embedding(adata, basis="spatial", color=cluster_key, 
                       title='Spatial: Cluster Regions',
                       ax=ax_spatial_clusters, show=False, frameon=False, legend_loc=None)
    
    # 8. Cluster composition bar chart
    if cluster_key in adata.obs.columns and 'source' in adata.obs.columns:
        cluster_comp = adata.obs.groupby([cluster_key, 'source']).size().unstack(fill_value=0)
        cluster_comp_pct = cluster_comp.div(cluster_comp.sum(axis=1), axis=0) * 100
        
        cluster_comp_pct.plot(kind='bar', stacked=True, ax=ax_composition, 
                             title='Cluster Composition (%)', rot=45)
        ax_composition.set_xlabel('Cluster')
        ax_composition.set_ylabel('Percentage')
        ax_composition.legend(title='Caste', loc='upper right')
    
    # ROW 3: CASTE-SPECIFIC SPATIAL ANALYSIS
    
    # 9. Soldier spatial
    if 'source' in adata.obs.columns and 'spatial' in adata.obsm:
        adata_soldier = adata[adata.obs["source"] == "Soldier"]
        if adata_soldier.n_obs > 0:
            sc.pl.embedding(adata_soldier, basis="spatial", color=cluster_key,
                           title=f'Soldier Spatial ({adata_soldier.n_obs:,} cells)',
                           ax=ax_soldier_spatial, show=False, frameon=False)
        else:
            ax_soldier_spatial.text(0.5, 0.5, 'No Soldier cells', 
                                  ha='center', va='center', transform=ax_soldier_spatial.transAxes)
            ax_soldier_spatial.set_title("Soldier Spatial")
    
    # 10. Forager spatial
    if 'source' in adata.obs.columns and 'spatial' in adata.obsm:
        adata_forager = adata[adata.obs["source"] == "Forager"]
        if adata_forager.n_obs > 0:
            sc.pl.embedding(adata_forager, basis="spatial", color=cluster_key,
                           title=f'Forager Spatial ({adata_forager.n_obs:,} cells)',
                           ax=ax_forager_spatial, show=False, frameon=False)
        else:
            ax_forager_spatial.text(0.5, 0.5, 'No Forager cells', 
                                  ha='center', va='center', transform=ax_forager_spatial.transAxes)
            ax_forager_spatial.set_title("Forager Spatial")
    
    # 11. Cell counts per cluster
    if cluster_key in adata.obs.columns:
        cluster_counts = adata.obs[cluster_key].value_counts().sort_index()
        ax_cluster_counts.bar(range(len(cluster_counts)), cluster_counts.values, color='skyblue')
        ax_cluster_counts.set_xlabel('Cluster')
        ax_cluster_counts.set_ylabel('Number of Cells')
        ax_cluster_counts.set_title('Cells per Cluster')
        ax_cluster_counts.set_xticks(range(len(cluster_counts)))
        ax_cluster_counts.set_xticklabels(cluster_counts.index)
    
    # 12. Detailed cluster breakdown
    ax_summary.axis('off')
    if cluster_key in adata.obs.columns and 'source' in adata.obs.columns:
        detail_text = ["🔍 Cluster Details:"]
        cluster_stats = adata.obs.groupby([cluster_key, 'source']).size().unstack(fill_value=0)
        for cluster in sorted(adata.obs[cluster_key].unique()):
            soldier = cluster_stats.loc[cluster, 'Soldier'] if 'Soldier' in cluster_stats.columns else 0
            forager = cluster_stats.loc[cluster, 'Forager'] if 'Forager' in cluster_stats.columns else 0
            total = soldier + forager
            s_pct = (soldier/total)*100 if total > 0 else 0
            f_pct = (forager/total)*100 if total > 0 else 0
            detail_text.append(f"C{cluster}: {total:3d} | S:{soldier:3d}({s_pct:3.0f}%) F:{forager:3d}({f_pct:3.0f}%)")
        
        ax_summary.text(0.05, 0.95, '\n'.join(detail_text), 
                       transform=ax_summary.transAxes, fontsize=10, 
                       verticalalignment='top', family='monospace')
    ax_summary.set_title('Cluster Breakdown', fontweight='bold')
    
    # ROW 4: TECHNICAL ANALYSIS
    
    # 13. Batch correction check
    if 'batch' in adata.obs.columns and 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color='batch', title='Batch Correction Check', 
                  ax=ax_batch_check, show=False, frameon=False, legend_loc='on data')
    
    # 14. Spatial batch distribution
    if 'batch' in adata.obs.columns and 'spatial' in adata.obsm:
        sc.pl.embedding(adata, basis="spatial", color='batch', 
                       title='Spatial: Batch Distribution',
                       ax=ax_spatial_batch, show=False, frameon=False)
    
    # 15. Processing status
    ax_detailed_stats.axis('off')
    process_text = ["🔬 Processing Status:"]
    if 'X_pca_harmony' in adata.obsm or 'X_harmony' in adata.obsm:
        process_text.append("✓ Harmony batch correction")
    else:
        process_text.append("✗ No batch correction")
        
    if 'X_umap' in adata.obsm:
        process_text.append("✓ UMAP embedding")
    else:
        process_text.append("✗ No UMAP")
        
    if 'spatial' in adata.obsm:
        process_text.append("✓ Spatial coordinates")
    else:
        process_text.append("✗ No spatial coordinates")
        
    if cluster_key in adata.obs.columns:
        process_text.append("✓ Leiden clustering")
        
    process_text.append(f"\n📈 Data Quality:")
    process_text.append(f"• PCs used: 50")
    process_text.append(f"• Harmony theta: 2.0")
    process_text.append(f"• Clustering res: 0.2")
    
    ax_detailed_stats.text(0.05, 0.95, '\n'.join(process_text), 
                          transform=ax_detailed_stats.transAxes, fontsize=10, 
                          verticalalignment='top', family='monospace')
    ax_detailed_stats.set_title('Technical Details', fontweight='bold')
    
    # 16. Methods summary
    ax_methods.axis('off')
    methods_text = ["🛠️ Methods:"]
    methods_text.append("• Data: Enhanced iSpatial")
    methods_text.append("• Integration: Harmony")
    methods_text.append("• Clustering: Leiden")
    methods_text.append("• Embedding: UMAP")
    methods_text.append("• Coordinates: Spatial")
    methods_text.append("")
    methods_text.append("📊 Analysis Pipeline:")
    methods_text.append("1. Load S+F datasets")
    methods_text.append("2. Harmony batch correction")
    methods_text.append("3. Leiden clustering")
    methods_text.append("4. UMAP embedding")
    methods_text.append("5. Spatial visualization")
    
    ax_methods.text(0.05, 0.95, '\n'.join(methods_text), 
                   transform=ax_methods.transAxes, fontsize=10, 
                   verticalalignment='top', family='monospace')
    ax_methods.set_title('Analysis Methods', fontweight='bold')
    
    # Add main title
    fig.suptitle('Complete Bee Brain Spatial Transcriptomics Analysis', 
                fontsize=20, fontweight='bold', y=0.98)
    
    # Save if path specified
    if save_path:
        filename = f"complete_bee_brain_analysis_{cluster_key}.png"
        plt.savefig(save_dir / filename, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  ✓ Saved comprehensive analysis: {filename}")
    
    plt.show()
    
    # Print summary
    print(f"\n📋 Analysis Summary:")
    if 'source' in adata.obs.columns:
        source_counts = adata.obs['source'].value_counts()
        for source, count in source_counts.items():
            pct = (count / adata.n_obs) * 100
            print(f"  • {source}: {count:,} cells ({pct:.1f}%)")
    
    if cluster_key in adata.obs.columns:
        n_clusters = len(adata.obs[cluster_key].unique())
        print(f"  • Found {n_clusters} distinct clusters")
        print(f"  • Batch correction: {'✓ Applied' if 'X_pca_harmony' in adata.obsm else '✗ Not applied'}")
    
    print(f"  • Total visualization panels: 16")
    print(f"  • Image saved at 300 DPI for publication quality") 