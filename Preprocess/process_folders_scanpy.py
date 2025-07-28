#!/usr/bin/env python3
"""
Script to process two folders containing scanpy_coordinates.csv and scanpy_counts.csv
Each folder is processed independently to create normalized, clustered h5ad files.
"""

import os
import pandas as pd
import numpy as np
import scanpy as sc
from anndata import AnnData
from datetime import datetime

def process_folder(folder_path, folder_name):
    """
    Process a single folder containing scanpy_coordinates.csv and scanpy_counts.csv
    
    Parameters:
    folder_path (str): Path to the folder containing the CSV files
    folder_name (str): Name of the folder (used as source identifier)
    
    Returns:
    AnnData: Processed and clustered AnnData object
    """
    print(f"Processing folder: {folder_name}")
    
    # 1) Build paths to CSV files
    coord_path = os.path.join(folder_path, "scanpy_coordinates.csv")
    counts_path = os.path.join(folder_path, "scanpy_counts.csv")
    
    # Check if files exist
    if not os.path.exists(coord_path):
        raise FileNotFoundError(f"scanpy_coordinates.csv not found in {folder_path}")
    if not os.path.exists(counts_path):
        raise FileNotFoundError(f"scanpy_counts.csv not found in {folder_path}")
    
    # 2) Load coordinates (cells × [x,y,size,…])
    print(f"  Loading coordinates from {coord_path}")
    coords = pd.read_csv(coord_path, index_col="cell_id")
    
    # 3) Load counts
    print(f"  Loading counts from {counts_path}")
    raw_counts = pd.read_csv(counts_path, index_col=0)
    
    # Check if we need to transpose (similar logic from original code)
    if raw_counts.index.intersection(coords.index).size > 0:
        counts = raw_counts
    else:
        # Assume rows=genes, cols=cells, so transpose
        counts = raw_counts.T
    
    # 4) Align on the same set of cells
    common = counts.index.intersection(coords.index)
    if common.size == 0:
        raise ValueError("No overlap between count rows and coordinate rows. "
                         "Check whether you need to transpose the counts.")
    
    print(f"  Found {len(common)} common cells")
    counts = counts.loc[common]
    coords = coords.loc[common]
    
    # 5) Build AnnData
    adata = AnnData(
        X=counts.values,
        obs=pd.DataFrame(index=common),
        var=pd.DataFrame(index=counts.columns),
        obsm={"spatial": coords[["x", "y"]].values}
    )
    
    # Add metadata
    adata.obs["size"] = coords["size"].values
    adata.obs["total_counts"] = adata.X.sum(axis=1)
    adata.obs["source"] = folder_name
    
    # 6) Optional filter (keep cells with > 2 total counts)
    adata = adata[adata.obs["total_counts"] > 2].copy()
    print(f"  After filtering: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # 7) Normalize counts by size
    print("  Normalizing counts by cell size")
    cell_sizes = adata.obs['size'].values
    normalized_counts = (adata.X / cell_sizes[:, None]) * np.mean(cell_sizes)
    adata.X = normalized_counts
    adata.obs['n_counts'] = normalized_counts.sum(axis=1)
    
    # 8) Log1p transformation
    print("  Applying log1p transformation")
    sc.pp.log1p(adata)
    
    # 9) PCA, neighbors, UMAP, and clustering
    print("  Computing PCA")
    sc.pp.pca(adata, n_comps=15)
    
    print("  Computing neighbors")
    sc.pp.neighbors(adata)
    
    print("  Computing UMAP")
    sc.tl.umap(adata, min_dist=0.2)
    
    print("  Computing Leiden clustering")
    sc.tl.leiden(adata, key_added="clusters", resolution=0.5)
    
    return adata

def main():
    """
    Main function to process two folders and save results
    """
    # Get today's date for filename
    today = datetime.now().strftime("%Y%m%d")
    
    # Define folder paths - modify these as needed
    base_path = "/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/data/SvsF_col116"
    folder_names = ["f11", "s31"]  # Modify these folder names as needed
    
    # Process each folder independently
    for folder_name in folder_names:
        folder_path = os.path.join(base_path, folder_name)
        
        if not os.path.exists(folder_path):
            print(f"Warning: Folder {folder_path} does not exist. Skipping.")
            continue
        
        try:
            # Process the folder
            adata = process_folder(folder_path, folder_name)
            
            # Save the processed h5ad file with today's date
            output_filename = f"{folder_name}_processed_{today}.h5ad"
            output_path = os.path.join(folder_path, output_filename)
            
            print(f"  Saving processed data to {output_path}")
            adata.write_h5ad(output_path)
            
            print(f"Successfully processed {folder_name}")
            print(f"  Final data shape: {adata.n_obs} cells × {adata.n_vars} genes")
            print(f"  Clusters found: {len(adata.obs['clusters'].unique())}")
            print("")
            
        except Exception as e:
            print(f"Error processing folder {folder_name}: {str(e)}")
            print("")

if __name__ == "__main__":
    # Set scanpy settings
    sc.settings.verbosity = 1  # Reduce scanpy verbosity
    sc.settings.set_figure_params(dpi=80, facecolor='white')
    
    main() 