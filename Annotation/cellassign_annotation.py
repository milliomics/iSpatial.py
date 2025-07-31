import pandas as pd
import scanpy as sc
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Import scvi-tools for CellAssign
try:
    import scvi
    import torch
    from scvi.external import CellAssign  # Correct import from scvi.external
    print("✅ scvi-tools imported successfully")
except ImportError as e:
    print("❌ scvi-tools not installed. Please install with:")
    print("   pip install scvi-tools")
    print("   or: conda install -c conda-forge scvi-tools")
    raise e

def prepare_marker_matrix(markers_df: pd.DataFrame, 
                         adata_genes: list,
                         cell_type_col: str = 'general_type') -> tuple:
    """
    Prepare marker gene matrix for CellAssign from markers.xlsx.
    
    Parameters
    ----------
    markers_df : pd.DataFrame
        Marker dataframe from markers.xlsx
    adata_genes : list
        List of gene names from adata.var_names
    cell_type_col : str
        Column to use for cell types ('general_type' or 'specific_type')
        
    Returns
    -------
    tuple
        (marker_matrix, cell_types, marker_genes)
    """
    print(f"🧬 Preparing marker matrix using '{cell_type_col}' column...")
    
    # Get unique cell types
    cell_types = sorted(markers_df[cell_type_col].dropna().unique())
    print(f"   Found {len(cell_types)} cell types: {cell_types}")
    
    # Get marker genes that exist in adata
    marker_genes_in_data = []
    marker_genes_not_found = []
    
    for gene_id in markers_df['gene_id'].unique():
        if gene_id in adata_genes:
            marker_genes_in_data.append(gene_id)
        else:
            marker_genes_not_found.append(gene_id)
    
    print(f"   Marker genes found in data: {len(marker_genes_in_data)}/{len(markers_df['gene_id'].unique())}")
    if marker_genes_not_found:
        print(f"   Genes not found: {marker_genes_not_found[:5]}{'...' if len(marker_genes_not_found) > 5 else ''}")
    
    # Create binary marker matrix: genes × cell_types
    marker_matrix = np.zeros((len(marker_genes_in_data), len(cell_types)))
    
    for i, gene_id in enumerate(marker_genes_in_data):
        # Find which cell types this gene is a marker for
        gene_rows = markers_df[markers_df['gene_id'] == gene_id]
        for _, row in gene_rows.iterrows():
            cell_type = row[cell_type_col]
            if pd.notna(cell_type) and cell_type in cell_types:
                j = cell_types.index(cell_type)
                marker_matrix[i, j] = 1
    
    # Convert to DataFrame for easier handling
    marker_matrix_df = pd.DataFrame(
        marker_matrix,
        index=marker_genes_in_data,
        columns=cell_types
    )
    
    print(f"   Marker matrix shape: {marker_matrix_df.shape}")
    print(f"   Markers per cell type:")
    
    # Filter out cell types with no markers
    valid_cell_types = []
    for cell_type in cell_types:
        n_markers = (marker_matrix_df[cell_type] > 0).sum()
        print(f"     {cell_type}: {n_markers} markers")
        if n_markers > 0:
            valid_cell_types.append(cell_type)
        else:
            print(f"       ⚠️  Removing {cell_type} (no markers found)")
    
    # Filter marker matrix to only valid cell types
    marker_matrix_filtered = marker_matrix_df[valid_cell_types].copy()
    
    print(f"   Final marker matrix shape: {marker_matrix_filtered.shape}")
    print(f"   Using {len(valid_cell_types)} cell types with markers")
    
    return marker_matrix_filtered, valid_cell_types, marker_genes_in_data

def run_cellassign(adata_path: str,
                  markers_path: str,
                  output_path: str = None,
                  cell_type_col: str = 'general_type',
                  max_epochs: int = 100) -> sc.AnnData:
    """
    Run CellAssign for cell type annotation.
    
    Parameters
    ----------
    adata_path : str
        Path to anndata object (.h5ad file)
    markers_path : str
        Path to markers.xlsx file
    output_path : str, optional
        Path to save results
    cell_type_col : str
        Cell type column to use ('general_type' or 'specific_type')
    max_epochs : int
        Maximum training epochs
        
    Returns
    -------
    sc.AnnData
        Annotated adata object with CellAssign predictions
    """
    
    print("🔬 CellAssign Cell Type Annotation")
    print("=" * 60)
    
    # Load data
    print(f"📂 Loading data...")
    adata = sc.read_h5ad(adata_path)
    markers_df = pd.read_excel(markers_path)
    
    print(f"   AnnData: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print(f"   Markers: {len(markers_df)} marker entries")
    
    # Prepare marker matrix
    marker_matrix, cell_types, marker_genes = prepare_marker_matrix(
        markers_df, adata.var_names.tolist(), cell_type_col
    )
    
    if len(marker_genes) == 0:
        raise ValueError("No marker genes found in the data! Check gene IDs.")
    
    # CellAssign requires adata to contain ONLY marker genes
    print(f"\n🎯 Subsetting data to marker genes only...")
    marker_genes_in_data = [g for g in marker_genes if g in adata.var_names]
    print(f"   Using {len(marker_genes_in_data)} marker genes out of {len(marker_genes)} total")
    
    # Subset adata to marker genes only (CellAssign requirement)
    adata_subset = adata[:, marker_genes_in_data].copy()
    
    # Prepare marker matrix to match adata gene order exactly
    marker_matrix_final = marker_matrix.loc[marker_genes_in_data].copy()
    
    print(f"   Final data shape: {adata_subset.shape}")
    print(f"   Final marker matrix: {marker_matrix_final.shape}")
    print(f"   Non-zero entries: {(marker_matrix_final > 0).sum().sum()}")
    
    # Ensure data is raw counts (not log-transformed)
    if adata_subset.X.max() < 50:  # Likely log-transformed
        print("   Warning: Data appears to be log-transformed. Converting back to linear scale...")
        # More careful conversion to avoid negative values
        adata_subset.X = np.expm1(adata_subset.X)
        # Ensure no negative values (set to small positive value)
        adata_subset.X[adata_subset.X < 0] = 1e-6
        print(f"   Data range after conversion: {adata_subset.X.min():.6f} - {adata_subset.X.max():.1f}")
    
    # Calculate size factors (required by CellAssign)
    print(f"\n📏 Calculating size factors...")
    library_size = np.array(adata_subset.X.sum(1)).flatten()
    # Ensure library sizes are positive
    library_size[library_size <= 0] = 1e-6
    size_factors = library_size / np.mean(library_size)
    # Ensure size factors are positive and reasonable
    size_factors = np.clip(size_factors, 1e-6, 1000)
    adata_subset.obs["size_factor"] = size_factors
    print(f"   Size factor range: {adata_subset.obs['size_factor'].min():.3f} - {adata_subset.obs['size_factor'].max():.3f}")
    
    # Validate data before CellAssign
    print(f"\n🔍 Validating data for CellAssign...")
    print(f"   Data shape: {adata_subset.shape}")
    print(f"   Data type: {adata_subset.X.dtype}")
    print(f"   Data range: {adata_subset.X.min():.6f} - {adata_subset.X.max():.1f}")
    print(f"   Any NaN in data: {np.isnan(adata_subset.X).any()}")
    print(f"   Any infinite in data: {np.isinf(adata_subset.X).any()}")
    print(f"   Marker matrix shape: {marker_matrix_final.shape}")
    print(f"   Any NaN in markers: {marker_matrix_final.isna().any().any()}")
    
    # Ensure data is sparse if it's a large matrix
    if not hasattr(adata_subset.X, 'toarray'):
        print("   Converting to sparse matrix...")
        from scipy.sparse import csr_matrix
        adata_subset.X = csr_matrix(adata_subset.X)
    
    # Setup CellAssign
    print(f"\n⚙️  Setting up CellAssign model...")
    scvi.settings.verbosity = 1
    
    # Force CPU usage (no GPU)
    import torch
    if torch.cuda.is_available():
        print("   GPU available but forcing CPU usage...")
    scvi.settings.device = "cpu"
    
    # Setup anndata with size factor
    CellAssign.setup_anndata(adata_subset, size_factor_key="size_factor")
    
    # Create model with proper marker matrix format
    model = CellAssign(adata_subset, marker_matrix_final)
    
    print(f"   Model created with {len(cell_types)} cell types")
    
    # Train model
    print(f"\n🚀 Training CellAssign model...")
    print(f"   Max epochs: {max_epochs}")
    
    try:
        # Use more conservative training parameters
        model.train(
            max_epochs=max_epochs, 
            check_val_every_n_epoch=50
        )
    except Exception as train_error:
        print(f"   ❌ Training failed: {train_error}")
        print("   🔄 Trying with fewer epochs...")
        try:
            model.train(
                max_epochs=min(100, max_epochs),
                check_val_every_n_epoch=25
            )
        except Exception as second_error:
            print(f"   ❌ Second attempt failed: {second_error}")
            raise second_error
    
    # Get predictions
    print(f"\n📊 Getting cell type predictions...")
    try:
        # Get hard predictions (cell type labels)
        predictions = model.predict()
        print(f"   Prediction result type: {type(predictions)}")
        print(f"   Prediction shape: {predictions.shape if hasattr(predictions, 'shape') else len(predictions)}")
        
        # Handle different prediction formats
        if isinstance(predictions, pd.DataFrame):
            print(f"   Prediction columns: {predictions.columns.tolist()}")
            
            if predictions.shape[1] == len(cell_types):
                # This is a probability matrix - convert to labels by taking argmax
                print("   Converting probabilities to cell type labels...")
                max_prob_indices = predictions.values.argmax(axis=1)
                cell_type_predictions = [predictions.columns[i] for i in max_prob_indices]
                
                # Also get the max probabilities for confidence
                max_probs = predictions.values.max(axis=1)
                print(f"   Confidence range: {max_probs.min():.3f} - {max_probs.max():.3f}")
                print(f"   Mean confidence: {max_probs.mean():.3f}")
                
                # Show prediction distribution
                pred_counts = pd.Series(cell_type_predictions).value_counts()
                print(f"   Prediction distribution:")
                for cell_type, count in pred_counts.items():
                    pct = count / len(cell_type_predictions) * 100
                    print(f"     {cell_type}: {count} cells ({pct:.1f}%)")
                    
            elif predictions.shape[1] == 1:
                # Single column DataFrame
                cell_type_predictions = predictions.iloc[:, 0]
            else:
                # Multiple columns - try to find a 'prediction' or 'celltype' column
                if 'prediction' in predictions.columns:
                    cell_type_predictions = predictions['prediction']
                elif 'celltype' in predictions.columns:
                    cell_type_predictions = predictions['celltype']
                else:
                    # Use the first column as predictions
                    cell_type_predictions = predictions.iloc[:, 0]
        else:
            # Assume it's already a series or array
            cell_type_predictions = predictions
        
        # Add predictions to original adata
        cell_type_col_name = f'cellassign_{cell_type_col}'
        adata.obs[cell_type_col_name] = cell_type_predictions
        
        # Store prediction probabilities (these are the same as the hard predictions for CellAssign)
        try:
            if isinstance(predictions, pd.DataFrame) and predictions.shape[1] == len(cell_types):
                # The predictions DataFrame already contains probabilities
                print("   Storing prediction probabilities...")
                for cell_type in cell_types:
                    if cell_type in predictions.columns:
                        adata.obs[f'{cell_type_col_name}_prob_{cell_type}'] = predictions[cell_type].values
                        
                # Also store the maximum probability as confidence score
                max_probs = predictions.values.max(axis=1)
                adata.obs[f'{cell_type_col_name}_confidence'] = max_probs
            else:
                print("   Predictions not in expected probability format")
        except Exception as prob_error:
            print(f"   Could not store prediction probabilities: {prob_error}")
            
    except Exception as pred_error:
        print(f"   ❌ Prediction failed: {pred_error}")
        raise pred_error
    
    # Results summary
    print(f"\n✅ CellAssign annotation completed!")
    print(f"   Results in: adata.obs['{cell_type_col_name}']")
    
    prediction_counts = adata.obs[cell_type_col_name].value_counts()
    print(f"\n📈 Cell type distribution:")
    for cell_type, count in prediction_counts.items():
        pct = count / len(adata) * 100
        print(f"   {cell_type}: {count:,} cells ({pct:.1f}%)")
    
    # Save results
    if output_path:
        print(f"\n💾 Saving results to: {output_path}")
        adata.write_h5ad(output_path)
        
        # Save detailed results
        results_dir = Path(output_path).parent
        results_df = adata.obs[[cell_type_col_name]].copy()
        results_df.to_excel(results_dir / f"cellassign_{cell_type_col}_results.xlsx")
        
        # Save marker matrix used
        marker_matrix_final.to_excel(results_dir / f"cellassign_{cell_type_col}_marker_matrix.xlsx")
        
        print("   ✓ Results saved")
    
    return adata

def main():
    """Example usage"""
    
    # ==================== CONFIGURE YOUR PATHS HERE ====================
    
    # Input paths
    adata_path = "/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/code/output_ispatial/f11_s31/f11_s31_col116_w0.2.h5ad"
    markers_path = "/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/code/annotation/markers_Mu.xlsx"
    
    # Output path
    output_path = "/Users/farah/Library/CloudStorage/GoogleDrive-qianluf2@illinois.edu/My Drive/Han_lab_Drive/p5_SvsF/code/annotation/annotation_output/f11_s31_col116_cellassign_20250730.h5ad"
    
    # ====================================================================
    
    # Check if files exist
    if not Path(adata_path).exists():
        print(f"❌ AnnData file not found: {adata_path}")
        return
        
    if not Path(markers_path).exists():
        print(f"❌ Markers file not found: {markers_path}")
        return
    
    # Create output directory
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    
    # User selection for annotation type
    print("\n🎯 CellAssign Annotation Type Selection")
    print("=" * 50)
    print("1. General cell types only")
    print("2. Specific cell types only")
    print("3. Both general and specific types")
    print("=" * 50)
    
    while True:
        try:
            choice = input("Please select an option (1, 2, or 3): ").strip()
            if choice in ['1', '2', '3']:
                break
            else:
                print("❌ Invalid choice. Please enter 1, 2, or 3.")
        except KeyboardInterrupt:
            print("\n❌ Operation cancelled by user.")
            return
        except Exception:
            print("❌ Invalid input. Please enter 1, 2, or 3.")
    
    # Run selected annotation(s)
    if choice in ['1', '3']:  # General types
        print("\n🎯 Running CellAssign for GENERAL cell types...")
        try:
            adata_general = run_cellassign(
                adata_path=adata_path,
                markers_path=markers_path,
                output_path=output_path.replace('.h5ad', '_general.h5ad'),
                cell_type_col='general_type',
                max_epochs=100  # Reduced for CPU training
            )
            print("✅ General type annotation completed!")
            
        except Exception as e:
            print(f"❌ Error in general type annotation: {e}")
    
    if choice in ['2', '3']:  # Specific types
        if choice == '3':
            print("\n" + "="*60)
        print("🎯 Running CellAssign for SPECIFIC cell types...")
        try:
            adata_specific = run_cellassign(
                adata_path=adata_path,
                markers_path=markers_path,
                output_path=output_path.replace('.h5ad', '_specific.h5ad'),
                cell_type_col='specific_type',
                max_epochs=100  # Reduced for CPU training
            )
            print("✅ Specific type annotation completed!")
            
        except Exception as e:
            print(f"❌ Error in specific type annotation: {e}")
    
    print(f"\n🎉 All selected annotations completed!")
    print(f"📁 Results saved in: {Path(output_path).parent}")

if __name__ == "__main__":
    main() 