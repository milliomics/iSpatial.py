#!/usr/bin/env python3
"""
Cell Type Annotation Toolkit for Bee Brain Spatial Data
========================================================

Comprehensive toolkit for annotating cell types in enhanced spatial transcriptomics
data using marker genes from literature.

Key Features:
- Marker gene-based scoring for multiple cell types
- Spatial and UMAP visualization
- Confidence-based assignment
- Quality control and validation
- Support for hierarchical annotations (broad -> specific)

Author: Bee Brain Analysis Pipeline
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import warnings
warnings.filterwarnings('ignore')

# Set scanpy settings
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=80, facecolor='white')

class BeebrainCellTypeAnnotator:
    """
    Comprehensive cell type annotation for bee brain spatial data.
    """
    
    def __init__(self, adata: ad.AnnData, marker_genes: Dict[str, List[str]] = None):
        """
        Initialize the annotator.
        
        Parameters
        ----------
        adata : anndata.AnnData
            Enhanced spatial data from iSpatial pipeline
        marker_genes : dict, optional
            Dictionary of cell type -> marker genes. If None, uses default bee brain markers.
        """
        self.adata = adata.copy()
        self.marker_genes = marker_genes or self._get_default_bee_markers()
        self.annotation_results = {}
        
        # Check data quality
        self._validate_data()
        
    def _validate_data(self):
        """Validate input data quality."""
        print(f"📊 Data validation:")
        print(f"  - Cells: {self.adata.n_obs:,}")
        print(f"  - Genes: {self.adata.n_vars:,}")
        
        # Check for required coordinates
        if 'spatial' in self.adata.obsm:
            print(f"  ✓ Spatial coordinates available")
        else:
            print(f"  ⚠️  No spatial coordinates found")
            
        if 'X_umap' in self.adata.obsm:
            print(f"  ✓ UMAP coordinates available")
        else:
            print(f"  ⚠️  No UMAP coordinates found")
            
        # Check expression range
        if hasattr(self.adata.X, 'max'):
            max_expr = self.adata.X.max()
        else:
            max_expr = self.adata.X.toarray().max()
            
        print(f"  - Expression range: 0 to {max_expr:.2f}")
        if max_expr > 20:
            print(f"    ⚠️  High values suggest raw counts - consider log transformation")
        elif max_expr < 10:
            print(f"    ✓ Values suggest log-transformed data")
            
    def _get_default_bee_markers(self) -> Dict[str, List[str]]:
        """
        Default bee brain marker genes from literature.
        
        Returns
        -------
        dict
            Cell type -> marker genes mapping
        """
        # These are example markers - you should replace with your literature-curated list
        default_markers = {
            # Broad cell types
            'Neurons': [
                'elav', 'nSyb', 'Syt1', 'CaMKII', 'futsch', 'brp',
                'nrv2', 'Syt4', 'cac', 'para', 'shab', 'Shaw'
            ],
            'Glia': [
                'repo', 'Gs2', 'Eaat1', 'wrapper', 'nrv2', 'Gli',
                'NetA', 'NetB', 'wrapper', 'moody'
            ],
            'Hemocytes': [
                'Hml', 'srp', 'gcm', 'lz', 'Pxn', 'NimC1',
                'PGRP-LC', 'Dscam', 'eater', 'crq'
            ],
            
            # Specific Kenyon cell types
            'sKC': [  # Small Kenyon cells
                'mb247', 'FasII', 'arm', 'Lac', 'Synapsin',
                'OK107', '201Y', 'c739', 'MB-DsRed'
            ],
            'lKC': [  # Large Kenyon cells  
                'mb247', 'FasII', 'VGlut', 'ChAT', 'Gad1',
                'OK107', 'c305a', 'MB-DsRed', 'ple'
            ],
            'mKC': [  # Medium Kenyon cells
                'mb247', 'FasII', 'Synapsin', 'VGlut', 'ChAT',
                'OK107', 'Gad1', 'MB-DsRed'
            ],
            
            # Olfactory system
            'OLC': [  # Olfactory lobe cells
                'Orco', 'Ir8a', 'Ir25a', 'Ir76b', 'Or22a',
                'Or47b', 'Or88a', 'GluClalpha', 'Rdl'
            ],
            'PN': [   # Projection neurons
                'GH146', 'Mz19', 'NP225', 'NP2631', 'NP3056',
                'GluClalpha', 'Rdl', 'ChAT', 'VGlut'
            ],
            
            # Other brain regions
            'VisualSystem': [
                'rh1', 'rh3', 'rh4', 'rh5', 'rh6', 'norpA',
                'trp', 'inaC', 'ninaE', 'Rh1'
            ],
            'MotorNeurons': [
                'ChAT', 'VGlut', 'Gad1', 'OK371', 'eve',
                'en', 'hb', 'ftz', 'even-skipped'
            ],
            'Interneurons': [
                'Gad1', 'VGlut', 'ChAT', 'DvGlut', 'Rdl',
                'GluClalpha', 'nAChRalpha', 'Cha'
            ]
        }
        
        print(f"📚 Loaded {len(default_markers)} default cell type marker sets")
        print(f"   Cell types: {', '.join(default_markers.keys())}")
        print(f"   Note: Replace with your literature-curated markers!")
        
        return default_markers
    
    def load_marker_genes(self, marker_file: Union[str, Dict[str, List[str]]]):
        """
        Load marker genes from file or dictionary.
        
        Parameters
        ----------
        marker_file : str or dict
            Path to CSV file or dictionary of cell type -> marker genes
        """
        if isinstance(marker_file, dict):
            self.marker_genes = marker_file
        elif isinstance(marker_file, str):
            # Load from CSV file
            df = pd.read_csv(marker_file)
            if 'cell_type' in df.columns and 'marker_genes' in df.columns:
                self.marker_genes = {}
                for _, row in df.iterrows():
                    cell_type = row['cell_type']
                    markers = row['marker_genes'].split(',')
                    markers = [m.strip() for m in markers]
                    self.marker_genes[cell_type] = markers
            else:
                raise ValueError("CSV must have 'cell_type' and 'marker_genes' columns")
        
        print(f"✓ Loaded {len(self.marker_genes)} cell type marker sets")
        
    def check_marker_availability(self, verbose: bool = True) -> Dict[str, Dict[str, float]]:
        """
        Check how many marker genes are available in the dataset.
        
        Parameters
        ----------
        verbose : bool
            Whether to print detailed results
            
        Returns
        -------
        dict
            Cell type -> marker availability statistics
        """
        availability = {}
        all_genes = set(self.adata.var_names)
        
        if verbose:
            print(f"\n🔍 Marker Gene Availability Analysis:")
            print(f"{'Cell Type':<20} {'Available':<10} {'Total':<8} {'Percentage':<12} {'Missing Genes'}")
            print("-" * 80)
        
        for cell_type, markers in self.marker_genes.items():
            available_markers = [m for m in markers if m in all_genes]
            missing_markers = [m for m in markers if m not in all_genes]
            
            availability[cell_type] = {
                'available': len(available_markers),
                'total': len(markers),
                'percentage': len(available_markers) / len(markers) * 100,
                'available_genes': available_markers,
                'missing_genes': missing_markers
            }
            
            if verbose:
                pct = availability[cell_type]['percentage']
                missing_str = ', '.join(missing_markers[:3])  # Show first 3 missing
                if len(missing_markers) > 3:
                    missing_str += f"... (+{len(missing_markers)-3} more)"
                    
                print(f"{cell_type:<20} {len(available_markers):<10} {len(markers):<8} {pct:<11.1f}% {missing_str}")
        
        return availability
    
    def score_cell_types(self, method: str = 'mean_expression', 
                        score_name_suffix: str = '_score') -> None:
        """
        Score cells for each cell type based on marker gene expression.
        
        Parameters
        ----------
        method : str
            Scoring method: 'mean_expression', 'scanpy_score', 'median_expression'
        score_name_suffix : str
            Suffix for score column names
        """
        print(f"\n🧮 Computing cell type scores using method: {method}")
        
        # Check marker availability first
        availability = self.check_marker_availability(verbose=False)
        
        for cell_type, markers in self.marker_genes.items():
            available_markers = availability[cell_type]['available_genes']
            
            if len(available_markers) == 0:
                print(f"  ⚠️  {cell_type}: No markers available - skipping")
                continue
                
            print(f"  - {cell_type}: {len(available_markers)}/{len(markers)} markers available")
            
            score_name = f"{cell_type}{score_name_suffix}"
            
            if method == 'mean_expression':
                # Simple mean of marker gene expression
                marker_expr = self.adata[:, available_markers].X
                if hasattr(marker_expr, 'toarray'):
                    marker_expr = marker_expr.toarray()
                scores = np.mean(marker_expr, axis=1)
                self.adata.obs[score_name] = scores
                
            elif method == 'scanpy_score':
                # Use scanpy's gene scoring (like Seurat's AddModuleScore)
                sc.tl.score_genes(self.adata, available_markers, score_name=score_name)
                
            elif method == 'median_expression':
                # Median of marker gene expression (more robust to outliers)
                marker_expr = self.adata[:, available_markers].X
                if hasattr(marker_expr, 'toarray'):
                    marker_expr = marker_expr.toarray()
                scores = np.median(marker_expr, axis=1)
                self.adata.obs[score_name] = scores
                
            else:
                raise ValueError(f"Unknown scoring method: {method}")
        
        # Store scoring info
        self.annotation_results['scoring_method'] = method
        self.annotation_results['score_columns'] = [f"{ct}{score_name_suffix}" 
                                                   for ct in self.marker_genes.keys()]
        
        print(f"✓ Cell type scoring completed")
    
    def assign_cell_types(self, score_suffix: str = '_score', 
                         confidence_threshold: float = 0.1,
                         assign_method: str = 'winner_takes_all') -> None:
        """
        Assign cell types based on scores.
        
        Parameters
        ----------
        score_suffix : str
            Suffix of score columns to use
        confidence_threshold : float
            Minimum score difference for confident assignment
        assign_method : str
            Assignment method: 'winner_takes_all', 'threshold_based'
        """
        print(f"\n🎯 Assigning cell types using method: {assign_method}")
        
        # Get score columns
        score_cols = [col for col in self.adata.obs.columns if col.endswith(score_suffix)]
        
        if len(score_cols) == 0:
            print(f"  ⚠️  No score columns found with suffix '{score_suffix}'")
            print(f"  Available columns: {list(self.adata.obs.columns)}")
            return
            
        print(f"  Using {len(score_cols)} cell type scores")
        
        # Extract scores
        scores_df = self.adata.obs[score_cols].copy()
        scores_df.columns = [col.replace(score_suffix, '') for col in score_cols]
        
        if assign_method == 'winner_takes_all':
            # Assign to highest scoring cell type
            max_scores = scores_df.max(axis=1)
            assignments = scores_df.idxmax(axis=1)
            
            # Calculate confidence (difference between top 2 scores)
            sorted_scores = np.sort(scores_df.values, axis=1)
            if scores_df.shape[1] > 1:
                confidence = sorted_scores[:, -1] - sorted_scores[:, -2]  # Top - Second
            else:
                confidence = sorted_scores[:, -1]  # Only one cell type
            
            # Low confidence assignments
            low_confidence = confidence < confidence_threshold
            assignments[low_confidence] = 'Unassigned'
            
            # Very low scores
            very_low = max_scores < np.percentile(max_scores, 10)  # Bottom 10%
            assignments[very_low] = 'Unassigned'
            
        else:
            raise ValueError(f"Unknown assignment method: {assign_method}")
        
        # Store results
        self.adata.obs['cell_type_assignment'] = assignments.astype('category')
        self.adata.obs['assignment_confidence'] = confidence
        self.adata.obs['max_score'] = max_scores
        
        # Summary
        assignment_counts = assignments.value_counts()
        print(f"\n📊 Assignment Summary:")
        for cell_type, count in assignment_counts.items():
            pct = count / len(assignments) * 100
            print(f"  {cell_type:<20}: {count:5d} cells ({pct:5.1f}%)")
        
        unassigned_pct = (assignments == 'Unassigned').sum() / len(assignments) * 100
        print(f"\n  Confidence threshold: {confidence_threshold}")
        print(f"  Unassigned cells: {unassigned_pct:.1f}%")
        
        self.annotation_results['assignments'] = assignments
        self.annotation_results['confidence'] = confidence
        
    def plot_score_distributions(self, score_suffix: str = '_score', 
                               figsize: Tuple[int, int] = (15, 10),
                               save_path: Optional[str] = None) -> None:
        """
        Plot distribution of cell type scores.
        """
        score_cols = [col for col in self.adata.obs.columns if col.endswith(score_suffix)]
        
        if len(score_cols) == 0:
            print(f"No score columns found with suffix '{score_suffix}'")
            return
            
        n_types = len(score_cols)
        n_cols = 4
        n_rows = (n_types + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
        axes = axes.flatten() if n_rows > 1 else [axes]
        
        for i, score_col in enumerate(score_cols):
            cell_type = score_col.replace(score_suffix, '')
            scores = self.adata.obs[score_col]
            
            axes[i].hist(scores, bins=50, alpha=0.7, edgecolor='black')
            axes[i].set_title(f'{cell_type} Score Distribution')
            axes[i].set_xlabel('Score')
            axes[i].set_ylabel('Number of Cells')
            axes[i].grid(True, alpha=0.3)
        
        # Hide unused subplots
        for i in range(n_types, len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved score distributions to {save_path}")
        
        plt.show()
    
    def plot_umap_annotations(self, color_by: str = 'cell_type_assignment',
                            figsize: Tuple[int, int] = (12, 10),
                            save_path: Optional[str] = None) -> None:
        """
        Plot cell type annotations on UMAP.
        """
        if 'X_umap' not in self.adata.obsm:
            print("⚠️  UMAP coordinates not available")
            return
            
        if color_by not in self.adata.obs.columns:
            print(f"⚠️  Column '{color_by}' not found in adata.obs")
            return
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Plot UMAP with annotations
        sc.pl.umap(self.adata, color=color_by, ax=ax, show=False, 
                  legend_loc='right margin', frameon=False)
        
        plt.title(f'UMAP: {color_by.replace("_", " ").title()}', fontsize=14)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved UMAP plot to {save_path}")
        
        plt.show()
    
    def plot_spatial_annotations(self, color_by: str = 'cell_type_assignment',
                               figsize: Tuple[int, int] = (12, 10),
                               save_path: Optional[str] = None,
                               point_size: Optional[float] = None) -> None:
        """
        Plot cell type annotations on spatial coordinates.
        """
        if 'spatial' not in self.adata.obsm:
            print("⚠️  Spatial coordinates not available")
            return
            
        if color_by not in self.adata.obs.columns:
            print(f"⚠️  Column '{color_by}' not found in adata.obs")
            return
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Plot spatial with annotations
        sc.pl.embedding(self.adata, basis='spatial', color=color_by, 
                       ax=ax, show=False, legend_loc='right margin', 
                       frameon=False, size=point_size)
        
        plt.title(f'Spatial: {color_by.replace("_", " ").title()}', fontsize=14)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved spatial plot to {save_path}")
        
        plt.show()
    
    def plot_comprehensive_annotation_results(self, figsize: Tuple[int, int] = (20, 16),
                                           save_path: Optional[str] = None) -> None:
        """
        Create comprehensive visualization of annotation results.
        """
        if 'cell_type_assignment' not in self.adata.obs.columns:
            print("⚠️  No cell type assignments found. Run assign_cell_types() first.")
            return
            
        # Create subplot grid
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(3, 4, hspace=0.4, wspace=0.3)
        
        # 1. UMAP with assignments
        if 'X_umap' in self.adata.obsm:
            ax1 = fig.add_subplot(gs[0, 0])
            sc.pl.umap(self.adata, color='cell_type_assignment', ax=ax1, 
                      show=False, frameon=False, legend_loc=None)
            ax1.set_title('UMAP: Cell Types')
        
        # 2. Spatial with assignments
        if 'spatial' in self.adata.obsm:
            ax2 = fig.add_subplot(gs[0, 1])
            sc.pl.embedding(self.adata, basis='spatial', color='cell_type_assignment',
                           ax=ax2, show=False, frameon=False, legend_loc=None)
            ax2.set_title('Spatial: Cell Types')
        
        # 3. Assignment confidence
        if 'assignment_confidence' in self.adata.obs.columns:
            ax3 = fig.add_subplot(gs[0, 2])
            if 'X_umap' in self.adata.obsm:
                sc.pl.umap(self.adata, color='assignment_confidence', ax=ax3,
                          show=False, frameon=False, legend_loc=None)
                ax3.set_title('UMAP: Confidence')
        
        # 4. Cell type proportions
        ax4 = fig.add_subplot(gs[0, 3])
        assignment_counts = self.adata.obs['cell_type_assignment'].value_counts()
        wedges, texts, autotexts = ax4.pie(assignment_counts.values, 
                                          labels=assignment_counts.index,
                                          autopct='%1.1f%%', startangle=90)
        ax4.set_title('Cell Type Proportions')
        
        # 5-8. Top scoring genes for each major cell type (first 4)
        score_cols = [col for col in self.adata.obs.columns if col.endswith('_score')][:4]
        
        for i, score_col in enumerate(score_cols):
            ax = fig.add_subplot(gs[1, i])
            cell_type = score_col.replace('_score', '')
            
            if 'spatial' in self.adata.obsm:
                sc.pl.embedding(self.adata, basis='spatial', color=score_col,
                               ax=ax, show=False, frameon=False, legend_loc=None)
                ax.set_title(f'{cell_type} Score')
        
        # 9. Score distribution
        ax9 = fig.add_subplot(gs[2, 0:2])
        score_data = []
        for score_col in score_cols:
            cell_type = score_col.replace('_score', '')
            scores = self.adata.obs[score_col]
            score_data.extend([(cell_type, score) for score in scores])
        
        if score_data:
            score_df = pd.DataFrame(score_data, columns=['Cell_Type', 'Score'])
            sns.boxplot(data=score_df, x='Cell_Type', y='Score', ax=ax9)
            ax9.set_title('Score Distributions by Cell Type')
            ax9.tick_params(axis='x', rotation=45)
        
        # 10. Assignment statistics
        ax10 = fig.add_subplot(gs[2, 2:])
        ax10.axis('off')
        
        # Create summary text
        summary_text = ["📊 Annotation Summary:\n"]
        summary_text.append(f"Total cells: {self.adata.n_obs:,}")
        summary_text.append(f"Cell types identified: {len(assignment_counts)}")
        
        if 'assignment_confidence' in self.adata.obs.columns:
            mean_conf = self.adata.obs['assignment_confidence'].mean()
            summary_text.append(f"Mean confidence: {mean_conf:.3f}")
        
        summary_text.append(f"\nCell type breakdown:")
        for cell_type, count in assignment_counts.head(10).items():
            pct = count / len(self.adata) * 100
            summary_text.append(f"  {cell_type}: {count:,} ({pct:.1f}%)")
        
        ax10.text(0.05, 0.95, '\n'.join(summary_text), 
                 transform=ax10.transAxes, fontsize=10, 
                 verticalalignment='top', family='monospace')
        
        plt.suptitle('Comprehensive Cell Type Annotation Results', fontsize=16, y=0.95)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved comprehensive results to {save_path}")
        
        plt.show()
    
    def validate_spatial_assignments(self) -> Dict[str, float]:
        """
        Validate spatial coherence of cell type assignments.
        
        Returns
        -------
        dict
            Validation metrics
        """
        if 'cell_type_assignment' not in self.adata.obs.columns:
            print("⚠️  No assignments found")
            return {}
            
        if 'spatial' not in self.adata.obsm:
            print("⚠️  No spatial coordinates found")
            return {}
        
        # TODO: Implement spatial coherence metrics
        # - Local purity scores
        # - Neighborhood composition analysis
        # - Spatial autocorrelation
        
        print("🔍 Spatial validation not yet implemented")
        return {}
    
    def save_annotations(self, output_path: str) -> None:
        """
        Save annotated data to file.
        
        Parameters
        ----------
        output_path : str
            Path to save annotated h5ad file
        """
        # Add annotation metadata
        self.adata.uns['annotation_results'] = self.annotation_results
        self.adata.uns['marker_genes'] = self.marker_genes
        
        # Save
        self.adata.write_h5ad(output_path)
        print(f"✓ Saved annotated data to {output_path}")
        
    def get_annotation_summary(self) -> pd.DataFrame:
        """
        Get summary table of annotation results.
        
        Returns
        -------
        pandas.DataFrame
            Summary of cell type assignments
        """
        if 'cell_type_assignment' not in self.adata.obs.columns:
            print("⚠️  No assignments found")
            return pd.DataFrame()
        
        # Create summary
        assignment_counts = self.adata.obs['cell_type_assignment'].value_counts()
        summary_df = pd.DataFrame({
            'Cell_Type': assignment_counts.index,
            'Count': assignment_counts.values,
            'Percentage': (assignment_counts.values / self.adata.n_obs) * 100
        })
        
        # Add confidence stats if available
        if 'assignment_confidence' in self.adata.obs.columns:
            conf_stats = self.adata.obs.groupby('cell_type_assignment')['assignment_confidence'].agg(['mean', 'std'])
            summary_df['Mean_Confidence'] = summary_df['Cell_Type'].map(conf_stats['mean'])
            summary_df['Std_Confidence'] = summary_df['Cell_Type'].map(conf_stats['std'])
        
        return summary_df.sort_values('Count', ascending=False)


def create_marker_gene_template(output_path: str) -> None:
    """
    Create a template CSV file for marker genes.
    
    Parameters
    ----------
    output_path : str
        Path to save the template CSV file
    """
    template_data = {
        'cell_type': [
            'Neurons', 'Glia', 'Hemocytes', 'sKC', 'lKC', 'mKC', 'OLC', 'PN',
            'VisualSystem', 'MotorNeurons', 'Interneurons'
        ],
        'marker_genes': [
            'elav,nSyb,Syt1,CaMKII,futsch,brp',
            'repo,Gs2,Eaat1,wrapper,nrv2,Gli',
            'Hml,srp,gcm,lz,Pxn,NimC1',
            'mb247,FasII,arm,Lac,Synapsin',
            'mb247,FasII,VGlut,ChAT,Gad1',
            'mb247,FasII,Synapsin,VGlut,ChAT',
            'Orco,Ir8a,Ir25a,Ir76b,Or22a',
            'GH146,Mz19,NP225,NP2631,NP3056',
            'rh1,rh3,rh4,rh5,rh6,norpA',
            'ChAT,VGlut,Gad1,OK371,eve',
            'Gad1,VGlut,ChAT,DvGlut,Rdl'
        ],
        'description': [
            'Pan-neuronal markers',
            'Glial cell markers', 
            'Immune cell markers',
            'Small Kenyon cells',
            'Large Kenyon cells',
            'Medium Kenyon cells',
            'Olfactory lobe cells',
            'Projection neurons',
            'Visual system cells',
            'Motor neurons',
            'Interneurons'
        ]
    }
    
    df = pd.DataFrame(template_data)
    df.to_csv(output_path, index=False)
    print(f"✓ Created marker gene template at {output_path}")
    print("Edit this file with your literature-curated marker genes!")


if __name__ == "__main__":
    # Example usage
    print("🧠 Bee Brain Cell Type Annotation Toolkit")
    print("This module provides comprehensive cell type annotation capabilities.")
    print("Use the BeebrainCellTypeAnnotator class for your analysis.")
    
    # Create template if run directly
    template_path = "marker_genes_template.csv"
    create_marker_gene_template(template_path) 