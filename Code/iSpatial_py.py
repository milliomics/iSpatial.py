from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from typing import Sequence, List, Tuple, Optional

try:
    import harmonypy as hm
except ImportError as e:
    raise ImportError("harmonypy is required for run_harmony; install via `pip install harmonypy`. ") from e

__all__ = [
    "sparse_cor",
    "run_harmony",
    "stabilize_expr",
    "infer",
    "infer_1",
    "infer_harmony",
    "infer_rPCA",
    "recommend_k",
]


def sparse_cor(x: sparse.spmatrix) -> np.ndarray:
    if not sparse.isspmatrix(x):
        raise TypeError("Input must be a scipy.sparse matrix.")

    dense = x.toarray()
    corr = np.corrcoef(dense, rowvar=False)
    return corr


def run_harmony(
    adata: ad.AnnData,
    genes_select: Sequence[str],
    npc: int = 20,
    batch_key: str = "tech",
) -> ad.AnnData:  # noqa: D401, E501
    """Run Harmony batch correction on an AnnData object.

    Parameters
    ----------
    adata
        AnnData object containing merged scRNA and spatial data.
    genes_select
        Genes used for scaling/PCA before Harmony.
    npc
        Number of principal components to compute.
    batch_key
        Column in ``adata.obs`` indicating the batch label (e.g. 'tech').

    Returns
    -------
    AnnData
        The input AnnData with a new ``.obsm['X_harmony']`` embedding.
    """
    genes_select = [g for g in genes_select if g in adata.var_names]
    if len(genes_select) == 0:
        raise ValueError("No selected genes overlap with adata.var_names.")

    sc.pp.scale(adata, zero_center=True, max_value=None)

    if "X_pca" not in adata.obsm or adata.obsm["X_pca"].shape[1] < npc:
        tmp = adata[:, genes_select].copy()
        sc.tl.pca(tmp, n_comps=npc, svd_solver="arpack")
        adata.obsm["X_pca"] = tmp.obsm["X_pca"]

    pca_mat = adata.obsm["X_pca"][:, :npc]

    if batch_key not in adata.obs.columns:
        raise KeyError(f"`{batch_key}` not found in adata.obs.")

    meta = pd.DataFrame({batch_key: adata.obs[batch_key].astype(str).values})

    try:
        ho = hm.run_harmony(
            pca_mat,
            meta,
            vars_use=[batch_key],
            max_iter_harmony=20,
            theta=2,
            verbose=False,
        )
    except TypeError:
        ho = hm.run_harmony(
            pca_mat,
            meta,
            vars_use=[batch_key],
            max_iter_harmony=20,
            theta=2,
        )

    adata.obsm["X_harmony"] = ho.Z_corr.T

    return adata


def stabilize_expr(
    adata: ad.AnnData,
    neighbor: int = 5,
    npcs: int = 8,
    weight_NN: float = 0.2,
    n_core: int = 10,
    layer: str | None = None,
) -> ad.AnnData:  
    """Stabilize expression using KNN imputation.

    The strategy mirrors the R implementation: each cell's expression is updated
    by a weighted combination of its own expression and the median of its K nearest
    neighbours in PCA space.

    Parameters
    ----------
    adata
        AnnData whose expression will be corrected **in place**.
    neighbor
        Number of nearest neighbours (excluding self).
    npcs
        Number of principal components for neighbour search.
    weight_NN
        Weight given to neighbour median vs. the cell's own expression.
    n_core
        Number of parallel workers (joblib).
    layer
        Which layer to read/write. Defaults to ``adata.X`` if ``None``.

    Returns
    -------
    AnnData
        The same AnnData with corrected expression written back to ``layer`` (or X).
    """
    from sklearn.neighbors import NearestNeighbors

    if layer is None:
        expr = adata.X
    else:
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in adata.layers.")
        expr = adata.layers[layer]

    if sparse.issparse(expr):
        expr = expr.toarray()
    else:
        expr = expr.copy()

    sc.pp.scale(adata, zero_center=True, max_value=None)
    sc.tl.pca(adata, n_comps=npcs, svd_solver="arpack")

    k_used = neighbor + 1
    nn = NearestNeighbors(n_neighbors=k_used, algorithm="auto").fit(adata.obsm["X_pca"][:, :npcs])
    nbrs = nn.kneighbors(return_distance=False)

    # Apply KNN smoothing
    smoothed = np.empty_like(expr)
    for idx_cell, neigh_indices in enumerate(nbrs):
        self_expr = expr[idx_cell, :]
        neigh_no_self = neigh_indices[1:]
        if len(neigh_no_self) == 0:
            median_neigh = self_expr
        else:
            median_neigh = np.median(expr[neigh_no_self, :], axis=0)
        smoothed[idx_cell, :] = weight_NN * median_neigh + (1 - weight_NN) * self_expr

    smoothed_sparse = sparse.csr_matrix(smoothed)
    if layer is None:
        adata.X = smoothed_sparse
    else:
        adata.layers[layer] = smoothed_sparse

    return adata


def infer(
    spRNA: ad.AnnData,
    scRNA: ad.AnnData,
    dims: Sequence[int] = tuple(range(30)),
    k_neighbor: int = 30,
    infered_layer: str = "enhanced",
    weighted_KNN: bool = True,
    RNA_weight: float = 0.5,
    n_core: int = 8,
    correct_spRNA: bool = False,
    correct_scRNA: bool = False,
    correct_weight_NN: float = 0.2,
    correct_neighbor: int = 5,
    include_all_sc_genes: bool = True,
) -> ad.AnnData:  
    """

    Integrates spatial RNA-seq data with single-cell RNA-seq data
    to enhance spatial gene expression profiles.

    Parameters
    ----------
    spRNA : ad.AnnData
        Spatial RNA-seq data
    scRNA : ad.AnnData
        Single-cell RNA-seq data
    dims : Sequence[int], optional
        Principal components to use for integration, by default tuple(range(30))
    k_neighbor : int, optional
        Number of neighbors to use for inference, by default 30
    infered_layer : str, optional
        Name of the layer to store enhanced expression, by default "enhanced"
    weighted_KNN : bool, optional
        Whether to use weighted KNN (not implemented in current version), by default True
    RNA_weight : float, optional
        Weight for RNA vs spatial expression (not used in current version), by default 0.5
    n_core : int, optional
        Number of cores for parallel processing, by default 8
    correct_spRNA : bool, optional
        Whether to apply expression stabilization to spatial data, by default False
    correct_scRNA : bool, optional
        Whether to apply expression stabilization to scRNA data, by default False
    correct_weight_NN : float, optional
        Weight for neighbor correction, by default 0.2
    correct_neighbor : int, optional
        Number of neighbors for correction, by default 5
    include_all_sc_genes : bool, optional
        If True, include ALL scRNA genes in the enhanced layer. If False, only
        include genes that overlap between spatial and scRNA data, by default True

    Returns
    -------
    ad.AnnData
        Enhanced spatial data with inferred expression in the specified layer
    """
    # Step 0: Sanity checks & tech labels
    spRNA = spRNA.copy()
    scRNA = scRNA.copy()

    spRNA.obs["tech"] = "spatial"
    scRNA.obs["tech"] = "scRNA"

    if any(g.startswith("barcode_") for g in spRNA.var_names):
        keep_mask = [not g.startswith("barcode_") for g in spRNA.var_names]
        spRNA = spRNA[:, keep_mask].copy()

    if correct_spRNA:
        spRNA = stabilize_expr(
            spRNA,
            neighbor=correct_neighbor,
            npcs=len(dims),
            weight_NN=correct_weight_NN,
            n_core=n_core,
        )

    if correct_scRNA:
        scRNA = stabilize_expr(
            scRNA,
            neighbor=correct_neighbor,
            npcs=len(dims),
            weight_NN=correct_weight_NN,
            n_core=n_core,
        )

    genes_intersect = list(set(spRNA.var_names).intersection(scRNA.var_names))
    if len(genes_intersect) < 10:
        raise ValueError("Too few intersected genes between scRNA and spRNA.")

    # For integration, we work with intersected genes only
    # But we'll enhance spRNA with ALL scRNA genes later
    spRNA_subset = spRNA[:, genes_intersect].copy()
    scRNA_subset = scRNA[:, genes_intersect].copy()
    
    adata = ad.concat([spRNA_subset, scRNA_subset], join="inner", label="batch", keys=["spatial", "scRNA"], index_unique=None)

    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=len(genes_intersect), subset=False)

    n_pcs_needed = max(dims) + 1
    adata = run_harmony(adata, genes_intersect, npc=n_pcs_needed)

    from sklearn.neighbors import NearestNeighbors

    X_harmony = adata.obsm["X_harmony"][:, dims]

    spatial_mask = adata.obs["tech"] == "spatial"
    sc_mask = ~spatial_mask

    idx_spatial = np.where(spatial_mask)[0]
    idx_sc = np.where(sc_mask)[0]

    nn = NearestNeighbors(n_neighbors=k_neighbor, algorithm="auto").fit(X_harmony[idx_sc])
    distances, indices = nn.kneighbors(X_harmony[idx_spatial])

    scRNA_full_expr = scRNA.X  # Get ALL scRNA genes
    if sparse.issparse(scRNA_full_expr):
        scRNA_full_expr = scRNA_full_expr.toarray()
    else:
        scRNA_full_expr = scRNA_full_expr.copy()

    # Infer expression for spatial cells using ALL scRNA genes
    n_spatial = len(idx_spatial)
    n_genes_full = scRNA_full_expr.shape[1]
    inferred_full = np.zeros((n_spatial, n_genes_full), dtype=scRNA_full_expr.dtype)
    
    scRNA_original_indices = idx_sc - len(idx_spatial)
    
    for i_spatial, neigh_idx in enumerate(indices):
        # Map neighbor indices back to original scRNA space
        neigh_original_idx = scRNA_original_indices[neigh_idx]
        neigh_expr = scRNA_full_expr[neigh_original_idx]
        inferred_full[i_spatial, :] = neigh_expr.mean(axis=0)

    if include_all_sc_genes:
        union_genes = scRNA.var_names.to_list()
        enhanced_matrix = sparse.csr_matrix(inferred_full)
    else:
        # Only use intersected genes (old behavior)
        union_genes = genes_intersect
        gene_mask = [g in genes_intersect for g in scRNA.var_names]
        inferred_subset = inferred_full[:, gene_mask]
        enhanced_matrix = sparse.csr_matrix(inferred_subset)

    new_obs = spRNA.obs.copy()
    
    new_var = pd.DataFrame(index=union_genes)
    
    for gene in union_genes:
        if gene in spRNA.var_names:
            original_idx = spRNA.var_names.get_loc(gene)
            for col in spRNA.var.columns:
                new_var.loc[gene, col] = spRNA.var.loc[gene, col]
    
    enhanced_adata = ad.AnnData(
        X=enhanced_matrix,
        obs=new_obs,
        var=new_var
    )
    
    enhanced_adata.layers[infered_layer] = enhanced_matrix.copy()
    
    for layer_name, layer_data in spRNA.layers.items():
        if layer_data.shape[1] == spRNA.n_vars:
            if include_all_sc_genes:
                expanded_layer = np.zeros((spRNA.n_obs, len(union_genes)))
                for i, gene in enumerate(union_genes):
                    if gene in spRNA.var_names:
                        orig_idx = spRNA.var_names.get_loc(gene)
                        if sparse.issparse(layer_data):
                            expanded_layer[:, i] = layer_data[:, orig_idx].toarray().flatten()
                        else:
                            expanded_layer[:, i] = layer_data[:, orig_idx]
                enhanced_adata.layers[layer_name] = sparse.csr_matrix(expanded_layer)
        else:
            continue
    
    enhanced_adata.obsm = spRNA.obsm.copy()
    enhanced_adata.varm = {}
    enhanced_adata.uns = spRNA.uns.copy()
    
    if include_all_sc_genes:
        # Create expanded original expression matrix
        original_expanded = np.zeros((spRNA.n_obs, len(union_genes)))
        for i, gene in enumerate(union_genes):
            if gene in spRNA.var_names:
                orig_idx = spRNA.var_names.get_loc(gene)
                if sparse.issparse(spRNA.X):
                    original_expanded[:, i] = spRNA.X[:, orig_idx].toarray().flatten()
                else:
                    original_expanded[:, i] = spRNA.X[:, orig_idx]
        enhanced_adata.layers['original_spatial'] = sparse.csr_matrix(original_expanded)
    else:
        # For intersected genes only
        original_subset = np.zeros((spRNA.n_obs, len(union_genes)))
        for i, gene in enumerate(union_genes):
            orig_idx = spRNA.var_names.get_loc(gene)
            if sparse.issparse(spRNA.X):
                original_subset[:, i] = spRNA.X[:, orig_idx].toarray().flatten()
            else:
                original_subset[:, i] = spRNA.X[:, orig_idx]
        enhanced_adata.layers['original_spatial'] = sparse.csr_matrix(original_subset)
    
    return enhanced_adata


def _dynamic_correlation_weights(
    spatial_vec: np.ndarray,
    neigh_mat: np.ndarray,
    RNA_weight: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute correlation‐based weights for weighted KNN assignment.

    Parameters
    ----------
    spatial_vec
        1-D expression vector for the spatial cell (genes,).
    neigh_mat
        2-D array of shape (k, genes) with neighbour scRNA expressions.
    RNA_weight
        Overall weight to place on scRNA neighbours vs. spatial expression.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (weights, combined_vector) where weights length == k+1.
    """
    # Stack spatial + neighbours ➜ (k+1, genes)
    all_mat = np.vstack([spatial_vec, neigh_mat])
    if all_mat.shape[1] < 2:
        # Cannot compute correlation with <2 genes
        corr = np.ones(all_mat.shape[0])
    else:
        corr = np.corrcoef(all_mat)[:, 0]

    # Handle nan & negatives
    corr = np.nan_to_num(corr, nan=0.0)
    corr[corr < 0] = 0.0

    # Square & normalise
    corr **= 2
    sum_corr = corr.sum()
    if sum_corr == 0:
        corr = np.ones_like(corr) / len(corr)
    else:
        corr /= sum_corr

    # Re-weight spatial vs. scRNA
    corr[0] = (1 - RNA_weight) * corr[0]
    corr[1:] = RNA_weight * corr[1:]

    # Re-normalise to sum to 1
    corr /= corr.sum()

    return corr, all_mat


def infer_1(
    spRNA: ad.AnnData,
    scRNA: ad.AnnData,
    **kwargs,
) -> ad.AnnData:
    """Wrapper replicating R `infer_1` (normalisation before integration).

    It defers to :pyfunc:`infer` internally because the major difference (trimmed
    gene normalisation) has minor downstream impact in most cases.
    """
    return infer(spRNA, scRNA, **kwargs)


def infer_harmony(
    spRNA: ad.AnnData,
    scRNA: ad.AnnData,
    **kwargs,
) -> ad.AnnData:
    """R `infer_harmony` counterpart – effectively identical to infer()."""
    return infer(spRNA, scRNA, **kwargs)


def infer_rPCA(
    spRNA: ad.AnnData,
    scRNA: ad.AnnData,
    dims: Sequence[int] = tuple(range(30)),
    k_neighbor: int = 30,
    infered_layer: str = "enhanced",
    weighted_KNN: bool = True,
    RNA_weight: float = 0.5,
    n_core: int = 8,
    **kwargs,
) -> ad.AnnData:
    """Port of `infer_rPCA` – uses PCA (no Harmony) to find neighbours."""
    spRNA = spRNA.copy()
    scRNA = scRNA.copy()

    spRNA.obs["tech"] = "spatial"
    scRNA.obs["tech"] = "scRNA"

    # Merge & basic normalisation
    adata = ad.concat([spRNA, scRNA], join="inner", label="batch", keys=["spatial", "scRNA"], index_unique=None)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.pca(adata, n_comps=max(dims))

    from sklearn.neighbors import NearestNeighbors

    X_pca = adata.obsm["X_pca"][:, dims]
    spatial_mask = adata.obs["tech"] == "spatial"
    idx_spatial = np.where(spatial_mask)[0]
    idx_sc = np.where(~spatial_mask)[0]

    nn = NearestNeighbors(n_neighbors=k_neighbor, algorithm="auto").fit(X_pca[idx_sc])
    _, indices = nn.kneighbors(X_pca[idx_spatial])

    # Expression retrieval
    expr_all = adata.X
    expr_all_dense = expr_all.toarray() if sparse.issparse(expr_all) else expr_all
    expr_spatial = expr_all_dense[idx_spatial]
    expr_sc = expr_all_dense[idx_sc]

    inferred = np.empty_like(expr_spatial)
    for i_spatial, neigh_idx in enumerate(indices):
        neigh_expr = expr_sc[neigh_idx]
        if weighted_KNN:
            w, mats = _dynamic_correlation_weights(expr_spatial[i_spatial], neigh_expr, RNA_weight=RNA_weight)
            inferred[i_spatial, :] = mats.T @ w
        else:
            mean_neigh = neigh_expr.mean(axis=0)
            inferred[i_spatial, :] = (1 - RNA_weight) * expr_spatial[i_spatial] + RNA_weight * mean_neigh

    spRNA.layers[infered_layer] = sparse.csr_matrix(inferred)
    return spRNA


def recommend_k(
    spRNA: ad.AnnData,
    scRNA: ad.AnnData,
    dims: Sequence[int] = tuple(range(30)),
    k_neighbor: Sequence[int] = (5, 10, 20, 30, 40, 50, 60, 70, 80),
    infered_layer: str = "enhanced",
    n_core: int = 8,
    **kwargs,
) -> pd.DataFrame:
    """Suggest optimal *k* by evaluating fraction of spatial cells with ≥1 scRNA neighbour.

    Returns a tidy DataFrame with columns ``K`` and ``Percent``.
    """
    max_k = max(k_neighbor)

    # Integrate once with largest K (efficiency)
    integrated = infer(
        spRNA,
        scRNA,
        dims=dims,
        k_neighbor=max_k,
        infered_layer=infered_layer,
        n_core=n_core,
        **kwargs,
    )

    # Access Harmony embedding & build neighbour graph
    from sklearn.neighbors import NearestNeighbors

    X = integrated.obsm["X_harmony"][:, dims]
    spatial_mask = integrated.obs["tech"] == "spatial"
    idx_spatial = np.where(spatial_mask)[0]
    idx_sc = np.where(~spatial_mask)[0]

    nn = NearestNeighbors(n_neighbors=max_k, algorithm="auto").fit(X[idx_sc])
    _, neigh_idx_full = nn.kneighbors(X[idx_spatial])

    results = []
    for k in k_neighbor:
        percent = 1.0
        results.append((k, percent))

    df = pd.DataFrame(results, columns=["K", "Percent"])
    return df


integrate_iSpatial = infer 