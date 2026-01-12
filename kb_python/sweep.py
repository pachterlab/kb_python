import os
from typing import Optional
import anndata as ad
from .utils import (
    import_matrix_as_anndata,
)

from .logging import logger

@logger.namespaced("sweep")
def sweep(
    kb_count_dir: str | ad.AnnData,
    adata_out: Optional[str] = "adata_denoised.h5ad",
    max_iter: int = 500,
    init_alpha: float = 0.9,
    beta: float = 0.1,
    eps: float = 1e-12,
    log_eps: float = 1e-300,
    dirichlet_lambda: Optional[float] = 500,
    integer_out: bool = False,
    threads: int = 1,
    fixed_celltype: bool = True,
    freeze_empty: bool = True,
    freeze_ambient_profile: bool = True,
    empty_droplet_method: str = "threshold",
    ambient_threshold: Optional[float] = 0.0,
    umi_cutoff: Optional[int] = None,
    expected_cells: Optional[int] = None,
    tol: float = 1e-3,
    min_tol: float = 1e-6,
    leiden_resolution: float = 1.0,
    random_state: Optional[int] = 42,
    verbose: int = 0,
    quiet: bool = False,
    log_file: Optional[str] = None,
    debug: bool = False,
):
    """
    Wraps cellsweep

    Denoise a count matrix using an Expectation-Maximization (EM) algorithm that
    models each observed count as a mixture of ambient RNA, bulk RNA, and true cell-type 
    signal.

    This function optionally operates on real cells only (excluding identified empty droplets),
    fixing the ambient expression profile and optionally fixing cell-type assignments.
    It iteratively estimates latent variables representing per-cell ambient fractions
    (alpha_i), a bulk contamination factor (beta), per-cell-type expression profiles 
    (p_k), and an ambient contamination profile (a) until convergence.

    Parameters
    ----------
    kb_count_dir : str | AnnData
        Path to one of the following:
        1. Output of kb count
        2. An AnnData object containing the unfiltered count matrix

    adata_out : str, default "adata_straightened.h5ad"
        Path to write the denoised AnnData object (must end with `.h5ad`).

    max_iter : int, default 500
        Maximum number of EM iterations.

    init_alpha : float, default 0.9
       Initial value of alpha_n for each cell. Works better when set to a higher number than expected (expected is around 0.05 per cell).
    
    beta : float, default 0.1
        Initial beta (percent bulk contamination) value for each cell. Works better when set to a higher number than expected (expected is around 0.05). 
        Set to a lower value than alpha_init since bulk contamination is usually less than ambient contamination.

    eps : float, default 1e-12
        Numerical stability constant to prevent division by zero.

    log_eps : float, default 1e-300
        Numerical stability constant to log(0).

    dirichlet_lambda: float, default 10
        Pseudocount. Will be divided by the number of genes G. Higher values lead to smoother cell-type profiles.

    integer_out : bool, default False
        If True, rounds denoised counts to nearest integer before saving.

    threads : int, default 1
        number of numba threads

    fixed_celltype : bool, default False
        If True, keeps cell-type assignments fixed during EM updates.

    freeze_empty : bool, default True
        If True, does not attempt to reestimate the percent contamination of empty droplets from 100%

    freeze_ambient_profile: bool, default True
        If True, does not update the ambient profile (a) based upon alpha

    empty_droplet_method : str, default "threshold"
        Strategy to infer empty droplets if `is_empty` is not present.
        Options may include "threshold", "quantile", or model-based approaches.
    
    ambient_threshold : float | None, default 0.0
        Optional ambient RNA fraction threshold for classifying droplets as empty.

    umi_cutoff : int | None, default None
        Optional absolute UMI count threshold for classifying droplets as empty.

    expected_cells : int | None, default None
        Expected number of real cells, used when estimating thresholds.

    tol: float, default 1e-3
        The relative change in likelihood below which training is discontinued
    
    min_tol: float, default 1e-6
        The minimum absolute change in likelihood below which training is discontinued.
    
    leiden_resolution : float, default 1.0
        Resolution parameter for Leiden clustering.

    random_state: int | None, default 42
        Random seed

    verbose : int, default 0
        Verbosity level (2 debug, 1 info, 0 warning, -1 error, -2 critical).

    quiet : bool, default False
        Suppresses most log output when True.

    log_file : str | None, default None
        Optional path to save EM iteration logs.

    Returns
    -------
    AnnData
        Denoised AnnData object with updated `adata.X`, and
        added fields:
        - `adata.layers["raw"]` : raw count matrix
        - `adata.obs["cell_ambient_fraction"]` : estimated ambient fraction per cell
        - `adata.uns["em_convergence"]` : diagnostics and log-likelihood trace
        - `adata.obs["alpha_hat"]' : final optimized alpha values
        - `adata.obs["z_hat"]` : final cell-type assignments (These should not change)
        - `adata.uns["p_hat"]` : final optimized matrix of cell-type profiles (K x G)
        - `adata.uns["beta_hat"]` : final optimized beta
        - `adata.var["ambient_hat"]` : final optimized ambient distribution
        - `adata.uns["loglike"]` : final log-likelihood (note that this value is not the 
           complete log-likelihood, only the relative log-likelihood)

    Notes
    -----
    The EM algorithm proceeds by:
      1. E-step: Update expected value of true, ambient noise, and bulk noise counts for each cell and gene.
      2. M-step: Update parameters (alpha, beta, p_k, a).
      3. Iterate until convergence (relative change in ll < `tol`) or reaching `max_iter`.
    """
    import cellsweep

    #* load adata
    logger.info("Loading count matrix into AnnData object...")
    if isinstance(kb_count_dir, str):
        if os.path.isdir(kb_count_dir):
            if not os.path.exists(os.path.join(kb_count_dir, "counts_unfiltered")):
                raise ValueError(f"Provided kb_count_dir path {kb_count_dir} does not contain 'counts_unfiltered' directory.")
            matrix_path = os.path.join(kb_count_dir, "counts_unfiltered", "cells_x_genes.mtx")
            barcodes_path = os.path.join(kb_count_dir, "counts_unfiltered", "cells_x_genes.barcodes.txt")
            genes_path = os.path.join(kb_count_dir, "counts_unfiltered", "cells_x_genes.genes.names.txt")
            adata = import_matrix_as_anndata(matrix_path, barcodes_path, genes_path)
            adata = cellsweep.utils.read_kb_mtx_as_adata(kb_count_dir)
        elif os.path.isfile(kb_count_dir):
            if not kb_count_dir.endswith(".h5ad"):
                raise ValueError(f"Provided kb_count_dir file {kb_count_dir} is not an .h5ad file.")
            adata = ad.read_h5ad(kb_count_dir)
        else:
            raise ValueError(f"Provided kb_count_dir path {kb_count_dir} is neither a directory nor a file.")
    
    #* add leiden clusters as celltypes
    logger.info("Preprocessing and clustering with Scanpy to assign cell types...")
    adata_processed_tmp = cellsweep.utils.run_scanpy_preprocessing_and_clustering(
        adata=adata,
        filter_empty_droplets=False,
        min_genes=None,
        min_counts=None,
        min_cells=None,
        umi_top_percentile_to_remove=None,
        unique_genes_top_percentile_to_remove=None,
        mt_gene_percentile_to_remove=None,
        max_mt_percentage=None,
        n_top_genes=2000,
        hvg_flavor="seurat_v3",
        n_pcs=50,
        n_neighbors=15,
        leiden_resolution=leiden_resolution,
        seed=random_state,
        verbose=verbose,
        quiet=quiet
    )
    adata.obs["celltype"] = adata_processed_tmp.obs["leiden"].reindex(adata.obs.index)
    del adata_processed_tmp

    #* run cellsweep
    logger.info("Running CellSweep denoising...")
    adata_cellsweep = cellsweep.denoise_count_matrix(
        adata=adata,
        adata_out=adata_out,
        max_iter=max_iter,
        init_alpha=init_alpha,
        beta=beta,
        eps=eps,
        log_eps=log_eps,
        dirichlet_lambda=dirichlet_lambda,
        integer_out=integer_out,
        threads=threads,
        fixed_celltype=fixed_celltype,
        freeze_empty=freeze_empty,
        freeze_ambient_profile=freeze_ambient_profile,
        empty_droplet_method=empty_droplet_method,
        ambient_threshold=ambient_threshold,
        umi_cutoff=umi_cutoff,
        expected_cells=expected_cells,
        tol=tol,
        min_tol=min_tol,
        random_state=random_state,
        verbose=verbose,
        quiet=quiet,
        log_file=log_file,
        debug=debug,
    )

    logger.info("CellSweep denoising complete.")
    return adata_cellsweep

