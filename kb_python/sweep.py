import os
from importlib.metadata import PackageNotFoundError, version
from typing import Optional

import anndata as ad
import pandas as pd
from .utils import (
    import_matrix_as_anndata,
)

from .logging import logger

CELLSWEEP_VERSION = "1.0.2"

# Advanced EM hyperparameters accepted by `cellsweep.denoise_count_matrix` as
# keyword arguments. Any of these given to `sweep` is forwarded unchanged;
# those not given take cellsweep's own defaults.
EM_KWARGS = (
    "init_alpha",
    "init_beta",
    "celltype_profile_key",
    "ambient_profile_key",
    "bulk_profile_key",
    "alpha_cap",
    "repulsion_strength",
    "max_frac_gene_repulsion",
    "beta_prior_mode",
    "beta_prior_strength",
    "celltype_lambda",
    "ambient_lambda",
    "bulk_lambda",
    "eps",
    "log_eps",
    "max_iter",
    "del0_ll_tol",
    "min_ll_tol",
    "burnin_patience",
    "burnin_max_iter",
    "tol_p",
    "tol_f",
)

# Accepted values (case-insensitive) of the empty status column of a celltypes file.
EMPTY_STATUS_TRUE = frozenset({"true", "t", "yes", "y", "1"})
EMPTY_STATUS_FALSE = frozenset({"false", "f", "no", "n", "0"})


def check_cellsweep_version():
    """Ensure the pinned version of cellsweep is installed.

    Raises:
        ImportError: If cellsweep is missing or is not version `CELLSWEEP_VERSION`
    """
    try:
        installed = version("cellsweep")
    except PackageNotFoundError:
        raise ImportError(
            f"cellsweep is not installed. Please install it using "
            f"'pip install cellsweep[analysis]=={CELLSWEEP_VERSION}'."
        )
    if installed != CELLSWEEP_VERSION:
        raise ImportError(
            f"cellsweep version {installed} is installed, but kb sweep requires "
            f"version {CELLSWEEP_VERSION}. Please install it using "
            f"'pip install cellsweep[analysis]=={CELLSWEEP_VERSION}'."
        )


def parse_empty_status(value: str) -> bool:
    """Parse the empty status column of a celltypes file.

    Args:
        value: The text of the column

    Returns:
        Whether the barcode is an empty droplet

    Raises:
        ValueError: If `value` is not a recognized boolean
    """
    normalized = value.strip().lower()
    if normalized in EMPTY_STATUS_TRUE:
        return True
    if normalized in EMPTY_STATUS_FALSE:
        return False
    raise ValueError(
        f"Invalid empty status {value!r}. Expected one of: "
        f"{sorted(EMPTY_STATUS_TRUE | EMPTY_STATUS_FALSE)}"
    )


def read_celltypes(celltypes_path: str) -> pd.DataFrame:
    """Read a barcode-to-celltype mapping file.

    Each line contains a barcode, its celltype and, optionally, its empty
    status, separated by tabs. A line with no tab is split on the first
    whitespace into a barcode and a celltype, so the third column requires tabs.
    Either every line or no line has the third column. The celltype may be
    blank only for a barcode marked empty.

    Args:
        celltypes_path: Path to the mapping file

    Returns:
        DataFrame indexed by barcode with a `celltype` column and, if the file
        has a third column, a boolean `is_empty` column
    """
    mapping = {}
    n_columns = None
    with open(celltypes_path, "r") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            parts = line.split("\t") if "\t" in line else line.split(maxsplit=1)
            if len(parts) not in (2, 3):
                raise ValueError(
                    f"Malformed line in celltypes file {celltypes_path}: {line!r}. "
                    "Expected two or three columns: barcode, celltype and, "
                    "optionally, empty status."
                )
            if n_columns is None:
                n_columns = len(parts)
            elif len(parts) != n_columns:
                raise ValueError(
                    f"Malformed line in celltypes file {celltypes_path}: {line!r}. "
                    f"Expected {n_columns} columns, as on the first line."
                )
            barcode, celltype = parts[0].strip(), parts[1].strip()
            row = [celltype or None]
            if n_columns == 3:
                try:
                    row.append(parse_empty_status(parts[2]))
                except ValueError as e:
                    raise ValueError(
                        f"Malformed line in celltypes file {celltypes_path}: {line!r}. {e}"
                    )
            if not celltype and not (n_columns == 3 and row[1]):
                raise ValueError(
                    f"Malformed line in celltypes file {celltypes_path}: {line!r}. "
                    "Only a barcode marked empty may have a blank celltype."
                )
            mapping[barcode] = row
    columns = ["celltype", "is_empty"] if n_columns == 3 else ["celltype"]
    return pd.DataFrame.from_dict(mapping, orient="index", columns=columns)


def assign_celltypes_from_file(
    adata: ad.AnnData, celltypes_path: str, celltype_key: str = "celltype"
) -> bool:
    """Set the celltypes and empty droplets of an AnnData from a celltypes file.

    Barcodes of `adata` missing from the file are treated as empty droplets, and
    barcodes in the file missing from `adata` are ignored. `adata.obs["is_empty"]`
    is left alone only if the file has no empty status column and lists every
    barcode of `adata`, as it then says nothing about which are empty.

    Args:
        adata: AnnData to modify in place
        celltypes_path: Path to the mapping file (see `read_celltypes`)
        celltype_key: Column of `adata.obs` to store the celltypes in

    Returns:
        Whether `adata.obs["is_empty"]` was set from the file

    Raises:
        ValueError: If no barcode is in both, or no barcode is a real cell
    """
    if celltype_key in adata.obs.columns:
        logger.warning(f"Existing column '{celltype_key}' of the input AnnData will be replaced by the celltypes from {celltypes_path}")
    table = read_celltypes(celltypes_path)

    listed = adata.obs_names.isin(table.index)
    if not listed.any():
        raise ValueError(f"No barcodes in the count matrix were found in celltypes file {celltypes_path}.")
    n_ignored = int((~table.index.isin(adata.obs_names)).sum())
    if n_ignored > 0:
        logger.warning(f"Ignoring {n_ignored} of {len(table)} barcodes in {celltypes_path} that are not in the count matrix.")
    n_missing = int((~listed).sum())
    if n_missing > 0:
        logger.warning(f"{n_missing} of {adata.n_obs} barcodes in the count matrix are not in {celltypes_path} and will be treated as empty droplets.")

    celltypes = table["celltype"].reindex(adata.obs_names)
    has_status = "is_empty" in table.columns
    if not has_status and n_missing == 0:
        adata.obs[celltype_key] = celltypes.values
        return False

    is_empty = ~listed
    if has_status:
        is_empty[listed] = table.loc[adata.obs_names[listed], "is_empty"].values.astype(bool)
    if is_empty.all():
        raise ValueError(f"Every barcode is empty according to celltypes file {celltypes_path}, so there are no cells to denoise.")
    if "is_empty" in adata.obs.columns:
        logger.warning(f"Existing column 'is_empty' of the input AnnData will be replaced based on {celltypes_path}")
    # cellsweep averages celltype profiles by label alone, so an empty droplet
    # must not carry a celltype.
    celltypes[is_empty] = None
    adata.obs[celltype_key] = celltypes.values
    adata.obs["is_empty"] = is_empty
    return True


@logger.namespaced("sweep")
def sweep(
    kb_count_dir: str,
    out: Optional[str] = None,
    h5ad: bool = False,
    celltypes_path: Optional[str] = None,
    celltype_column: str = "celltype",
    leiden_resolution: Optional[float] = None,
    round_X: bool = False,
    keep_empties: bool = False,
    threads: int = 1,
    freeze_ambient_profile: bool = True,
    empty_droplet_method: str = "kneed",
    umi_cutoff: Optional[int] = None,
    expected_cells: Optional[int] = None,
    random_state: Optional[int] = 42,
    verbose: int = 0,
    quiet: bool = False,
    log_file: Optional[str] = None,
    overwrite: bool = False,
    **em_kwargs,
):
    """
    Wraps cellsweep

    Denoise a count matrix using an Expectation-Maximization (EM) algorithm that
    models each observed count as a mixture of ambient RNA, bulk RNA, and true cell-type
    signal.

    cellsweep infers the empty droplets, computes a mean expression profile for each
    celltype, and then iteratively estimates per-cell ambient fractions (alpha_i), a
    bulk contamination factor (beta), per-cell-type expression profiles (p_k), and the
    ambient contamination profile (a) until convergence.

    Parameters
    ----------
    kb_count_dir : str
        Path to one of the following:
        1. Output directory of kb count
        2. An AnnData (.h5ad) file containing the unfiltered count matrix

    out : str, default <kb_count_dir>/counts_unfiltered/adata_denoised.h5ad
        Path to write the denoised AnnData object (must end with `.h5ad`).

    h5ad : bool, default False
        If True and `kb_count_dir` is a directory, read the count matrix from
        `counts_unfiltered/adata.h5ad` instead of the .mtx and associated files.

    celltypes_path : str | None, default None
        Path to a tab-separated text file with columns barcode, celltype and,
        optionally, empty status (true/false). Takes precedence over `celltype_column`
        and `leiden_resolution`. The celltypes are stored in `adata.obs["celltype"]`.
        Barcodes of the count matrix missing from the file are treated as empty
        droplets (with a warning), and barcodes in the file missing from the count
        matrix are ignored (with a warning). The empty status is stored in
        `adata.obs["is_empty"]`, unless the file has no third column and lists every
        barcode, in which case empty droplets are inferred as usual.

    celltype_column : str, default "celltype"
        Column of the input AnnData's `obs` holding celltypes (h5ad input only).
        Used if `celltypes_path` is not provided; cellsweep reads the celltypes
        directly from this column.

    leiden_resolution : float | None, default None
        Resolution parameter for Leiden clustering. Celltypes are assigned by
        Leiden clustering only if neither `celltypes_path` nor `celltype_column`
        (h5ad input) supplies them. The clusters are stored in `adata.obs["celltype"]`.

    round_X : bool, default False
        If True, rounds denoised counts to nearest integer (stochastic rounding
        seeded by `random_state`) before saving.

    keep_empties : bool, default False
        If False, the empty droplets are removed from the output after denoising, so
        it contains only the real cells. If True, all input barcodes are kept; empty
        droplets then have `alpha_hat` = 1 and `z_hat` = -1. The empty droplets are
        always used to fit the model regardless of this setting.

    threads : int, default 1
        number of numba threads

    freeze_ambient_profile: bool, default True
        If True, the ambient profile (a) is anchored on the empty droplets and
        re-estimated from them only. If False, the ambient profile is modeled as a
        mixture of cell-type profiles whose weights are updated during training.

    empty_droplet_method : str, default "kneed"
        Strategy to infer empty droplets if `adata.obs["is_empty"]` is not present.
        "kneed" and "mx_filter" (see https://github.com/cellatlas/mx, requires
        scikit-learn) detect the UMI cutoff automatically. "manual" requires
        `umi_cutoff` or `expected_cells`, and is selected automatically when
        either is given.

    umi_cutoff : int | None, default None
        Optional absolute UMI count threshold for classifying droplets as empty.

    expected_cells : int | None, default None
        Expected number of real cells, used when estimating thresholds.

    random_state: int | None, default 42
        Random seed for stochastic rounding. Only used if `round_X=True`.

    verbose : int, default 0
        Verbosity level (2 debug, 1 info, 0 warning, -1 error, -2 critical).

    quiet : bool, default False
        Suppresses most log output when True.

    log_file : str | None, default None
        Optional path to save EM iteration logs.

    overwrite : bool, default False
        If True, overwrite `out` if it already exists.

    **em_kwargs
        Advanced EM hyperparameters forwarded unchanged to
        `cellsweep.denoise_count_matrix` (see `EM_KWARGS` and the cellsweep
        documentation): `init_alpha`, `init_beta`, `max_iter`, `del0_ll_tol`,
        `min_ll_tol`, `tol_p`, `tol_f`, `burnin_patience`, `burnin_max_iter`,
        `alpha_cap`, `repulsion_strength`, `max_frac_gene_repulsion`,
        `beta_prior_mode`, `beta_prior_strength`, `celltype_lambda`,
        `ambient_lambda`, `bulk_lambda`, `eps`, `log_eps`, `celltype_profile_key`,
        `ambient_profile_key`, `bulk_profile_key`. Any not given take cellsweep's
        defaults.

    Returns
    -------
    AnnData
        Denoised AnnData object with updated `adata.X`, and added fields. Unless
        `keep_empties=True`, it holds only the real cells (empty droplets removed).
        - `adata.layers["raw"]` : raw count matrix
        - `adata.obs["is_empty"]` : whether each barcode was treated as an empty droplet
        - `adata.obs["alpha_hat"]` : final optimized alpha (ambient fraction) values
        - `adata.obs["contamination_fraction"]` : total contamination fraction per cell,
          (1 - beta) * alpha + beta
        - `adata.obs["z_hat"]` : final cell-type assignments
        - `adata.var["ambient_hat"]` : final optimized ambient distribution
        - `adata.var["bulk_hat"]` : global (bulk) noise distribution
        - `adata.uns["p_hat"]` : final optimized matrix of cell-type profiles (K x G)
        - `adata.uns["beta_hat"]` : final optimized beta
        - `adata.uns["loglike"]` : final log-likelihood (note that this value is not the
           complete log-likelihood, only the relative log-likelihood)

    Notes
    -----
    The EM algorithm proceeds by:
      1. E-step: Update expected value of true, ambient noise, and bulk noise counts for each cell and gene.
      2. M-step: Update parameters (alpha, beta, p_k, a).
      3. Iterate until convergence (relative change in ll < `del0_ll_tol`) or reaching `max_iter`.
    """
    check_cellsweep_version()
    import cellsweep

    unknown = sorted(set(em_kwargs) - set(EM_KWARGS))
    if unknown:
        raise TypeError(
            f"sweep() got unexpected keyword argument(s): {unknown}. "
            f"Supported EM keyword arguments are: {list(EM_KWARGS)}"
        )

    #* load adata
    logger.info("Loading count matrix into AnnData object...")
    if not isinstance(kb_count_dir, str):
        raise ValueError(f"kb_count_dir must be a string path to a directory or .h5ad file, but got {type(kb_count_dir)}")
    if os.path.isdir(kb_count_dir):
        counts_dir = os.path.join(kb_count_dir, "counts_unfiltered")
        if not os.path.exists(counts_dir):
            raise ValueError(f"Provided kb_count_dir path {kb_count_dir} does not contain 'counts_unfiltered' directory.")
        if h5ad:
            h5ad_path = os.path.join(counts_dir, "adata.h5ad")
            if not os.path.isfile(h5ad_path):
                raise ValueError(f"--h5ad was specified but {h5ad_path} does not exist.")
            adata = ad.read_h5ad(h5ad_path)
        else:
            matrix_path = os.path.join(counts_dir, "cells_x_genes.mtx")
            barcodes_path = os.path.join(counts_dir, "cells_x_genes.barcodes.txt")
            genes_path = os.path.join(counts_dir, "cells_x_genes.genes.names.txt")
            adata = import_matrix_as_anndata(matrix_path, barcodes_path, genes_path)

        if out is None:
            out = os.path.join(counts_dir, "adata_denoised.h5ad")
    elif os.path.isfile(kb_count_dir):
        if not kb_count_dir.endswith(".h5ad"):
            raise ValueError(f"Provided kb_count_dir file {kb_count_dir} is not an .h5ad file.")
        h5ad = True
        adata = ad.read_h5ad(kb_count_dir)

        if out is None:
            out = kb_count_dir[:-len(".h5ad")] + "_denoised.h5ad"
    else:
        raise ValueError(f"Provided kb_count_dir path {kb_count_dir} is neither a directory nor a file.")

    if not out.endswith(".h5ad"):
        raise ValueError(f"Output path {out} must end with '.h5ad'.")
    if os.path.exists(out):
        if overwrite:
            logger.warning(f"Output file {out} already exists and will be overwritten.")
        else:
            raise FileExistsError(f"Output file {out} already exists. Set overwrite=True to overwrite it.")

    #* assign celltypes: (1) celltypes file, (2) h5ad obs column, (3) leiden clustering
    # `celltype_key` is the obs column cellsweep reads the celltypes from.
    if celltypes_path is not None:
        logger.info(f"Assigning celltypes from {celltypes_path}")
        celltype_key = "celltype"
        empties_from_file = assign_celltypes_from_file(adata, celltypes_path, celltype_key)
        if empties_from_file and (umi_cutoff is not None or expected_cells is not None):
            logger.warning(f"Ignoring umi_cutoff and expected_cells, since {celltypes_path} determines the empty droplets.")
    elif h5ad and celltype_column in adata.obs.columns:
        logger.info(f"Using celltypes from column '{celltype_column}' of the input AnnData")
        celltype_key = celltype_column
    elif leiden_resolution is not None:
        logger.info("Preprocessing and clustering with Scanpy to assign cell types...")
        celltype_key = "celltype"
        # Pass a copy: the preprocessing function mutates adata in place
        # (normalize_total, log1p, gene filtering). We only want the leiden
        # labels, so the raw counts in `adata` must be preserved for denoising.
        adata_processed_tmp = cellsweep.utils.run_scanpy_preprocessing_and_clustering(
            adata=adata.copy(),
            min_genes=None,
            min_cells=3,
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
        adata.obs[celltype_key] = adata_processed_tmp.obs["leiden"].reindex(adata.obs.index)
        del adata_processed_tmp
    else:
        if h5ad:
            raise ValueError(
                f"No celltypes available: column '{celltype_column}' is not in the input AnnData. "
                "Provide a celltypes file with -c, a different --celltype-column, or --leiden-resolution to cluster."
            )
        raise ValueError(
            "No celltypes available. Provide a celltypes file with -c, "
            "an h5ad input with a celltype column (--h5ad), or --leiden-resolution to cluster."
        )

    #* run cellsweep
    logger.info("Running CellSweep denoising...")
    adata_cellsweep = cellsweep.denoise_count_matrix(
        adata=adata,
        adata_out=out,
        round_X=round_X,
        keep_empties=keep_empties,
        threads=threads,
        freeze_ambient_profile=freeze_ambient_profile,
        empty_droplet_method=empty_droplet_method,
        umi_cutoff=umi_cutoff,
        expected_cells=expected_cells,
        random_state=random_state,
        inplace=True,
        verbose=verbose,
        quiet=quiet,
        log_file=log_file,
        celltype_key=celltype_key,
        **em_kwargs,
    )

    logger.info("CellSweep denoising complete.")
    return adata_cellsweep
