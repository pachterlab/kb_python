import os
from unittest import skipUnless, TestCase

import anndata as ad
import numpy as np
import scipy.io
import scipy.sparse as sp

import kb_python.sweep as sweep
from tests.mixins import TestMixin

try:
    import cellsweep  # noqa: F401
    CELLSWEEP_INSTALLED = True
except ImportError:
    CELLSWEEP_INSTALLED = False


class TestSweep(TestMixin, TestCase):

    def make_matrix(self, n_real, n_empty, n_genes, seed=0):
        """Build a small synthetic unfiltered count matrix.

        Real cells get high counts and empty droplets near-zero counts, so
        threshold-based empty droplet detection is well-defined.
        """
        rng = np.random.default_rng(seed)
        real = rng.poisson(3.0, size=(n_real, n_genes))
        empty = rng.poisson(0.1, size=(n_empty, n_genes))
        X = sp.csr_matrix(np.vstack([real, empty]).astype(np.float32))
        barcodes = [f'bc{i}' for i in range(n_real + n_empty)]
        genes = [f'g{i}' for i in range(n_genes)]
        return X, barcodes, genes

    @skipUnless(CELLSWEEP_INSTALLED, 'cellsweep is not installed')
    def test_sweep_h5ad(self):
        # Pre-assign a `celltype` column so `sweep` skips the heavy Scanpy
        # clustering step and exercises the cellsweep denoising path directly.
        X, barcodes, genes = self.make_matrix(60, 40, 30)
        adata = ad.AnnData(X=X)
        adata.obs_names = barcodes
        adata.var_names = genes
        adata.obs['celltype'] = ['A', 'B'] * (adata.n_obs // 2)

        in_path = os.path.join(self.temp_dir, 'in.h5ad')
        out_path = os.path.join(self.temp_dir, 'out.h5ad')
        adata.write_h5ad(in_path)

        result = sweep.sweep(
            in_path,
            out=out_path,
            max_iter=5,
            threads=1,
            expected_cells=60,
            quiet=True,
        )

        # output written and matches the input dimensions
        self.assertTrue(os.path.exists(out_path))
        self.assertEqual((100, 30), result.shape)

        # documented fields are populated
        written = ad.read_h5ad(out_path)
        self.assertEqual((100, 30), written.shape)
        self.assertIn('raw', written.layers)
        for col in ('is_empty', 'cell_ambient_fraction', 'alpha_hat', 'z_hat'):
            self.assertIn(col, written.obs.columns)
        for key in ('p_hat', 'beta_hat', 'loglike'):
            self.assertIn(key, written.uns)

    @skipUnless(CELLSWEEP_INSTALLED, 'cellsweep is not installed')
    def test_sweep_kb_count_dir(self):
        # No `celltype` column here, so this additionally exercises matrix
        # loading from a kb count directory and the Scanpy clustering path.
        # Use a larger matrix so clustering (HVG / PCA / Leiden) is stable.
        X, barcodes, genes = self.make_matrix(400, 300, 300)
        counts_dir = os.path.join(self.temp_dir, 'counts_unfiltered')
        os.makedirs(counts_dir)
        scipy.io.mmwrite(os.path.join(counts_dir, 'cells_x_genes.mtx'), X)
        with open(
            os.path.join(counts_dir, 'cells_x_genes.barcodes.txt'), 'w'
        ) as f:
            f.write('\n'.join(barcodes) + '\n')
        with open(
            os.path.join(counts_dir, 'cells_x_genes.genes.names.txt'), 'w'
        ) as f:
            f.write('\n'.join(genes) + '\n')

        out_path = os.path.join(self.temp_dir, 'out.h5ad')
        result = sweep.sweep(
            self.temp_dir,
            out=out_path,
            max_iter=3,
            threads=1,
            expected_cells=400,
            quiet=True,
        )

        self.assertTrue(os.path.exists(out_path))
        self.assertEqual((700, 300), result.shape)
        self.assertIn('raw', result.layers)
