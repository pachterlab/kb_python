import os
import sys
from unittest import mock, skipUnless, TestCase

import anndata as ad
import numpy as np
import scipy.io
import scipy.sparse as sp

import kb_python.sweep as sweep
from tests.mixins import TestMixin

try:
    sweep.check_cellsweep_version()
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

    def write_counts_dir(self, X, barcodes, genes, adata=None):
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
        if adata is not None:
            adata.write_h5ad(os.path.join(counts_dir, 'adata.h5ad'))
        return counts_dir

    def run_mocked_sweep(self, *args, **kwargs):
        """Run `sweep` with cellsweep mocked out and return the AnnData that
        would have been passed to `cellsweep.denoise_count_matrix`."""
        fake_cellsweep = mock.MagicMock()
        fake_cellsweep.denoise_count_matrix.side_effect = (
            lambda adata, **kw: adata
        )
        clustered = mock.MagicMock()
        fake_cellsweep.utils.run_scanpy_preprocessing_and_clustering \
            .side_effect = lambda adata, **kw: setattr(
                clustered, 'obs', adata.obs.assign(leiden='0')
            ) or clustered
        with mock.patch.object(sweep, 'check_cellsweep_version'), \
                mock.patch.dict(sys.modules, {'cellsweep': fake_cellsweep}):
            result = sweep.sweep(*args, quiet=True, **kwargs)
        return result, fake_cellsweep

    def make_h5ad_dir(self, with_column=True):
        X, barcodes, genes = self.make_matrix(4, 2, 3)
        adata = ad.AnnData(X=X)
        adata.obs_names = barcodes
        adata.var_names = genes
        if with_column:
            adata.obs['ct'] = ['h5'] * adata.n_obs
        self.write_counts_dir(X, barcodes, genes, adata=adata)

    def write_celltypes(self):
        path = os.path.join(self.temp_dir, 'celltypes.txt')
        with open(path, 'w') as f:
            f.write('bc0\tT cell\nbc1\tB cell\n')
        return path

    def test_sweep_celltypes_file_takes_precedence(self):
        self.make_h5ad_dir()
        result, fake = self.run_mocked_sweep(
            self.temp_dir,
            h5ad=True,
            celltypes_path=self.write_celltypes(),
            celltype_column='ct',
            leiden_resolution=1.0,
        )
        self.assertEqual('T cell', result.obs.loc['bc0', 'celltype'])
        self.assertEqual('B cell', result.obs.loc['bc1', 'celltype'])
        self.assertTrue(result.obs['celltype'].iloc[2:].isna().all())
        _, kwargs = fake.denoise_count_matrix.call_args
        self.assertEqual('celltype', kwargs['celltype_key'])

    def test_sweep_celltype_column(self):
        self.make_h5ad_dir()
        result, fake = self.run_mocked_sweep(
            self.temp_dir, h5ad=True, celltype_column='ct',
            leiden_resolution=1.0,
        )
        # the column is passed to cellsweep as-is, not copied to 'celltype'
        self.assertTrue((result.obs['ct'] == 'h5').all())
        self.assertNotIn('celltype', result.obs.columns)
        _, kwargs = fake.denoise_count_matrix.call_args
        self.assertEqual('ct', kwargs['celltype_key'])
        fake.utils.run_scanpy_preprocessing_and_clustering.assert_not_called()

    def test_sweep_leiden(self):
        self.make_h5ad_dir(with_column=False)
        result, fake = self.run_mocked_sweep(
            self.temp_dir, h5ad=True, celltype_column='ct',
            leiden_resolution=0.5,
        )
        self.assertTrue((result.obs['celltype'] == '0').all())
        _, kwargs = fake.utils.run_scanpy_preprocessing_and_clustering \
            .call_args
        self.assertEqual(0.5, kwargs['leiden_resolution'])
        _, kwargs = fake.denoise_count_matrix.call_args
        self.assertEqual('celltype', kwargs['celltype_key'])

    def test_sweep_mtx_ignores_h5ad_column(self):
        # without --h5ad, the .mtx files are read, so there is no obs column
        self.make_h5ad_dir()
        with self.assertRaises(ValueError):
            self.run_mocked_sweep(self.temp_dir, celltype_column='ct')

    def test_sweep_no_celltypes_h5ad(self):
        self.make_h5ad_dir(with_column=False)
        with self.assertRaises(ValueError):
            self.run_mocked_sweep(self.temp_dir, h5ad=True)

    def test_sweep_h5ad_missing(self):
        X, barcodes, genes = self.make_matrix(4, 2, 3)
        self.write_counts_dir(X, barcodes, genes)
        with self.assertRaises(ValueError):
            self.run_mocked_sweep(self.temp_dir, h5ad=True, leiden_resolution=1.0)

    def test_sweep_em_kwargs_forwarded(self):
        self.make_h5ad_dir()
        _, fake = self.run_mocked_sweep(
            self.temp_dir, h5ad=True, celltype_column='ct',
            max_iter=7, init_beta=0.02, round_X=True, keep_empties=True,
        )
        _, kwargs = fake.denoise_count_matrix.call_args
        self.assertEqual(7, kwargs['max_iter'])
        self.assertEqual(0.02, kwargs['init_beta'])
        self.assertTrue(kwargs['round_X'])
        self.assertTrue(kwargs['keep_empties'])
        # EM kwargs that were not given are left to cellsweep's defaults
        self.assertNotIn('init_alpha', kwargs)
        self.assertNotIn('del0_ll_tol', kwargs)

    def test_sweep_unknown_kwarg(self):
        self.make_h5ad_dir()
        # arguments from cellsweep 0.1.0 that no longer exist
        for bad in ({'tol': 1e-3}, {'dirichlet_lambda': 500},
                    {'fixed_celltype': True}, {'integer_out': True}):
            with self.assertRaises(TypeError):
                self.run_mocked_sweep(
                    self.temp_dir, h5ad=True, celltype_column='ct', **bad
                )

    def test_sweep_out_not_h5ad(self):
        self.make_h5ad_dir()
        with self.assertRaises(ValueError):
            self.run_mocked_sweep(
                self.temp_dir, h5ad=True, celltype_column='ct',
                out=os.path.join(self.temp_dir, 'out.txt'),
            )

    def test_check_cellsweep_version(self):
        with mock.patch.object(sweep, 'version', return_value='0.1.0'):
            with self.assertRaises(ImportError):
                sweep.check_cellsweep_version()
        with mock.patch.object(sweep, 'version', return_value='1.0.0'):
            sweep.check_cellsweep_version()

    @skipUnless(CELLSWEEP_INSTALLED, 'cellsweep 1.0.0 is not installed')
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

        # output written; by default only the real cells are kept
        self.assertTrue(os.path.exists(out_path))
        self.assertEqual(30, result.n_vars)
        self.assertLess(result.n_obs, 100)
        self.assertFalse(result.obs['is_empty'].any())

        # documented fields are populated
        written = ad.read_h5ad(out_path)
        self.assertEqual(result.shape, written.shape)
        self.assertIn('raw', written.layers)
        for col in ('is_empty', 'contamination_fraction', 'alpha_hat',
                    'z_hat'):
            self.assertIn(col, written.obs.columns)
        for col in ('ambient_hat', 'bulk_hat'):
            self.assertIn(col, written.var.columns)
        for key in ('p_hat', 'beta_hat', 'loglike'):
            self.assertIn(key, written.uns)

    @skipUnless(CELLSWEEP_INSTALLED, 'cellsweep 1.0.0 is not installed')
    def test_sweep_h5ad_keep_empties_custom_column(self):
        # celltypes in a non-default column, read by cellsweep directly; with
        # keep_empties the output has every input barcode.
        X, barcodes, genes = self.make_matrix(60, 40, 30)
        adata = ad.AnnData(X=X)
        adata.obs_names = barcodes
        adata.var_names = genes
        adata.obs['annot'] = ['A', 'B'] * (adata.n_obs // 2)

        in_path = os.path.join(self.temp_dir, 'in.h5ad')
        out_path = os.path.join(self.temp_dir, 'out.h5ad')
        adata.write_h5ad(in_path)

        result = sweep.sweep(
            in_path,
            out=out_path,
            celltype_column='annot',
            keep_empties=True,
            round_X=True,
            max_iter=5,
            threads=1,
            expected_cells=60,
            quiet=True,
        )
        self.assertEqual((100, 30), result.shape)
        self.assertNotIn('celltype', result.obs.columns)
        self.assertEqual(40, int(result.obs['is_empty'].sum()))
        self.assertTrue((result.obs.loc[result.obs['is_empty'], 'z_hat'] == -1).all())
        # rounded output is integer-valued
        self.assertTrue(np.allclose(result.X.data, np.round(result.X.data)))

    @skipUnless(CELLSWEEP_INSTALLED, 'cellsweep 1.0.0 is not installed')
    def test_sweep_kb_count_dir(self):
        # No `celltype` column here, so this additionally exercises matrix
        # loading from a kb count directory and the Scanpy clustering path.
        # Use a larger matrix so clustering (HVG / PCA / Leiden) is stable.
        X, barcodes, genes = self.make_matrix(400, 300, 300)
        self.write_counts_dir(X, barcodes, genes)

        out_path = os.path.join(self.temp_dir, 'out.h5ad')
        result = sweep.sweep(
            self.temp_dir,
            out=out_path,
            max_iter=3,
            threads=1,
            expected_cells=400,
            leiden_resolution=1.0,
            quiet=True,
        )

        self.assertTrue(os.path.exists(out_path))
        self.assertEqual(300, result.n_vars)
        self.assertLess(result.n_obs, 700)
        self.assertIn('raw', result.layers)
