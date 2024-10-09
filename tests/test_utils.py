import os
import subprocess as sp
from unittest import mock, TestCase
from unittest.mock import call
from unittest.mock import ANY

import anndata
import numpy as np
import pandas as pd

import kb_python.utils as utils
from kb_python.config import UnsupportedOSError
from tests.mixins import TestMixin


@utils.restore_cwd
def change_cwd_func(temp_dir):
    path = os.path.join(temp_dir, 'test')
    os.makedirs(path)
    os.chdir(path)


class TestUtils(TestMixin, TestCase):

    def test_update_filename(self):
        self.assertEqual(
            'output.s.c.bus', utils.update_filename('output.s.bus', 'c')
        )

    def test_make_directory(self):
        with mock.patch('kb_python.utils.os.makedirs') as makedirs:
            utils.make_directory('path')
            makedirs.assert_called_once_with('path', exist_ok=True)

    def test_remove_directory(self):
        with mock.patch('kb_python.utils.shutil.rmtree') as rmtree:
            utils.remove_directory('path')
            rmtree.assert_called_once_with('path', ignore_errors=True)

    def test_run_executable(self):
        with mock.patch('kb_python.utils.STATS') as STATS:
            p, stdout, stderr = utils.run_executable(['echo', 'TEST'],
                                                     stream=False)
            self.assertEqual(stdout, 'TEST\n')
            STATS.command.assert_called_once_with(['echo', 'TEST'], runtime=ANY)

    def test_run_exectuable_raises_exception(self):
        with self.assertRaises(sp.SubprocessError):
            utils.run_executable(['bash', 'nonexistent option'])

    def test_run_exectuable_with_returncode(self):
        utils.run_executable(['bash', 'nonexistent option'], returncode=127)

    def test_run_executable_no_wait(self):
        with mock.patch('kb_python.utils.sp') as sp_mock:
            sp_mock.Popen().returncode = 0
            utils.run_executable(['echo', 'TEST'], wait=False)
            sp_mock.Popen().poll.assert_not_called()

    def test_run_executable_with_stream(self):
        with mock.patch('kb_python.utils.logger.debug') as debug_mock:
            utils.run_executable(['echo', 'TEST'], stream=True)
            debug_mock.assert_has_calls([call('TEST')])

    def test_get_kallisto_version(self):
        with mock.patch('kb_python.utils.run_executable') as run_executable:
            run_executable.return_value = None, 'kallisto 1.2.3', None
            self.assertEqual((1, 2, 3), utils.get_kallisto_version())

    def test_get_bustools_version(self):
        with mock.patch('kb_python.utils.run_executable') as run_executable:
            run_executable.return_value = None, 'bustools 1.2.3', None
            self.assertEqual((1, 2, 3), utils.get_bustools_version())

    def test_parse_technologies(self):
        lines = [
            'short name       description',
            '----------       -----------',
            '10xv1            10x version 1 chemistry',
            '10xv2            10x version 2 chemistry',
        ]
        self.assertEqual(utils.parse_technologies(lines), {'10xv1', '10xv2'})

    def test_get_supported_technologies(self):
        with mock.patch('kb_python.utils.run_executable') as run_executable,\
           mock.patch('kb_python.utils.parse_technologies') as parse_technologies:
            run_executable.return_value = None, 'TEST', None
            utils.get_supported_technologies()
            parse_technologies.assert_called_once_with('TEST')

    def test_whitelist_provided(self):
        self.assertTrue(utils.whitelist_provided('10xv2'))
        self.assertFalse(utils.whitelist_provided('UNSUPPORTED'))

    def test_stream_file(self):
        with mock.patch('kb_python.utils.PLATFORM', 'linux'),\
            mock.patch('kb_python.utils.os') as os,\
            mock.patch('kb_python.utils.threading') as threading,\
            mock.patch('kb_python.utils.urlretrieve') as urlretrieve:
            url = mock.MagicMock()
            path = mock.MagicMock()
            utils.stream_file(url, path)
            os.mkfifo.assert_called_once_with(path)
            threading.Thread.assert_called_once_with(
                target=urlretrieve, args=(url, path), daemon=True
            )
            threading.Thread().start.assert_called_once_with()

    def test_stream_file_windows(self):
        with mock.patch('kb_python.utils.PLATFORM', 'windows'),\
            mock.patch('kb_python.utils.os') as os,\
            mock.patch('kb_python.utils.threading') as threading,\
            mock.patch('kb_python.utils.urlretrieve') as urlretrieve:
            url = mock.MagicMock()
            path = mock.MagicMock()
            with self.assertRaises(UnsupportedOSError):
                utils.stream_file(url, path)
            os.mkfifo.assert_not_called()
            threading.thread.assert_not_called()
            urlretrieve.assert_not_called()

    def test_read_t2g(self):
        t2g = utils.read_t2g(self.t2g_path)
        self.assertEqual(16457, len(t2g))
        self.assertIn('ENSMUST00000158488.1', t2g)
        self.assertEqual(('ENSMUSG00000089113.1', 'Gm22902'),
                         t2g['ENSMUST00000158488.1'])

    def test_import_tcc_matrix_as_anndata(self):
        adata = utils.import_tcc_matrix_as_anndata(
            self.tcc_matrix_path, self.tcc_barcodes_path, self.tcc_ec_path,
            self.tcc_txnames_path
        )
        self.assertIsInstance(adata, anndata.AnnData)
        self.assertEqual({'transcript_ids'}, set(adata.var))
        self.assertIn(';', adata.var.iloc[-1]['transcript_ids'])
        self.assertEqual(set(), set(adata.obs))
        self.assertEqual('ec', adata.var.index.name)
        self.assertEqual('barcode', adata.obs.index.name)

    def test_import_matrix_as_anndata(self):
        adata = utils.import_matrix_as_anndata(
            self.matrix_path, self.barcodes_path, self.genes_path
        )
        self.assertIsInstance(adata, anndata.AnnData)
        self.assertEqual((29, 17), adata.shape)
        self.assertEqual(set(), set(adata.var))
        self.assertEqual(set(), set(adata.obs))
        self.assertEqual('gene_id', adata.var.index.name)
        self.assertEqual('barcode', adata.obs.index.name)
        self.assertEqual(1, adata.X[9, 0])
        self.assertEqual(1, adata.X[15, 0])

    def test_import_matrix_as_anndata_duplicated(self):
        adata = utils.import_matrix_as_anndata(
            self.matrix_duplicated_path,
            self.barcodes_path,
            self.genes_duplicated_path,
        )
        self.assertIsInstance(adata, anndata.AnnData)
        self.assertEqual(set(), set(adata.obs))
        self.assertEqual('gene_id', adata.var.index.name)
        self.assertEqual('barcode', adata.obs.index.name)

        self.assertEqual((29, 16), adata.shape)
        self.assertNotIn('ENSMUSG00000092572.7', adata.var.index)
        self.assertEqual(
            5, adata.X[15, adata.var.index.get_loc('ENSMUSG00000026034.17')]
        )

    # def test_import_matrix_as_anndata_with_t2g(self):
    #     adata = utils.import_matrix_as_anndata(
    #         self.matrix_path,
    #         self.barcodes_path,
    #         self.genes_path,
    #         t2g_path=self.t2g_path
    #     )
    #     self.assertIsInstance(adata, anndata.AnnData)
    #     self.assertEqual(set(), set(adata.obs))
    #     self.assertEqual('gene_id', adata.var.index.name)
    #     self.assertEqual('barcode', adata.obs.index.name)
    # 
    #     self.assertEqual([
    #         'Clk1', 'Serpinb10', 'Olfr421-ps1', 'Olfr335-ps', 'Olfr1001-ps1',
    #         'Olfr1010', 'Olfr1021-ps1', 'Olfr1038-ps', 'Olfr1077-ps1',
    #         'Olfr1083-ps', 'Olfr1117-ps1', 'Olfr1165-ps', 'Olfr475-ps1',
    #         'Olfr1267-ps1', 'Olfr1268-ps1', 'Olfr1273-ps', 'Olfr1300-ps1'
    #     ], list(adata.var.gene_id.values))
    # 
    # def test_import_matrix_as_anndata_with_t2g_no_gene_name(self):
    #     adata = utils.import_matrix_as_anndata(
    #         self.matrix_path,
    #         self.barcodes_path,
    #         self.genes_path,
    #         t2g_path=self.t2g_path2
    #     )
    #     self.assertIsInstance(adata, anndata.AnnData)
    #     self.assertEqual(set(), set(adata.obs))
    #     self.assertEqual('gene_id', adata.var.index.name)
    #     self.assertEqual('barcode', adata.obs.index.name)
    # 
    #     self.assertEqual([
    #         'Clk1', 'ENSMUSG00000092572.7', 'Olfr421-ps1', 'Olfr335-ps',
    #         'Olfr1001-ps1', 'Olfr1010', 'Olfr1021-ps1', 'Olfr1038-ps',
    #         'Olfr1077-ps1', 'Olfr1083-ps', 'Olfr1117-ps1', 'Olfr1165-ps',
    #         'Olfr475-ps1', 'Olfr1267-ps1', 'Olfr1268-ps1', 'Olfr1273-ps',
    #         'Olfr1300-ps1'
    #     ], list(adata.var.gene_id.values))

    def test_import_matrix_as_anndata_name(self):
        adata = utils.import_matrix_as_anndata(
            self.matrix_path, self.barcodes_path, self.genes_path, name='test'
        )
        self.assertIsInstance(adata, anndata.AnnData)
        self.assertEqual(set(), set(adata.var))
        self.assertEqual(set(), set(adata.obs))
        self.assertEqual('test_id', adata.var.index.name)
        self.assertEqual('barcode', adata.obs.index.name)

    def test_overlay_anndatas(self):
        adata_spliced = utils.import_matrix_as_anndata(
            self.spliced_matrix_path, self.spliced_barcodes_path,
            self.spliced_genes_path
        )
        adata_unspliced = utils.import_matrix_as_anndata(
            self.unspliced_matrix_path, self.unspliced_barcodes_path,
            self.unspliced_genes_path
        )
        adata = utils.overlay_anndatas(adata_spliced, adata_unspliced)
        self.assertEqual({'mature', 'nascent'}, set(adata.layers.keys()))

    def test_sum_anndatas(self):
        adata_spliced = utils.import_matrix_as_anndata(
            self.spliced_matrix_path, self.spliced_barcodes_path,
            self.spliced_genes_path
        )
        adata_unspliced = utils.import_matrix_as_anndata(
            self.unspliced_matrix_path, self.unspliced_barcodes_path,
            self.unspliced_genes_path
        )
        adata = utils.sum_anndatas(adata_spliced, adata_unspliced)
        self.assertEqual(2.0, adata.X[5, 15])

    def test_move_file(self):
        with mock.patch('kb_python.utils.shutil.move') as move:
            utils.move_file('source', 'destination')
            move.assert_called_once_with('source', 'destination')

    def test_copy_whitelist(self):
        whitelist_path = utils.copy_whitelist('10xv1', self.temp_dir)
        self.assertTrue(os.path.exists(whitelist_path))

    def test_restore_cwd(self):
        cwd = os.path.abspath(os.getcwd())
        change_cwd_func(self.temp_dir)
        self.assertEqual(cwd, os.path.abspath(os.getcwd()))

    def test_collapse_anndata_by_index(self):
        adata = anndata.AnnData(
            X=np.array([
                [0, 1, 2],
                [3, 4, 5],
            ]),
            layers={'layer': np.array([
                [6, 7, 8],
                [9, 10, 11],
            ])},
            obs=pd.DataFrame(index=pd.Series(['a', 'b'], name='cell')),
            var=pd.DataFrame(index=pd.Series(['c', 'c', 'd'], name='gene_id'))
        )

        collapsed = utils.collapse_anndata(adata)
        self.assertEqual((2, 2), collapsed.shape)
        pd.testing.assert_frame_equal(adata.obs, collapsed.obs)
        pd.testing.assert_index_equal(
            pd.Index(['c', 'd'], name='gene_id'), collapsed.var.index
        )
        np.testing.assert_array_equal(np.array([[1, 2], [7, 5]]), collapsed.X.A)
        np.testing.assert_array_equal(
            np.array([[13, 8], [19, 11]]), collapsed.layers['layer'].A
        )

    def test_collapse_anndata_by_column(self):
        adata = anndata.AnnData(
            X=np.array([
                [0, 1, 2],
                [3, 4, 5],
            ]),
            layers={'layer': np.array([
                [6, 7, 8],
                [9, 10, 11],
            ])},
            obs=pd.DataFrame(index=pd.Series(['a', 'b'], name='cell')),
            var=pd.DataFrame(
                index=pd.Series(['c', 'c', 'd'], name='gene_id'),
                data={'gene_name': ['e', 'f', 'f']}
            )
        )

        collapsed = utils.collapse_anndata(adata, by='gene_name')
        self.assertEqual((2, 2), collapsed.shape)
        pd.testing.assert_frame_equal(adata.obs, collapsed.obs)
        pd.testing.assert_index_equal(
            pd.Index(['e', 'f'], name='gene_name'), collapsed.var.index
        )
        np.testing.assert_array_equal(np.array([[0, 3], [3, 9]]), collapsed.X.A)
        np.testing.assert_array_equal(
            np.array([[6, 15], [9, 21]]), collapsed.layers['layer'].A
        )

#    def test_collapse_anndata_with_missing(self):
#        adata = anndata.AnnData(
#            X=np.array([
#                [0, 1, 2],
#                [3, 4, 5],
#            ]),
#            layers={'layer': np.array([
#                [6, 7, 8],
#                [9, 10, 11],
#            ])},
#            obs=pd.DataFrame(index=pd.Series(['a', 'b'], name='cell')),
#            var=pd.DataFrame(index=pd.Series(['c', None, 'c'], name='gene_id'))
#        )
#
#        collapsed = utils.collapse_anndata(adata)
#        self.assertEqual((2, 1), collapsed.shape)
#        pd.testing.assert_frame_equal(adata.obs, collapsed.obs)
#        pd.testing.assert_index_equal(
#            pd.Index(['c'], name='gene_id'), collapsed.var.index
#        )
#        np.testing.assert_array_equal(np.array([[2], [8]]), collapsed.X.A)
#        np.testing.assert_array_equal(
#            np.array([[14], [20]]), collapsed.layers['layer'].A
#        )

    def test_create_10x_feature_barcode_map(self):
        map_path = utils.create_10x_feature_barcode_map(
            os.path.join(self.temp_dir, 'map.txt')
        )
        self.assertTrue(os.path.exists(map_path))
        with open(map_path, 'r') as f:
            self.assertIn('\t', f.readline())
