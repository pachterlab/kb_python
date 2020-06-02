import os
import tempfile
from unittest import mock, TestCase
from unittest.mock import ANY, call

import kb_python.count as count
from kb_python.constants import (
    ADATA_PREFIX,
    BUS_CDNA_PREFIX,
    BUS_FILENAME,
    BUS_FILTERED_FILENAME,
    BUS_FILTERED_SUFFIX,
    BUS_INTRON_PREFIX,
    BUS_S_FILENAME,
    BUS_SC_FILENAME,
    BUS_UNFILTERED_FILENAME,
    BUS_UNFILTERED_SUFFIX,
    CELLRANGER_DIR,
    CELLRANGER_BARCODES,
    CELLRANGER_GENES,
    CELLRANGER_MATRIX,
    COUNTS_PREFIX,
    ECMAP_FILENAME,
    FEATURE_PREFIX,
    FILTER_WHITELIST_FILENAME,
    FILTERED_COUNTS_DIR,
    INSPECT_FILENAME,
    KALLISTO_INFO_FILENAME,
    KB_INFO_FILENAME,
    REPORT_HTML_FILENAME,
    REPORT_NOTEBOOK_FILENAME,
    TXNAMES_FILENAME,
    UNFILTERED_COUNTS_DIR,
    WHITELIST_FILENAME,
)
from tests.mixins import TestMixin


class TestCount(TestMixin, TestCase):

    def setUp(self):
        makedirs_mock = mock.patch('kb_python.count.make_directory')
        makedirs_mock.start()
        self.addCleanup(makedirs_mock.stop)

    def test_kallisto_bus(self):
        out_dir = tempfile.mkdtemp()
        result = count.kallisto_bus(
            self.fastqs, self.index_path, self.technology, out_dir, threads=1
        )
        self.assertEqual({
            'bus': os.path.join(out_dir, BUS_FILENAME),
            'ecmap': os.path.join(out_dir, ECMAP_FILENAME),
            'txnames': os.path.join(out_dir, TXNAMES_FILENAME),
            'info': os.path.join(out_dir, KALLISTO_INFO_FILENAME)
        }, result)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_merge(self):
        out_dir = tempfile.mkdtemp()
        result = count.bustools_merge(self.bus_split_paths, out_dir, threads=1)
        self.assertEqual({
            'bus': os.path.join(out_dir, BUS_FILENAME),
            'ecmap': os.path.join(out_dir, ECMAP_FILENAME),
            'txnames': os.path.join(out_dir, TXNAMES_FILENAME),
            'info': os.path.join(out_dir, KALLISTO_INFO_FILENAME)
        }, result)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_project(self):
        out_dir = tempfile.mkdtemp()
        out_path = os.path.join(out_dir, 'projected.bus')
        result = count.bustools_project(
            self.bus_path, out_path, self.kite_map_path, self.ecmap_path,
            self.txnames_path
        )
        self.assertEqual({'bus': out_path}, result)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_sort(self):
        out_dir = tempfile.mkdtemp()
        out_path = os.path.join(out_dir, BUS_S_FILENAME)
        result = count.bustools_sort(
            self.bus_path, out_path, threads=1, memory='1G'
        )
        self.assertEqual({'bus': out_path}, result)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_inspect(self):
        out_dir = tempfile.mkdtemp()
        out_path = os.path.join(out_dir, INSPECT_FILENAME)
        result = count.bustools_inspect(
            self.bus_s_path, out_path, self.whitelist_path, self.ecmap_path
        )
        self.assertEqual({'inspect': out_path}, result)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_correct(self):
        out_dir = tempfile.mkdtemp()
        out_path = os.path.join(out_dir, BUS_SC_FILENAME)
        result = count.bustools_correct(
            self.bus_s_path, out_path, self.whitelist_path
        )
        self.assertEqual({'bus': out_path}, result)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_count(self):
        out_dir = tempfile.mkdtemp()
        counts_path = os.path.join(out_dir, COUNTS_PREFIX)
        result = count.bustools_count(
            self.bus_scs_path, counts_path, self.t2g_path, self.ecmap_path,
            self.txnames_path
        )
        self.assertEqual({
            'mtx': '{}.mtx'.format(counts_path),
            'genes': '{}.genes.txt'.format(counts_path),
            'barcodes': '{}.barcodes.txt'.format(counts_path),
        }, result)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_capture(self):
        out_dir = tempfile.mkdtemp()
        out_path = os.path.join(out_dir, 'capture.bus')
        result = count.bustools_capture(
            self.lamanno_bus_scs_path, out_path, self.lamanno_cdna_t2c_path,
            self.lamanno_txnames_path, self.lamanno_txnames_path
        )
        self.assertEqual({'bus': out_path}, result)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_whitelist(self):
        out_dir = tempfile.mkdtemp()
        out_path = os.path.join(out_dir, 'whitelist.txt')
        result = count.bustools_whitelist(self.bus_s_path, out_path)
        self.assertEqual({'whitelist': out_path}, result)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_convert_matrix_loom(self):
        with mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata,\
            mock.patch('kb_python.count.import_tcc_matrix_as_anndata') as import_tcc_matrix_as_anndata:
            counts_dir = 'path/to/counts/dir'
            matrix_path = mock.MagicMock()
            barcodes_path = mock.MagicMock()
            genes_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            loom_path = os.path.join(counts_dir, '{}.loom'.format(ADATA_PREFIX))
            self.assertEqual({'loom': loom_path},
                             count.convert_matrix(
                                 counts_dir,
                                 matrix_path,
                                 barcodes_path,
                                 t2g_path=t2g_path,
                                 genes_path=genes_path,
                                 loom=True
                             ))
            import_matrix_as_anndata.assert_called_once_with(
                matrix_path,
                barcodes_path,
                genes_path,
                t2g_path=t2g_path,
                name='gene'
            )
            import_matrix_as_anndata.return_value.write_loom.assert_called_once_with(
                loom_path
            )
            import_tcc_matrix_as_anndata.assert_not_called()

    def test_convert_matrix_h5ad(self):
        with mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata,\
            mock.patch('kb_python.count.import_tcc_matrix_as_anndata') as import_tcc_matrix_as_anndata:
            counts_dir = 'path/to/counts/dir'
            matrix_path = mock.MagicMock()
            barcodes_path = mock.MagicMock()
            genes_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            h5ad_path = os.path.join(counts_dir, '{}.h5ad'.format(ADATA_PREFIX))
            self.assertEqual({'h5ad': h5ad_path},
                             count.convert_matrix(
                                 counts_dir,
                                 matrix_path,
                                 barcodes_path,
                                 t2g_path=t2g_path,
                                 genes_path=genes_path,
                                 h5ad=True
                             ))
            import_matrix_as_anndata.assert_called_once_with(
                matrix_path,
                barcodes_path,
                genes_path,
                t2g_path=t2g_path,
                name='gene'
            )
            import_matrix_as_anndata.return_value.write.assert_called_once_with(
                h5ad_path
            )
            import_tcc_matrix_as_anndata.assert_not_called()

    def test_convert_matrix_tcc(self):
        with mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata,\
            mock.patch('kb_python.count.import_tcc_matrix_as_anndata') as import_tcc_matrix_as_anndata:
            counts_dir = 'path/to/counts/dir'
            matrix_path = mock.MagicMock()
            barcodes_path = mock.MagicMock()
            ec_path = mock.MagicMock()
            txnames_path = mock.MagicMock()
            loom_path = os.path.join(counts_dir, '{}.loom'.format(ADATA_PREFIX))
            self.assertEqual({'loom': loom_path},
                             count.convert_matrix(
                                 counts_dir,
                                 matrix_path,
                                 barcodes_path,
                                 ec_path=ec_path,
                                 txnames_path=txnames_path,
                                 loom=True,
                                 tcc=True
                             ))
            import_tcc_matrix_as_anndata.assert_called_once_with(
                matrix_path, barcodes_path, ec_path, txnames_path, threads=8
            )
            import_tcc_matrix_as_anndata.return_value.write_loom.assert_called_once_with(
                loom_path
            )
            import_matrix_as_anndata.assert_not_called()

    def test_convert_matrices_loom(self):
        with mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata,\
            mock.patch('kb_python.count.import_tcc_matrix_as_anndata') as import_tcc_matrix_as_anndata,\
            mock.patch('kb_python.count.overlay_anndatas') as overlay_anndatas,\
            mock.patch('kb_python.count.sum_anndatas') as sum_anndatas:
            counts_dir = 'path/to/counts/dir'
            matrix_paths = [mock.MagicMock(), mock.MagicMock()]
            barcodes_paths = [mock.MagicMock(), mock.MagicMock()]
            genes_paths = [mock.MagicMock(), mock.MagicMock()]
            t2g_path = mock.MagicMock()
            ec_paths = [None, None]
            loom_path = os.path.join(counts_dir, '{}.loom'.format(ADATA_PREFIX))
            adatas = [mock.MagicMock(), mock.MagicMock()]
            import_matrix_as_anndata.side_effect = adatas
            self.assertEqual({'loom': loom_path},
                             count.convert_matrices(
                                 counts_dir,
                                 matrix_paths,
                                 barcodes_paths,
                                 t2g_path=t2g_path,
                                 genes_paths=genes_paths,
                                 ec_paths=ec_paths,
                                 loom=True
                             ))
            self.assertEqual(2, import_matrix_as_anndata.call_count)
            import_matrix_as_anndata.assert_has_calls([
                call(
                    matrix_path,
                    barcode_path,
                    genes_path,
                    t2g_path=t2g_path,
                    name='gene'
                ) for matrix_path, barcode_path, genes_path in
                zip(matrix_paths, barcodes_paths, genes_paths)
            ])
            import_tcc_matrix_as_anndata.assert_not_called()
            overlay_anndatas.assert_called_once_with(*adatas)
            sum_anndatas.assert_not_called()
            overlay_anndatas.return_value.write_loom.assert_called_once_with(
                loom_path
            )

    def test_convert_matrices_h5ad(self):
        with mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata,\
            mock.patch('kb_python.count.import_tcc_matrix_as_anndata') as import_tcc_matrix_as_anndata,\
            mock.patch('kb_python.count.overlay_anndatas') as overlay_anndatas,\
            mock.patch('kb_python.count.sum_anndatas') as sum_anndatas:
            counts_dir = 'path/to/counts/dir'
            matrix_paths = [mock.MagicMock(), mock.MagicMock()]
            barcodes_paths = [mock.MagicMock(), mock.MagicMock()]
            genes_paths = [mock.MagicMock(), mock.MagicMock()]
            t2g_path = mock.MagicMock()
            h5ad_path = os.path.join(counts_dir, '{}.h5ad'.format(ADATA_PREFIX))
            adatas = [mock.MagicMock(), mock.MagicMock()]
            import_matrix_as_anndata.side_effect = adatas
            self.assertEqual({'h5ad': h5ad_path},
                             count.convert_matrices(
                                 counts_dir,
                                 matrix_paths,
                                 barcodes_paths,
                                 t2g_path=t2g_path,
                                 genes_paths=genes_paths,
                                 h5ad=True
                             ))
            self.assertEqual(2, import_matrix_as_anndata.call_count)
            import_matrix_as_anndata.assert_has_calls([
                call(
                    matrix_path,
                    barcode_path,
                    genes_path,
                    t2g_path=t2g_path,
                    name='gene'
                ) for matrix_path, barcode_path, genes_path in
                zip(matrix_paths, barcodes_paths, genes_paths)
            ])
            import_tcc_matrix_as_anndata.assert_not_called()
            overlay_anndatas.assert_called_once_with(*adatas)
            sum_anndatas.assert_not_called()
            overlay_anndatas.return_value.write.assert_called_once_with(
                h5ad_path
            )

    def test_convert_matrices_tcc(self):
        with mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata,\
            mock.patch('kb_python.count.import_tcc_matrix_as_anndata') as import_tcc_matrix_as_anndata,\
            mock.patch('kb_python.count.overlay_anndatas') as overlay_anndatas,\
            mock.patch('kb_python.count.sum_anndatas') as sum_anndatas:
            counts_dir = 'path/to/counts/dir'
            matrix_paths = [mock.MagicMock(), mock.MagicMock()]
            barcodes_paths = [mock.MagicMock(), mock.MagicMock()]
            ec_paths = [mock.MagicMock(), mock.MagicMock()]
            txnames_path = mock.MagicMock()
            loom_path = os.path.join(counts_dir, '{}.loom'.format(ADATA_PREFIX))
            adatas = [mock.MagicMock(), mock.MagicMock()]
            import_tcc_matrix_as_anndata.side_effect = adatas
            self.assertEqual({'loom': loom_path},
                             count.convert_matrices(
                                 counts_dir,
                                 matrix_paths,
                                 barcodes_paths,
                                 ec_paths=ec_paths,
                                 txnames_path=txnames_path,
                                 loom=True,
                                 tcc=True
                             ))
            self.assertEqual(2, import_tcc_matrix_as_anndata.call_count)
            import_tcc_matrix_as_anndata.assert_has_calls([
                call(
                    matrix_path, barcode_path, ec_path, txnames_path, threads=8
                ) for matrix_path, barcode_path, ec_path in
                zip(matrix_paths, barcodes_paths, ec_paths)
            ])
            import_matrix_as_anndata.assert_not_called()
            overlay_anndatas.assert_called_once_with(*adatas)
            sum_anndatas.assert_not_called()
            overlay_anndatas.return_value.write_loom.assert_called_once_with(
                loom_path
            )

    def test_convert_matrices_nucleus(self):
        with mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata,\
            mock.patch('kb_python.count.import_tcc_matrix_as_anndata') as import_tcc_matrix_as_anndata,\
            mock.patch('kb_python.count.overlay_anndatas') as overlay_anndatas,\
            mock.patch('kb_python.count.sum_anndatas') as sum_anndatas:
            counts_dir = 'path/to/counts/dir'
            matrix_paths = [mock.MagicMock(), mock.MagicMock()]
            barcodes_paths = [mock.MagicMock(), mock.MagicMock()]
            genes_paths = [mock.MagicMock(), mock.MagicMock()]
            t2g_path = mock.MagicMock()
            loom_path = os.path.join(counts_dir, '{}.loom'.format(ADATA_PREFIX))
            adatas = [mock.MagicMock(), mock.MagicMock()]
            import_matrix_as_anndata.side_effect = adatas
            self.assertEqual({'loom': loom_path},
                             count.convert_matrices(
                                 counts_dir,
                                 matrix_paths,
                                 barcodes_paths,
                                 t2g_path=t2g_path,
                                 genes_paths=genes_paths,
                                 loom=True,
                                 nucleus=True
                             ))
            self.assertEqual(2, import_matrix_as_anndata.call_count)
            import_matrix_as_anndata.assert_has_calls([
                call(
                    matrix_path,
                    barcode_path,
                    genes_path,
                    t2g_path=t2g_path,
                    name='gene'
                ) for matrix_path, barcode_path, genes_path in
                zip(matrix_paths, barcodes_paths, genes_paths)
            ])
            import_tcc_matrix_as_anndata.assert_not_called()
            sum_anndatas.assert_called_once_with(*adatas)
            overlay_anndatas.assert_not_called()
            sum_anndatas.return_value.write_loom.assert_called_once_with(
                loom_path
            )

    def test_filter_with_bustools(self):
        with mock.patch('kb_python.count.bustools_whitelist') as bustools_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.os.makedirs'),\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
            mock.patch('kb_python.count.update_filename') as update_filename:
            bus_path = 'path/to/bus.bus'
            ecmap_path = mock.MagicMock()
            txnames_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            counts_dir = 'counts'
            threads = 99999
            memory = 'memory'
            counts_prefix = os.path.join(counts_dir, COUNTS_PREFIX)
            temp_dir = tempfile.mkdtemp()
            whitelist_path = mock.MagicMock()
            capture_path = mock.MagicMock()
            sort_path = 'busfile.c.bus'
            update_filename.return_value = 'bus.c.bus'
            bustools_whitelist.return_value = {'whitelist': whitelist_path}
            bustools_correct.return_value = {'bus': capture_path}
            bustools_sort.return_value = {'bus': sort_path}
            bustools_count.return_value = {
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
            }
            self.assertEqual({
                'whitelist': whitelist_path,
                'bus_scs': sort_path,
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
            },
                             count.filter_with_bustools(
                                 bus_path,
                                 ecmap_path,
                                 txnames_path,
                                 t2g_path,
                                 whitelist_path,
                                 sort_path,
                                 counts_prefix,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory
                             ))

            bustools_whitelist.assert_called_once_with(bus_path, whitelist_path)
            bustools_correct.assert_called_once_with(
                bus_path,
                os.path.join(temp_dir, 'bus.c.bus'),
                whitelist_path,
            )
            bustools_sort.assert_called_once_with(
                capture_path,
                sort_path,
                temp_dir=temp_dir,
                threads=threads,
                memory=memory
            )
            bustools_count.assert_called_once_with(
                sort_path,
                counts_prefix,
                t2g_path,
                ecmap_path,
                txnames_path,
                tcc=False,
                mm=False
            )
            convert_matrix.assert_not_called()

    def test_filter_with_bustools_convert(self):
        with mock.patch('kb_python.count.bustools_whitelist') as bustools_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.os.makedirs'),\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
            mock.patch('kb_python.count.update_filename') as update_filename:
            bus_path = 'path/to/bus.bus'
            ecmap_path = mock.MagicMock()
            txnames_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            counts_dir = 'counts'
            threads = 99999
            memory = 'memory'
            counts_prefix = os.path.join(counts_dir, COUNTS_PREFIX)
            temp_dir = tempfile.mkdtemp()
            whitelist_path = mock.MagicMock()
            capture_path = mock.MagicMock()
            sort_path = 'path/to/busfile.bus'
            loom_path = mock.MagicMock()
            update_filename.return_value = 'bus.c.bus'
            bustools_whitelist.return_value = {'whitelist': whitelist_path}
            bustools_correct.return_value = {'bus': capture_path}
            bustools_sort.return_value = {'bus': sort_path}
            bustools_count.return_value = {
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
            }
            convert_matrix.return_value = {'loom': loom_path}
            self.assertEqual({
                'whitelist': whitelist_path,
                'bus_scs': sort_path,
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
                'loom': loom_path
            },
                             count.filter_with_bustools(
                                 bus_path,
                                 ecmap_path,
                                 txnames_path,
                                 t2g_path,
                                 whitelist_path,
                                 sort_path,
                                 counts_prefix,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                                 loom=True
                             ))

            bustools_whitelist.assert_called_once_with(bus_path, whitelist_path)
            bustools_correct.assert_called_once_with(
                bus_path,
                os.path.join(temp_dir, 'bus.c.bus'),
                whitelist_path,
            )
            bustools_sort.assert_called_once_with(
                capture_path,
                sort_path,
                temp_dir=temp_dir,
                threads=threads,
                memory=memory
            )
            bustools_count.assert_called_once_with(
                sort_path,
                counts_prefix,
                t2g_path,
                ecmap_path,
                txnames_path,
                tcc=False,
                mm=False
            )
            convert_matrix.assert_called_once_with(
                counts_dir,
                '{}.mtx'.format(counts_prefix),
                '{}.barcodes.txt'.format(counts_prefix),
                genes_path='{}.genes.txt'.format(counts_prefix),
                t2g_path=t2g_path,
                ec_path=None,
                txnames_path=txnames_path,
                name='gene',
                loom=True,
                h5ad=False,
                tcc=False,
                threads=threads
            )

    def test_filter_with_bustools_dont_count(self):
        with mock.patch('kb_python.count.bustools_whitelist') as bustools_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.os.makedirs'),\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
            mock.patch('kb_python.count.update_filename') as update_filename:
            bus_path = 'path/to/bus.bus'
            ecmap_path = mock.MagicMock()
            txnames_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            counts_dir = 'counts'
            threads = 99999
            memory = 'memory'
            counts_prefix = os.path.join(counts_dir, COUNTS_PREFIX)
            temp_dir = tempfile.mkdtemp()
            whitelist_path = mock.MagicMock()
            capture_path = mock.MagicMock()
            sort_path = 'path/to/busfile.bus'
            update_filename.return_value = 'bus.c.bus'
            bustools_whitelist.return_value = {'whitelist': whitelist_path}
            bustools_correct.return_value = {'bus': capture_path}
            bustools_sort.return_value = {'bus': sort_path}
            bustools_count.return_value = {
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
            }
            self.assertEqual({
                'whitelist': whitelist_path,
                'bus_scs': sort_path,
            },
                             count.filter_with_bustools(
                                 bus_path,
                                 ecmap_path,
                                 txnames_path,
                                 t2g_path,
                                 whitelist_path,
                                 sort_path,
                                 counts_prefix,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                                 count=False,
                             ))

            bustools_whitelist.assert_called_once_with(bus_path, whitelist_path)
            bustools_correct.assert_called_once_with(
                bus_path,
                os.path.join(temp_dir, 'bus.c.bus'),
                whitelist_path,
            )
            bustools_sort.assert_called_once_with(
                capture_path,
                sort_path,
                temp_dir=temp_dir,
                threads=threads,
                memory=memory
            )
            bustools_count.assert_not_called()
            convert_matrix.assert_not_called()

    def test_filter_with_bustools_tcc(self):
        with mock.patch('kb_python.count.bustools_whitelist') as bustools_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.os.makedirs'),\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
            mock.patch('kb_python.count.update_filename') as update_filename:
            bus_path = 'path/to/bus.bus'
            ecmap_path = mock.MagicMock()
            txnames_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            counts_dir = 'counts'
            threads = 99999
            memory = 'memory'
            counts_prefix = os.path.join(counts_dir, COUNTS_PREFIX)
            temp_dir = tempfile.mkdtemp()
            whitelist_path = mock.MagicMock()
            capture_path = mock.MagicMock()
            sort_path = 'path/to/busfile.bus'
            update_filename.return_value = 'bus.c.bus'
            bustools_whitelist.return_value = {'whitelist': whitelist_path}
            bustools_correct.return_value = {'bus': capture_path}
            bustools_sort.return_value = {'bus': sort_path}
            bustools_count.return_value = {
                'mtx': '{}.mtx'.format(counts_prefix),
                'ec': '{}.ec.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
            }
            self.assertEqual({
                'whitelist': whitelist_path,
                'bus_scs': sort_path,
                'mtx': '{}.mtx'.format(counts_prefix),
                'ec': '{}.ec.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
            },
                             count.filter_with_bustools(
                                 bus_path,
                                 ecmap_path,
                                 txnames_path,
                                 t2g_path,
                                 whitelist_path,
                                 sort_path,
                                 counts_prefix,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                                 tcc=True,
                             ))

            bustools_whitelist.assert_called_once_with(bus_path, whitelist_path)
            bustools_correct.assert_called_once_with(
                bus_path,
                os.path.join(temp_dir, 'bus.c.bus'),
                whitelist_path,
            )
            bustools_sort.assert_called_once_with(
                capture_path,
                sort_path,
                temp_dir=temp_dir,
                threads=threads,
                memory=memory
            )
            bustools_count.assert_called_once_with(
                sort_path,
                counts_prefix,
                t2g_path,
                ecmap_path,
                txnames_path,
                tcc=True,
                mm=False
            )
            convert_matrix.assert_not_called()

    def test_filter_with_bustools_cellranger(self):
        with mock.patch('kb_python.count.bustools_whitelist') as bustools_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.os.makedirs'),\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
            mock.patch('kb_python.count.update_filename') as update_filename,\
            mock.patch('kb_python.count.matrix_to_cellranger') as matrix_to_cellranger:
            bus_path = 'path/to/bus.bus'
            ecmap_path = mock.MagicMock()
            txnames_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            counts_dir = 'counts'
            cellranger_dir = os.path.join(counts_dir, CELLRANGER_DIR)
            threads = 99999
            memory = 'memory'
            counts_prefix = os.path.join(counts_dir, COUNTS_PREFIX)
            temp_dir = tempfile.mkdtemp()
            whitelist_path = mock.MagicMock()
            capture_path = mock.MagicMock()
            sort_path = 'busfile.c.bus'
            update_filename.return_value = 'bus.c.bus'
            bustools_whitelist.return_value = {'whitelist': whitelist_path}
            bustools_correct.return_value = {'bus': capture_path}
            bustools_sort.return_value = {'bus': sort_path}
            bustools_count.return_value = {
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
            }
            matrix_to_cellranger.return_value = {
                'mtx': os.path.join(cellranger_dir, CELLRANGER_MATRIX),
                'genes': os.path.join(cellranger_dir, CELLRANGER_GENES),
                'barcodes': os.path.join(cellranger_dir, CELLRANGER_BARCODES),
            }
            self.assertEqual({
                'whitelist': whitelist_path,
                'bus_scs': sort_path,
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
                'cellranger': {
                    'mtx':
                        os.path.join(cellranger_dir, CELLRANGER_MATRIX),
                    'genes':
                        os.path.join(cellranger_dir, CELLRANGER_GENES),
                    'barcodes':
                        os.path.join(cellranger_dir, CELLRANGER_BARCODES),
                }
            },
                             count.filter_with_bustools(
                                 bus_path,
                                 ecmap_path,
                                 txnames_path,
                                 t2g_path,
                                 whitelist_path,
                                 sort_path,
                                 counts_prefix,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                                 cellranger=True
                             ))

            bustools_whitelist.assert_called_once_with(bus_path, whitelist_path)
            bustools_correct.assert_called_once_with(
                bus_path,
                os.path.join(temp_dir, 'bus.c.bus'),
                whitelist_path,
            )
            bustools_sort.assert_called_once_with(
                capture_path,
                sort_path,
                temp_dir=temp_dir,
                threads=threads,
                memory=memory
            )
            bustools_count.assert_called_once_with(
                sort_path,
                counts_prefix,
                t2g_path,
                ecmap_path,
                txnames_path,
                tcc=False,
                mm=False
            )
            convert_matrix.assert_not_called()
            matrix_to_cellranger.assert_called_once_with(
                '{}.mtx'.format(counts_prefix),
                '{}.barcodes.txt'.format(counts_prefix),
                '{}.genes.txt'.format(counts_prefix),
                t2g_path,
                cellranger_dir,
            )

    def test_stream_fastqs_local(self):
        with mock.patch('kb_python.count.stream_file') as stream_file:
            temp_dir = tempfile.mkdtemp()
            fastqs = ['path/to/file1.gz', 'path/to/file2.gz']
            stream_file.side_effect = ['FILE 1', 'FILE 2']
            self.assertEqual(
                fastqs, count.stream_fastqs(fastqs, temp_dir=temp_dir)
            )
            stream_file.assert_not_called()

    def test_stream_fastqs_remote(self):
        with mock.patch('kb_python.count.stream_file') as stream_file:
            temp_dir = tempfile.mkdtemp()
            fastqs = ['http://path/to/file1.gz', 'https://path/to/file2.gz']
            local_fastqs = [
                os.path.join(temp_dir, os.path.basename(fastq))
                for fastq in fastqs
            ]
            stream_file.side_effect = ['FILE 1', 'FILE 2']
            self.assertEqual(['FILE 1', 'FILE 2'],
                             count.stream_fastqs(fastqs, temp_dir=temp_dir))
            self.assertEqual(2, stream_file.call_count)
            stream_file.assert_has_calls([
                call(fastqs[0], local_fastqs[0]),
                call(fastqs[1], local_fastqs[1]),
            ])

    def test_copy_or_create_whitelist_provided(self):
        with mock.patch('kb_python.count.copy_whitelist') as copy_whitelist,\
            mock.patch('kb_python.count.bustools_whitelist') as bustools_whitelist:
            out_dir = tempfile.mkdtemp()
            count.copy_or_create_whitelist(
                self.technology, self.bus_s_path, out_dir
            )
            copy_whitelist.assert_called_once_with(self.technology, out_dir)
            bustools_whitelist.assert_not_called()

    def test_copy_or_create_whitelist_not_provided(self):
        with mock.patch('kb_python.count.copy_whitelist') as copy_whitelist,\
            mock.patch('kb_python.count.bustools_whitelist') as bustools_whitelist:
            out_dir = tempfile.mkdtemp()
            count.copy_or_create_whitelist(
                'UNSUPPORTED', self.bus_s_path, out_dir
            )
            copy_whitelist.assert_not_called()
            bustools_whitelist.assert_called_once_with(
                self.bus_s_path, os.path.join(out_dir, WHITELIST_FILENAME)
            )

    def test_matrix_to_cellranger(self):
        out_dir = tempfile.mkdtemp()
        result = count.matrix_to_cellranger(
            self.matrix_path, self.barcodes_path, self.genes_path,
            self.t2g_path, out_dir
        )

        with open(result['mtx'], 'r') as f, open(self.cr_matrix_path,
                                                 'r') as mtx:
            self.assertEqual(mtx.read(), f.read())
        with open(result['barcodes'], 'r') as f, open(self.cr_barcodes_path,
                                                      'r') as bc:
            self.assertEqual(bc.read(), f.read())
        with open(result['genes'], 'r') as f, open(self.cr_genes_path,
                                                   'r') as g:
            self.assertEqual(g.read(), f.read())

    def test_count_with_whitelist(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_merge') as bustools_merge,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
            mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report') as render_report,\
            mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata:
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_prefix = os.path.join(
                out_dir, UNFILTERED_COUNTS_DIR, COUNTS_PREFIX
            )
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path
            }
            bustools_sort.side_effect = [{
                'bus': bus_s_path
            }, {
                'bus': bus_scs_path
            }]
            bustools_inspect.return_value = {'inspect': inspect_path}
            bustools_correct.return_value = {'bus': bus_sc_path}
            bustools_count.return_value = {
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
            }
            STATS.save.return_value = 'stats'

            self.assertEqual({
                'stats': 'stats',
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'info': info_path,
                    'inspect': inspect_path,
                    'bus_scs': bus_scs_path,
                    'mtx': '{}.mtx'.format(counts_prefix),
                    'genes': '{}.genes.txt'.format(counts_prefix),
                    'barcodes': '{}.barcodes.txt'.format(counts_prefix),
                }
            },
                             count.count(
                                 self.index_path,
                                 self.t2g_path,
                                 self.technology,
                                 out_dir,
                                 self.fastqs,
                                 whitelist_path=self.whitelist_path,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                             ))

            stream_fastqs.assert_called_once_with(
                self.fastqs, temp_dir=temp_dir
            )
            kallisto_bus.assert_called_once_with(
                self.fastqs,
                self.index_path,
                self.technology,
                out_dir,
                threads=threads
            )
            self.assertEqual(bustools_sort.call_count, 2)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path, inspect_path, self.whitelist_path, ecmap_path
            )
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            bustools_count.assert_called_once_with(
                bus_scs_path,
                counts_prefix,
                self.t2g_path,
                ecmap_path,
                txnames_path,
                tcc=False,
                mm=False
            )
            convert_matrix.assert_not_called()
            filter_with_bustools.assert_not_called()
            bustools_merge.assert_not_called()

            STATS.start.assert_called_once()
            STATS.end.assert_called_once()
            STATS.save.assert_called_once_with(
                os.path.join(out_dir, KB_INFO_FILENAME)
            )
            import_matrix_as_anndata.assert_not_called()
            render_report.assert_not_called()

    def test_count_report(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_merge') as bustools_merge,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
            mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report') as render_report,\
            mock.patch('kb_python.count.import_matrix_as_anndata'):
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_prefix = os.path.join(
                out_dir, UNFILTERED_COUNTS_DIR, COUNTS_PREFIX
            )
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            report_html_path = os.path.join(out_dir, REPORT_HTML_FILENAME)
            report_nb_path = os.path.join(out_dir, REPORT_NOTEBOOK_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path
            }
            bustools_sort.side_effect = [{
                'bus': bus_s_path
            }, {
                'bus': bus_scs_path
            }]
            bustools_inspect.return_value = {'inspect': inspect_path}
            bustools_correct.return_value = {'bus': bus_sc_path}
            bustools_count.return_value = {
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
            }
            render_report.return_value = {
                'report_notebook': report_nb_path,
                'report_html': report_html_path
            }
            STATS.save.return_value = 'stats'

            self.assertEqual({
                'stats': 'stats',
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'info': info_path,
                    'inspect': inspect_path,
                    'report_notebook': report_nb_path,
                    'report_html': report_html_path,
                    'bus_scs': bus_scs_path,
                    'mtx': '{}.mtx'.format(counts_prefix),
                    'genes': '{}.genes.txt'.format(counts_prefix),
                    'barcodes': '{}.barcodes.txt'.format(counts_prefix),
                }
            },
                             count.count(
                                 self.index_path,
                                 self.t2g_path,
                                 self.technology,
                                 out_dir,
                                 self.fastqs,
                                 whitelist_path=self.whitelist_path,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                                 report=True
                             ))

            stream_fastqs.assert_called_once_with(
                self.fastqs, temp_dir=temp_dir
            )
            kallisto_bus.assert_called_once_with(
                self.fastqs,
                self.index_path,
                self.technology,
                out_dir,
                threads=threads
            )
            self.assertEqual(bustools_sort.call_count, 2)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path, inspect_path, self.whitelist_path, ecmap_path
            )
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            bustools_count.assert_called_once_with(
                bus_scs_path,
                counts_prefix,
                self.t2g_path,
                ecmap_path,
                txnames_path,
                tcc=False,
                mm=False
            )
            convert_matrix.assert_not_called()
            filter_with_bustools.assert_not_called()
            bustools_merge.assert_not_called()

            STATS.start.assert_called_once()
            STATS.end.assert_called_once()
            render_report.assert_called_once_with(
                STATS.save(),
                info_path,
                inspect_path,
                report_nb_path,
                report_html_path,
                '{}.mtx'.format(counts_prefix),
                '{}.barcodes.txt'.format(counts_prefix),
                '{}.genes.txt'.format(counts_prefix),
                self.t2g_path,
                temp_dir=temp_dir
            )

    def test_count_split(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_merge') as bustools_merge,\
            mock.patch('kb_python.count.move_file') as move_file,\
            mock.patch('kb_python.count.get_temporary_filename') as get_temporary_filename,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
            mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report'),\
            mock.patch('kb_python.count.import_matrix_as_anndata'):
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_prefix = os.path.join(
                out_dir, UNFILTERED_COUNTS_DIR, COUNTS_PREFIX
            )
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            stream_fastqs.return_value = self.fastqs
            bustools_merge.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path
            }
            bustools_sort.side_effect = [{
                'bus': 'temp1'
            }, {
                'bus': 'temp2'
            }, {
                'bus': bus_s_path
            }, {
                'bus': bus_scs_path
            }]
            bustools_inspect.return_value = {'inspect': inspect_path}
            bustools_correct.return_value = {'bus': bus_sc_path}
            bustools_count.return_value = {
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
            }
            STATS.save.return_value = 'stats'

            self.assertEqual({
                'stats': 'stats',
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'info': info_path,
                    'inspect': inspect_path,
                    'bus_scs': bus_scs_path,
                    'mtx': '{}.mtx'.format(counts_prefix),
                    'genes': '{}.genes.txt'.format(counts_prefix),
                    'barcodes': '{}.barcodes.txt'.format(counts_prefix),
                }
            },
                             count.count(['index.0', 'index.1'],
                                         self.t2g_path,
                                         self.technology,
                                         out_dir,
                                         self.fastqs,
                                         whitelist_path=self.whitelist_path,
                                         temp_dir=temp_dir,
                                         threads=threads,
                                         memory=memory))
            self.assertEqual(2, stream_fastqs.call_count)
            stream_fastqs.assert_has_calls([
                call(self.fastqs, temp_dir=temp_dir),
                call(self.fastqs, temp_dir=temp_dir)
            ])
            self.assertEqual(2, kallisto_bus.call_count)
            kallisto_bus.assert_has_calls([
                call(
                    self.fastqs,
                    'index.0',
                    self.technology,
                    os.path.join(temp_dir, 'bus_part0'),
                    threads=threads,
                    n=True
                ),
                call(
                    self.fastqs,
                    'index.1',
                    self.technology,
                    os.path.join(temp_dir, 'bus_part1'),
                    threads=threads,
                    n=True
                )
            ])
            self.assertEqual(2, move_file.call_count)
            move_file.assert_has_calls([
                call(
                    'temp1', os.path.join(temp_dir, 'bus_part0', 'output.bus')
                ),
                call(
                    'temp2', os.path.join(temp_dir, 'bus_part1', 'output.bus')
                )
            ])
            bustools_merge.assert_called_once_with([
                os.path.join(temp_dir, 'bus_part0'),
                os.path.join(temp_dir, 'bus_part1')
            ],
                                                   out_dir,
                                                   threads=threads)
            self.assertEqual(bustools_sort.call_count, 4)
            bustools_sort.assert_has_calls([
                call(
                    os.path.join(temp_dir, 'bus_part0', 'output.bus'),
                    get_temporary_filename.return_value,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                    flags=True
                ),
                call(
                    os.path.join(temp_dir, 'bus_part1', 'output.bus'),
                    get_temporary_filename.return_value,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                    flags=True
                ),
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path, inspect_path, self.whitelist_path, ecmap_path
            )
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            bustools_count.assert_called_once_with(
                bus_scs_path,
                counts_prefix,
                self.t2g_path,
                ecmap_path,
                txnames_path,
                tcc=False,
                mm=False
            )
            convert_matrix.assert_not_called()
            filter_with_bustools.assert_not_called()

    def test_count_convert(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
            mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report'),\
            mock.patch('kb_python.count.import_matrix_as_anndata'):
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_prefix = os.path.join(
                out_dir, UNFILTERED_COUNTS_DIR, COUNTS_PREFIX
            )
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            loom_path = mock.MagicMock()
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path
            }
            bustools_sort.side_effect = [{
                'bus': bus_s_path
            }, {
                'bus': bus_scs_path
            }]
            bustools_inspect.return_value = {'inspect': inspect_path}
            bustools_correct.return_value = {'bus': bus_sc_path}
            bustools_count.return_value = {
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
            }
            convert_matrix.return_value = {'loom': loom_path}
            STATS.save.return_value = 'stats'

            self.assertEqual({
                'stats': 'stats',
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'info': info_path,
                    'inspect': inspect_path,
                    'bus_scs': bus_scs_path,
                    'mtx': '{}.mtx'.format(counts_prefix),
                    'genes': '{}.genes.txt'.format(counts_prefix),
                    'barcodes': '{}.barcodes.txt'.format(counts_prefix),
                    'loom': loom_path,
                }
            },
                             count.count(
                                 self.index_path,
                                 self.t2g_path,
                                 self.technology,
                                 out_dir,
                                 self.fastqs,
                                 whitelist_path=self.whitelist_path,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                                 loom=True,
                             ))
            stream_fastqs.assert_called_once_with(
                self.fastqs, temp_dir=temp_dir
            )
            kallisto_bus.assert_called_once_with(
                self.fastqs,
                self.index_path,
                self.technology,
                out_dir,
                threads=threads
            )
            self.assertEqual(bustools_sort.call_count, 2)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path, inspect_path, self.whitelist_path, ecmap_path
            )
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            bustools_count.assert_called_once_with(
                bus_scs_path,
                counts_prefix,
                self.t2g_path,
                ecmap_path,
                txnames_path,
                tcc=False,
                mm=False
            )
            convert_matrix.assert_called_once_with(
                os.path.join(out_dir, UNFILTERED_COUNTS_DIR),
                '{}.mtx'.format(counts_prefix),
                '{}.barcodes.txt'.format(counts_prefix),
                genes_path='{}.genes.txt'.format(counts_prefix),
                t2g_path=self.t2g_path,
                ec_path=None,
                txnames_path=txnames_path,
                name='gene',
                loom=True,
                h5ad=False,
                tcc=False,
                threads=threads
            )
            filter_with_bustools.assert_not_called()

    def test_count_cellranger(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_merge') as bustools_merge,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
            mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report') as render_report,\
            mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata,\
            mock.patch('kb_python.count.matrix_to_cellranger') as matrix_to_cellranger:
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_prefix = os.path.join(
                out_dir, UNFILTERED_COUNTS_DIR, COUNTS_PREFIX
            )
            cellranger_dir = os.path.join(
                out_dir, UNFILTERED_COUNTS_DIR, CELLRANGER_DIR
            )
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path
            }
            bustools_sort.side_effect = [{
                'bus': bus_s_path
            }, {
                'bus': bus_scs_path
            }]
            bustools_inspect.return_value = {'inspect': inspect_path}
            bustools_correct.return_value = {'bus': bus_sc_path}
            bustools_count.return_value = {
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
            }
            matrix_to_cellranger.return_value = {
                'mtx': os.path.join(cellranger_dir, CELLRANGER_MATRIX),
                'genes': os.path.join(cellranger_dir, CELLRANGER_GENES),
                'barcodes': os.path.join(cellranger_dir, CELLRANGER_BARCODES),
            }
            STATS.save.return_value = 'stats'

            self.assertEqual({
                'stats': 'stats',
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'info': info_path,
                    'inspect': inspect_path,
                    'bus_scs': bus_scs_path,
                    'mtx': '{}.mtx'.format(counts_prefix),
                    'genes': '{}.genes.txt'.format(counts_prefix),
                    'barcodes': '{}.barcodes.txt'.format(counts_prefix),
                    'cellranger': {
                        'mtx':
                            os.path.join(cellranger_dir, CELLRANGER_MATRIX),
                        'genes':
                            os.path.join(cellranger_dir, CELLRANGER_GENES),
                        'barcodes':
                            os.path.join(cellranger_dir, CELLRANGER_BARCODES),
                    }
                }
            },
                             count.count(
                                 self.index_path,
                                 self.t2g_path,
                                 self.technology,
                                 out_dir,
                                 self.fastqs,
                                 whitelist_path=self.whitelist_path,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                                 cellranger=True
                             ))

            stream_fastqs.assert_called_once_with(
                self.fastqs, temp_dir=temp_dir
            )
            kallisto_bus.assert_called_once_with(
                self.fastqs,
                self.index_path,
                self.technology,
                out_dir,
                threads=threads
            )
            self.assertEqual(bustools_sort.call_count, 2)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path, inspect_path, self.whitelist_path, ecmap_path
            )
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            bustools_count.assert_called_once_with(
                bus_scs_path,
                counts_prefix,
                self.t2g_path,
                ecmap_path,
                txnames_path,
                tcc=False,
                mm=False
            )
            convert_matrix.assert_not_called()
            filter_with_bustools.assert_not_called()
            bustools_merge.assert_not_called()

            STATS.start.assert_called_once()
            STATS.end.assert_called_once()
            STATS.to_dict.assert_not_called()
            import_matrix_as_anndata.assert_not_called()
            render_report.assert_not_called()
            matrix_to_cellranger.assert_called_once_with(
                '{}.mtx'.format(counts_prefix),
                '{}.barcodes.txt'.format(counts_prefix),
                '{}.genes.txt'.format(counts_prefix), self.t2g_path,
                cellranger_dir
            )

    def test_count_filter(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.bustools_whitelist') as bustools_whitelist,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
            mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
            mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report'),\
            mock.patch('kb_python.count.import_matrix_as_anndata'):
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_prefix = os.path.join(
                out_dir, UNFILTERED_COUNTS_DIR, COUNTS_PREFIX
            )
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            filtered_counts_prefix = os.path.join(
                out_dir, FILTERED_COUNTS_DIR, COUNTS_PREFIX
            )
            filter_whitelist_path = os.path.join(
                out_dir, FILTER_WHITELIST_FILENAME
            )
            filtered_temp_bus_path = os.path.join(
                temp_dir, BUS_FILTERED_FILENAME
            )
            filtered_bus_path = os.path.join(out_dir, BUS_FILTERED_FILENAME)
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path
            }
            bustools_whitelist.return_value = {
                'whitelist': filter_whitelist_path
            }
            bustools_capture.return_value = {'bus': filtered_temp_bus_path}
            bustools_sort.side_effect = [{
                'bus': bus_s_path
            }, {
                'bus': bus_scs_path
            }, {
                'bus': filtered_bus_path
            }]
            bustools_inspect.return_value = {'inspect': inspect_path}
            bustools_correct.return_value = {'bus': bus_sc_path}
            bustools_count.side_effect = [{
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
            }, {
                'mtx': '{}.mtx'.format(filtered_counts_prefix),
                'genes': '{}.genes.txt'.format(filtered_counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(filtered_counts_prefix),
            }]

            filter_result = {
                'bus_scs': filtered_bus_path,
                'whitelist': filter_whitelist_path,
                'mtx': '{}.mtx'.format(filtered_counts_prefix),
                'genes': '{}.genes.txt'.format(filtered_counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(filtered_counts_prefix),
            }
            filter_with_bustools.return_value = filter_result
            STATS.save.return_value = 'stats'

            self.assertEqual({
                'stats': 'stats',
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'info': info_path,
                    'inspect': inspect_path,
                    'bus_scs': bus_scs_path,
                    'mtx': '{}.mtx'.format(counts_prefix),
                    'genes': '{}.genes.txt'.format(counts_prefix),
                    'barcodes': '{}.barcodes.txt'.format(counts_prefix),
                },
                'filtered': filter_result
            },
                             count.count(
                                 self.index_path,
                                 self.t2g_path,
                                 self.technology,
                                 out_dir,
                                 self.fastqs,
                                 whitelist_path=self.whitelist_path,
                                 filter='bustools',
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory
                             ))
            stream_fastqs.assert_called_once_with(
                self.fastqs, temp_dir=temp_dir
            )
            kallisto_bus.assert_called_once_with(
                self.fastqs,
                self.index_path,
                self.technology,
                out_dir,
                threads=threads
            )
            self.assertEqual(2, bustools_sort.call_count)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path, inspect_path, self.whitelist_path, ecmap_path
            )
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            self.assertEqual(1, bustools_count.call_count)
            bustools_count.assert_called_once_with(
                bus_scs_path,
                counts_prefix,
                self.t2g_path,
                ecmap_path,
                txnames_path,
                tcc=False,
                mm=False
            )
            filter_with_bustools.assert_called_once_with(
                bus_scs_path,
                ecmap_path,
                txnames_path,
                self.t2g_path,
                os.path.join(out_dir, FILTER_WHITELIST_FILENAME),
                os.path.join(out_dir, BUS_FILTERED_FILENAME),
                counts_prefix=os.path.join(
                    out_dir, FILTERED_COUNTS_DIR, COUNTS_PREFIX
                ),
                kite=False,
                temp_dir=temp_dir,
                threads=threads,
                memory=memory,
                loom=False,
                h5ad=False,
                tcc=False,
            )
            convert_matrix.assert_not_called()

    def test_count_without_whitelist(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
            mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report'),\
            mock.patch('kb_python.count.import_matrix_as_anndata'):
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_prefix = os.path.join(
                out_dir, UNFILTERED_COUNTS_DIR, COUNTS_PREFIX
            )
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path
            }
            bustools_sort.side_effect = [{
                'bus': bus_s_path
            }, {
                'bus': bus_scs_path
            }]
            bustools_inspect.return_value = {'inspect': inspect_path}
            copy_or_create_whitelist.return_value = self.whitelist_path
            bustools_correct.return_value = {'bus': bus_sc_path}
            bustools_count.return_value = {
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
            }
            STATS.save.return_value = 'stats'

            self.assertEqual({
                'stats': 'stats',
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'info': info_path,
                    'whitelist': self.whitelist_path,
                    'inspect': inspect_path,
                    'bus_scs': bus_scs_path,
                    'mtx': '{}.mtx'.format(counts_prefix),
                    'genes': '{}.genes.txt'.format(counts_prefix),
                    'barcodes': '{}.barcodes.txt'.format(counts_prefix),
                }
            },
                             count.count(
                                 self.index_path,
                                 self.t2g_path,
                                 self.technology,
                                 out_dir,
                                 self.fastqs,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory
                             ))
            stream_fastqs.assert_called_once_with(
                self.fastqs, temp_dir=temp_dir
            )
            kallisto_bus.assert_called_once_with(
                self.fastqs,
                self.index_path,
                self.technology,
                out_dir,
                threads=threads
            )
            self.assertEqual(bustools_sort.call_count, 2)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path, inspect_path, self.whitelist_path, ecmap_path
            )
            copy_or_create_whitelist.assert_called_once_with(
                self.technology, bus_s_path, out_dir
            )
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            bustools_count.assert_called_once_with(
                bus_scs_path,
                counts_prefix,
                self.t2g_path,
                ecmap_path,
                txnames_path,
                tcc=False,
                mm=False
            )
            convert_matrix.assert_not_called()
            filter_with_bustools.assert_not_called()

    def test_count_kite_convert(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
            mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report'),\
            mock.patch('kb_python.count.import_matrix_as_anndata'):
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_prefix = os.path.join(
                out_dir, UNFILTERED_COUNTS_DIR, FEATURE_PREFIX
            )
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            loom_path = mock.MagicMock()
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path
            }
            bustools_sort.side_effect = [{
                'bus': bus_s_path
            }, {
                'bus': bus_scs_path
            }]
            bustools_inspect.return_value = {'inspect': inspect_path}
            bustools_correct.return_value = {'bus': bus_sc_path}
            bustools_count.return_value = {
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
            }
            convert_matrix.return_value = {'loom': loom_path}
            STATS.save.return_value = 'stats'

            self.assertEqual({
                'stats': 'stats',
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'info': info_path,
                    'inspect': inspect_path,
                    'bus_scs': bus_scs_path,
                    'mtx': '{}.mtx'.format(counts_prefix),
                    'genes': '{}.genes.txt'.format(counts_prefix),
                    'barcodes': '{}.barcodes.txt'.format(counts_prefix),
                    'loom': loom_path,
                }
            },
                             count.count(
                                 self.index_path,
                                 self.t2g_path,
                                 self.technology,
                                 out_dir,
                                 self.fastqs,
                                 kite=True,
                                 whitelist_path=self.whitelist_path,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                                 loom=True,
                             ))
            stream_fastqs.assert_called_once_with(
                self.fastqs, temp_dir=temp_dir
            )
            kallisto_bus.assert_called_once_with(
                self.fastqs,
                self.index_path,
                self.technology,
                out_dir,
                threads=threads
            )
            self.assertEqual(bustools_sort.call_count, 2)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path, inspect_path, self.whitelist_path, ecmap_path
            )
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            bustools_count.assert_called_once_with(
                bus_scs_path,
                counts_prefix,
                self.t2g_path,
                ecmap_path,
                txnames_path,
                tcc=False,
                mm=False
            )
            convert_matrix.assert_called_once_with(
                os.path.join(out_dir, UNFILTERED_COUNTS_DIR),
                '{}.mtx'.format(counts_prefix),
                '{}.barcodes.txt'.format(counts_prefix),
                genes_path='{}.genes.txt'.format(counts_prefix),
                t2g_path=self.t2g_path,
                ec_path=None,
                txnames_path=txnames_path,
                name='feature',
                loom=True,
                h5ad=False,
                tcc=False,
                threads=threads
            )
            filter_with_bustools.assert_not_called()

    def test_count_kite_filter(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.bustools_whitelist') as bustools_whitelist,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
            mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
            mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report'),\
            mock.patch('kb_python.count.import_matrix_as_anndata'):
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_prefix = os.path.join(
                out_dir, UNFILTERED_COUNTS_DIR, FEATURE_PREFIX
            )
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            filtered_counts_prefix = os.path.join(
                out_dir, FILTERED_COUNTS_DIR, FEATURE_PREFIX
            )
            filter_whitelist_path = os.path.join(
                out_dir, FILTER_WHITELIST_FILENAME
            )
            filtered_temp_bus_path = os.path.join(
                temp_dir, BUS_FILTERED_FILENAME
            )
            filtered_bus_path = os.path.join(out_dir, BUS_FILTERED_FILENAME)
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path
            }
            bustools_whitelist.return_value = {
                'whitelist': filter_whitelist_path
            }
            bustools_capture.return_value = {'bus': filtered_temp_bus_path}
            bustools_sort.side_effect = [{
                'bus': bus_s_path
            }, {
                'bus': bus_scs_path
            }, {
                'bus': filtered_bus_path
            }]
            bustools_inspect.return_value = {'inspect': inspect_path}
            bustools_correct.return_value = {'bus': bus_sc_path}
            bustools_count.side_effect = [{
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
            }, {
                'mtx': '{}.mtx'.format(filtered_counts_prefix),
                'genes': '{}.genes.txt'.format(filtered_counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(filtered_counts_prefix),
            }]

            filter_result = {
                'bus_scs': filtered_bus_path,
                'whitelist': filter_whitelist_path,
                'mtx': '{}.mtx'.format(filtered_counts_prefix),
                'genes': '{}.genes.txt'.format(filtered_counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(filtered_counts_prefix),
            }
            filter_with_bustools.return_value = filter_result
            STATS.save.return_value = 'stats'

            self.assertEqual({
                'stats': 'stats',
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'info': info_path,
                    'inspect': inspect_path,
                    'bus_scs': bus_scs_path,
                    'mtx': '{}.mtx'.format(counts_prefix),
                    'genes': '{}.genes.txt'.format(counts_prefix),
                    'barcodes': '{}.barcodes.txt'.format(counts_prefix),
                },
                'filtered': filter_result
            },
                             count.count(
                                 self.index_path,
                                 self.t2g_path,
                                 self.technology,
                                 out_dir,
                                 self.fastqs,
                                 whitelist_path=self.whitelist_path,
                                 kite=True,
                                 filter='bustools',
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory
                             ))
            stream_fastqs.assert_called_once_with(
                self.fastqs, temp_dir=temp_dir
            )
            kallisto_bus.assert_called_once_with(
                self.fastqs,
                self.index_path,
                self.technology,
                out_dir,
                threads=threads
            )
            self.assertEqual(2, bustools_sort.call_count)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path, inspect_path, self.whitelist_path, ecmap_path
            )
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            self.assertEqual(1, bustools_count.call_count)
            bustools_count.assert_called_once_with(
                bus_scs_path,
                counts_prefix,
                self.t2g_path,
                ecmap_path,
                txnames_path,
                tcc=False,
                mm=False
            )
            filter_with_bustools.assert_called_once_with(
                bus_scs_path,
                ecmap_path,
                txnames_path,
                self.t2g_path,
                os.path.join(out_dir, FILTER_WHITELIST_FILENAME),
                os.path.join(out_dir, BUS_FILTERED_FILENAME),
                counts_prefix=os.path.join(
                    out_dir, FILTERED_COUNTS_DIR, FEATURE_PREFIX
                ),
                kite=True,
                temp_dir=temp_dir,
                threads=threads,
                memory=memory,
                loom=False,
                h5ad=False,
                tcc=False,
            )
            convert_matrix.assert_not_called()

    def test_count_kite_FB(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
            mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
            mock.patch('kb_python.count.bustools_project') as bustools_project,\
            mock.patch('kb_python.count.copy_map') as copy_map,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report'),\
            mock.patch('kb_python.count.import_matrix_as_anndata'):
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_prefix = os.path.join(
                out_dir, UNFILTERED_COUNTS_DIR, FEATURE_PREFIX
            )
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            map_path = os.path.join(out_dir, '10xv3_feature_barcode_map.txt')
            bus_s_path = os.path.join(temp_dir, 'output.s.bus')
            bus_sp_path = os.path.join(temp_dir, 'output.s.p.bus')
            bus_sps_path = os.path.join(temp_dir, 'output.s.p.s.bus')
            bus_spsc_path = os.path.join(temp_dir, 'output.s.p.s.c.bus')
            bus_spscs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path
            }
            bustools_sort.side_effect = [{
                'bus': bus_s_path
            }, {
                'bus': bus_sps_path
            }, {
                'bus': bus_spscs_path
            }]
            copy_map.return_value = map_path
            bustools_project.return_value = {'bus': bus_sp_path}
            bustools_inspect.return_value = {'inspect': inspect_path}
            bustools_correct.return_value = {'bus': bus_spsc_path}
            bustools_count.return_value = {
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
            }
            STATS.save.return_value = 'stats'

            self.assertEqual({
                'stats': 'stats',
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'info': info_path,
                    'inspect': inspect_path,
                    'bus_scs': bus_spscs_path,
                    'mtx': '{}.mtx'.format(counts_prefix),
                    'genes': '{}.genes.txt'.format(counts_prefix),
                    'barcodes': '{}.barcodes.txt'.format(counts_prefix),
                }
            },
                             count.count(
                                 self.index_path,
                                 self.t2g_path,
                                 self.technology,
                                 out_dir,
                                 self.fastqs,
                                 kite=True,
                                 FB=True,
                                 whitelist_path=self.whitelist_path,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory
                             ))
            stream_fastqs.assert_called_once_with(
                self.fastqs, temp_dir=temp_dir
            )
            kallisto_bus.assert_called_once_with(
                self.fastqs,
                self.index_path,
                self.technology,
                out_dir,
                threads=threads
            )
            self.assertEqual(3, bustools_sort.call_count)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_sp_path,
                    bus_sps_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_spsc_path,
                    bus_spscs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            copy_map.assert_called_once_with(self.technology, out_dir)
            bustools_project.assert_called_once_with(
                bus_s_path, bus_sp_path, map_path, ecmap_path, txnames_path
            )
            bustools_inspect.assert_called_once_with(
                bus_sps_path, inspect_path, self.whitelist_path, ecmap_path
            )
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_sps_path, bus_spsc_path, self.whitelist_path
            )
            bustools_count.assert_called_once_with(
                bus_spscs_path,
                counts_prefix,
                self.t2g_path,
                ecmap_path,
                txnames_path,
                tcc=False,
                mm=False
            )
            convert_matrix.assert_not_called()
            filter_with_bustools.assert_not_called()

    def test_count_velocity_with_whitelist(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_merge') as bustools_merge,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrices') as convert_matrices,\
            mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report') as render_report,\
            mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata:
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            inspect_cdna_path = os.path.join(
                out_dir, f'inspect.{BUS_CDNA_PREFIX}.json'
            )
            inspect_intron_path = os.path.join(
                out_dir, f'inspect.{BUS_INTRON_PREFIX}.json'
            )
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            cdna_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
            )
            intron_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
            )
            cdna_s_path = os.path.join(
                out_dir, '{}{}'.format(BUS_CDNA_PREFIX, BUS_UNFILTERED_SUFFIX)
            )
            intron_s_path = os.path.join(
                out_dir,
                '{}{}'.format(BUS_INTRON_PREFIX, BUS_UNFILTERED_SUFFIX)
            )
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path
            }
            bustools_sort.side_effect = [{
                'bus': bus_s_path
            }, {
                'bus': bus_scs_path
            }, {
                'bus': cdna_s_path
            }, {
                'bus': intron_s_path
            }]
            bustools_inspect.side_effect = [{
                'inspect': inspect_path
            }, {
                'inspect': inspect_cdna_path
            }, {
                'inspect': inspect_intron_path
            }]
            bustools_capture.side_effect = [{
                'bus': cdna_capture_path
            }, {
                'bus': intron_capture_path
            }]
            bustools_correct.return_value = {'bus': bus_sc_path}
            bustools_count.side_effect = [{
                'mtx':
                    '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ),
            }, {
                'mtx':
                    '{}.mtx'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
            }]
            STATS.save.return_value = 'stats'

            self.assertEqual({
                'stats': 'stats',
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'info': info_path,
                    'inspect': inspect_path,
                    'bus_scs': bus_scs_path,
                    BUS_CDNA_PREFIX: {
                        'bus':
                            cdna_s_path,
                        'mtx':
                            '{}.mtx'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'genes':
                            '{}.genes.txt'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'barcodes':
                            '{}.barcodes.txt'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'inspect':
                            inspect_cdna_path
                    },
                    BUS_INTRON_PREFIX: {
                        'bus':
                            intron_s_path,
                        'mtx':
                            '{}.mtx'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'genes':
                            '{}.genes.txt'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'barcodes':
                            '{}.barcodes.txt'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'inspect':
                            inspect_intron_path
                    }
                }
            },
                             count.count_velocity(
                                 self.index_path,
                                 self.t2g_path,
                                 cdna_t2c_path,
                                 intron_t2c_path,
                                 self.technology,
                                 out_dir,
                                 self.fastqs,
                                 whitelist_path=self.whitelist_path,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory
                             ))
            stream_fastqs.assert_called_once_with(
                self.fastqs, temp_dir=temp_dir
            )
            kallisto_bus.assert_called_once_with(
                self.fastqs,
                self.index_path,
                self.technology,
                out_dir,
                threads=threads
            )
            self.assertEqual(bustools_sort.call_count, 4)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    cdna_capture_path,
                    cdna_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    intron_capture_path,
                    intron_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            self.assertEqual(3, bustools_inspect.call_count)
            bustools_inspect.assert_has_calls([
                call(bus_s_path, inspect_path, self.whitelist_path, ecmap_path),
                call(
                    cdna_s_path, inspect_cdna_path, self.whitelist_path,
                    ecmap_path
                ),
                call(
                    intron_s_path, inspect_intron_path, self.whitelist_path,
                    ecmap_path
                )
            ])
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            self.assertEqual(2, bustools_count.call_count)
            bustools_count.assert_has_calls([
                call(
                    cdna_s_path,
                    os.path.join(counts_dir, BUS_CDNA_PREFIX),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                ),
                call(
                    intron_s_path,
                    os.path.join(counts_dir, BUS_INTRON_PREFIX),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                )
            ])
            filter_with_bustools.assert_not_called()
            convert_matrices.assert_not_called()
            bustools_merge.assert_not_called()

            STATS.start.assert_called_once()
            STATS.end.assert_called_once()
            STATS.to_dict.assert_not_called()
            import_matrix_as_anndata.assert_not_called()
            render_report.assert_not_called()

    def test_count_velocity_cellranger(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_merge') as bustools_merge,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrices') as convert_matrices,\
            mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report') as render_report,\
            mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata,\
            mock.patch('kb_python.count.matrix_to_cellranger') as matrix_to_cellranger:
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
            cellranger_cdna_dir = os.path.join(
                counts_dir, f'{CELLRANGER_DIR}_{BUS_CDNA_PREFIX}'
            )
            cellranger_intron_dir = os.path.join(
                counts_dir, f'{CELLRANGER_DIR}_{BUS_INTRON_PREFIX}'
            )
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            inspect_cdna_path = os.path.join(
                out_dir, f'inspect.{BUS_CDNA_PREFIX}.json'
            )
            inspect_intron_path = os.path.join(
                out_dir, f'inspect.{BUS_INTRON_PREFIX}.json'
            )
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            cdna_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
            )
            intron_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
            )
            cdna_s_path = os.path.join(
                out_dir, '{}{}'.format(BUS_CDNA_PREFIX, BUS_UNFILTERED_SUFFIX)
            )
            intron_s_path = os.path.join(
                out_dir,
                '{}{}'.format(BUS_INTRON_PREFIX, BUS_UNFILTERED_SUFFIX)
            )
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path
            }
            bustools_sort.side_effect = [{
                'bus': bus_s_path
            }, {
                'bus': bus_scs_path
            }, {
                'bus': cdna_s_path
            }, {
                'bus': intron_s_path
            }]
            bustools_inspect.side_effect = [{
                'inspect': inspect_path
            }, {
                'inspect': inspect_cdna_path
            }, {
                'inspect': inspect_intron_path
            }]
            bustools_capture.side_effect = [{
                'bus': cdna_capture_path
            }, {
                'bus': intron_capture_path
            }]
            bustools_correct.return_value = {'bus': bus_sc_path}
            bustools_count.side_effect = [{
                'mtx':
                    '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ),
            }, {
                'mtx':
                    '{}.mtx'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
            }]
            matrix_to_cellranger.side_effect = [{
                'mtx':
                    os.path.join(cellranger_cdna_dir, CELLRANGER_MATRIX),
                'genes':
                    os.path.join(cellranger_cdna_dir, CELLRANGER_GENES),
                'barcodes':
                    os.path.join(cellranger_cdna_dir, CELLRANGER_BARCODES),
            }, {
                'mtx':
                    os.path.join(cellranger_intron_dir, CELLRANGER_MATRIX),
                'genes':
                    os.path.join(cellranger_intron_dir, CELLRANGER_GENES),
                'barcodes':
                    os.path.join(cellranger_intron_dir, CELLRANGER_BARCODES),
            }]
            STATS.save.return_value = 'stats'

            self.assertEqual({
                'stats': 'stats',
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'info': info_path,
                    'inspect': inspect_path,
                    'bus_scs': bus_scs_path,
                    BUS_CDNA_PREFIX: {
                        'bus':
                            cdna_s_path,
                        'mtx':
                            '{}.mtx'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'genes':
                            '{}.genes.txt'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'barcodes':
                            '{}.barcodes.txt'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'inspect':
                            inspect_cdna_path,
                        'cellranger': {
                            'mtx':
                                os.path.join(
                                    cellranger_cdna_dir, CELLRANGER_MATRIX
                                ),
                            'genes':
                                os.path.join(
                                    cellranger_cdna_dir, CELLRANGER_GENES
                                ),
                            'barcodes':
                                os.path.join(
                                    cellranger_cdna_dir, CELLRANGER_BARCODES
                                ),
                        }
                    },
                    BUS_INTRON_PREFIX: {
                        'bus':
                            intron_s_path,
                        'mtx':
                            '{}.mtx'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'genes':
                            '{}.genes.txt'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'barcodes':
                            '{}.barcodes.txt'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'inspect':
                            inspect_intron_path,
                        'cellranger': {
                            'mtx':
                                os.path.join(
                                    cellranger_intron_dir, CELLRANGER_MATRIX
                                ),
                            'genes':
                                os.path.join(
                                    cellranger_intron_dir, CELLRANGER_GENES
                                ),
                            'barcodes':
                                os.path.join(
                                    cellranger_intron_dir, CELLRANGER_BARCODES
                                ),
                        }
                    }
                }
            },
                             count.count_velocity(
                                 self.index_path,
                                 self.t2g_path,
                                 cdna_t2c_path,
                                 intron_t2c_path,
                                 self.technology,
                                 out_dir,
                                 self.fastqs,
                                 whitelist_path=self.whitelist_path,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                                 cellranger=True
                             ))
            stream_fastqs.assert_called_once_with(
                self.fastqs, temp_dir=temp_dir
            )
            kallisto_bus.assert_called_once_with(
                self.fastqs,
                self.index_path,
                self.technology,
                out_dir,
                threads=threads
            )
            self.assertEqual(bustools_sort.call_count, 4)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    cdna_capture_path,
                    cdna_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    intron_capture_path,
                    intron_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            self.assertEqual(3, bustools_inspect.call_count)
            bustools_inspect.assert_has_calls([
                call(bus_s_path, inspect_path, self.whitelist_path, ecmap_path),
                call(
                    cdna_s_path, inspect_cdna_path, self.whitelist_path,
                    ecmap_path
                ),
                call(
                    intron_s_path, inspect_intron_path, self.whitelist_path,
                    ecmap_path
                )
            ])
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            self.assertEqual(2, bustools_count.call_count)
            bustools_count.assert_has_calls([
                call(
                    cdna_s_path,
                    os.path.join(counts_dir, BUS_CDNA_PREFIX),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                ),
                call(
                    intron_s_path,
                    os.path.join(counts_dir, BUS_INTRON_PREFIX),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                )
            ])
            filter_with_bustools.assert_not_called()
            convert_matrices.assert_not_called()
            bustools_merge.assert_not_called()

            STATS.start.assert_called_once()
            STATS.end.assert_called_once()
            STATS.to_dict.assert_not_called()
            import_matrix_as_anndata.assert_not_called()
            render_report.assert_not_called()
            self.assertEqual(2, matrix_to_cellranger.call_count)
            matrix_to_cellranger.assert_has_calls([
                call(
                    '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ), '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ), self.t2g_path, cellranger_cdna_dir
                ),
                call(
                    '{}.mtx'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ), '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ), '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ), self.t2g_path, cellranger_intron_dir
                ),
            ])

    def test_count_velocity_report(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_merge') as bustools_merge,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrices') as convert_matrices,\
            mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report') as render_report,\
            mock.patch('kb_python.count.import_matrix_as_anndata'):
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            inspect_cdna_path = os.path.join(
                out_dir, f'inspect.{BUS_CDNA_PREFIX}.json'
            )
            inspect_intron_path = os.path.join(
                out_dir, f'inspect.{BUS_INTRON_PREFIX}.json'
            )
            report_path = os.path.join(out_dir, REPORT_NOTEBOOK_FILENAME)
            report_cdna_path = os.path.join(
                out_dir, f'report.{BUS_CDNA_PREFIX}.html'
            )
            report_intron_path = os.path.join(
                out_dir, f'report.{BUS_INTRON_PREFIX}.html'
            )
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            cdna_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
            )
            intron_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
            )
            cdna_s_path = os.path.join(
                out_dir, '{}{}'.format(BUS_CDNA_PREFIX, BUS_UNFILTERED_SUFFIX)
            )
            intron_s_path = os.path.join(
                out_dir,
                '{}{}'.format(BUS_INTRON_PREFIX, BUS_UNFILTERED_SUFFIX)
            )
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path
            }
            bustools_sort.side_effect = [{
                'bus': bus_s_path
            }, {
                'bus': bus_scs_path
            }, {
                'bus': cdna_s_path
            }, {
                'bus': intron_s_path
            }]
            bustools_inspect.side_effect = [{
                'inspect': inspect_path
            }, {
                'inspect': inspect_cdna_path
            }, {
                'inspect': inspect_intron_path
            }]
            bustools_capture.side_effect = [{
                'bus': cdna_capture_path
            }, {
                'bus': intron_capture_path
            }]
            bustools_correct.return_value = {'bus': bus_sc_path}
            bustools_count.side_effect = [{
                'mtx':
                    '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ),
            }, {
                'mtx':
                    '{}.mtx'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
            }]
            render_report.side_effect = [{
                'report': report_path
            }, {
                'report': report_cdna_path
            }, {
                'report': report_intron_path
            }]
            STATS.save.return_value = 'stats'

            self.assertEqual({
                'stats': 'stats',
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'info': info_path,
                    'inspect': inspect_path,
                    'report': report_path,
                    'bus_scs': bus_scs_path,
                    BUS_CDNA_PREFIX: {
                        'bus':
                            cdna_s_path,
                        'mtx':
                            '{}.mtx'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'genes':
                            '{}.genes.txt'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'barcodes':
                            '{}.barcodes.txt'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'report':
                            report_cdna_path,
                        'inspect':
                            inspect_cdna_path
                    },
                    BUS_INTRON_PREFIX: {
                        'bus':
                            intron_s_path,
                        'mtx':
                            '{}.mtx'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'genes':
                            '{}.genes.txt'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'barcodes':
                            '{}.barcodes.txt'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'report':
                            report_intron_path,
                        'inspect':
                            inspect_intron_path
                    }
                }
            },
                             count.count_velocity(
                                 self.index_path,
                                 self.t2g_path,
                                 cdna_t2c_path,
                                 intron_t2c_path,
                                 self.technology,
                                 out_dir,
                                 self.fastqs,
                                 whitelist_path=self.whitelist_path,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                                 report=True
                             ))
            stream_fastqs.assert_called_once_with(
                self.fastqs, temp_dir=temp_dir
            )
            kallisto_bus.assert_called_once_with(
                self.fastqs,
                self.index_path,
                self.technology,
                out_dir,
                threads=threads
            )
            self.assertEqual(bustools_sort.call_count, 4)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    cdna_capture_path,
                    cdna_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    intron_capture_path,
                    intron_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            self.assertEqual(3, bustools_inspect.call_count)
            bustools_inspect.assert_has_calls([
                call(bus_s_path, inspect_path, self.whitelist_path, ecmap_path),
                call(
                    cdna_s_path, inspect_cdna_path, self.whitelist_path,
                    ecmap_path
                ),
                call(
                    intron_s_path, inspect_intron_path, self.whitelist_path,
                    ecmap_path
                )
            ])
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            self.assertEqual(2, bustools_count.call_count)
            bustools_count.assert_has_calls([
                call(
                    cdna_s_path,
                    os.path.join(counts_dir, BUS_CDNA_PREFIX),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                ),
                call(
                    intron_s_path,
                    os.path.join(counts_dir, BUS_INTRON_PREFIX),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                )
            ])
            filter_with_bustools.assert_not_called()
            convert_matrices.assert_not_called()
            bustools_merge.assert_not_called()

            STATS.start.assert_called_once()
            STATS.end.assert_called_once()
            self.assertEqual(3, render_report.call_count)
            render_report.assert_has_calls([
                call(
                    'stats',
                    info_path,
                    inspect_path,
                    ANY,
                    ANY,
                    temp_dir=temp_dir
                ),
                call(
                    'stats',
                    info_path,
                    inspect_cdna_path,
                    ANY,
                    ANY,
                    '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ),
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ),
                    self.t2g_path,
                    temp_dir=temp_dir
                ),
                call(
                    'stats',
                    info_path,
                    inspect_intron_path,
                    ANY,
                    ANY,
                    '{}.mtx'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                    self.t2g_path,
                    temp_dir=temp_dir
                )
            ])

    def test_count_velocity_split(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_merge') as bustools_merge,\
            mock.patch('kb_python.count.move_file') as move_file,\
            mock.patch('kb_python.count.get_temporary_filename') as get_temporary_filename,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrices') as convert_matrices,\
            mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report'),\
            mock.patch('kb_python.count.import_matrix_as_anndata'):
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            inspect_cdna_path = os.path.join(
                out_dir, f'inspect.{BUS_CDNA_PREFIX}.json'
            )
            inspect_intron_path = os.path.join(
                out_dir, f'inspect.{BUS_INTRON_PREFIX}.json'
            )
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            cdna_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
            )
            intron_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
            )
            cdna_s_path = os.path.join(
                out_dir, '{}{}'.format(BUS_CDNA_PREFIX, BUS_UNFILTERED_SUFFIX)
            )
            intron_s_path = os.path.join(
                out_dir,
                '{}{}'.format(BUS_INTRON_PREFIX, BUS_UNFILTERED_SUFFIX)
            )
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            stream_fastqs.return_value = self.fastqs
            bustools_merge.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path
            }
            bustools_sort.side_effect = [{
                'bus': 'temp1'
            }, {
                'bus': 'temp2'
            }, {
                'bus': bus_s_path
            }, {
                'bus': bus_scs_path
            }, {
                'bus': cdna_s_path
            }, {
                'bus': intron_s_path
            }]
            bustools_inspect.side_effect = [{
                'inspect': inspect_path
            }, {
                'inspect': inspect_cdna_path
            }, {
                'inspect': inspect_intron_path
            }]
            bustools_capture.side_effect = [{
                'bus': cdna_capture_path
            }, {
                'bus': intron_capture_path
            }]
            bustools_correct.return_value = {'bus': bus_sc_path}
            bustools_count.side_effect = [{
                'mtx':
                    '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ),
            }, {
                'mtx':
                    '{}.mtx'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
            }]
            STATS.save.return_value = 'stats'

            self.assertEqual(
                {
                    'stats': 'stats',
                    'unfiltered': {
                        'bus': bus_path,
                        'ecmap': ecmap_path,
                        'txnames': txnames_path,
                        'info': info_path,
                        'inspect': inspect_path,
                        'bus_scs': bus_scs_path,
                        BUS_CDNA_PREFIX: {
                            'bus':
                                cdna_s_path,
                            'mtx':
                                '{}.mtx'.format(
                                    os.path.join(counts_dir, BUS_CDNA_PREFIX)
                                ),
                            'genes':
                                '{}.genes.txt'.format(
                                    os.path.join(counts_dir, BUS_CDNA_PREFIX)
                                ),
                            'barcodes':
                                '{}.barcodes.txt'.format(
                                    os.path.join(counts_dir, BUS_CDNA_PREFIX)
                                ),
                            'inspect':
                                inspect_cdna_path
                        },
                        BUS_INTRON_PREFIX: {
                            'bus':
                                intron_s_path,
                            'mtx':
                                '{}.mtx'.format(
                                    os.path.join(counts_dir, BUS_INTRON_PREFIX)
                                ),
                            'genes':
                                '{}.genes.txt'.format(
                                    os.path.join(counts_dir, BUS_INTRON_PREFIX)
                                ),
                            'barcodes':
                                '{}.barcodes.txt'.format(
                                    os.path.join(counts_dir, BUS_INTRON_PREFIX)
                                ),
                            'inspect':
                                inspect_intron_path
                        }
                    }
                },
                count.count_velocity(['index.0', 'index.1'],
                                     self.t2g_path,
                                     cdna_t2c_path,
                                     intron_t2c_path,
                                     self.technology,
                                     out_dir,
                                     self.fastqs,
                                     whitelist_path=self.whitelist_path,
                                     temp_dir=temp_dir,
                                     threads=threads,
                                     memory=memory)
            )
            self.assertEqual(2, stream_fastqs.call_count)
            stream_fastqs.assert_has_calls([
                call(self.fastqs, temp_dir=temp_dir),
                call(self.fastqs, temp_dir=temp_dir)
            ])
            self.assertEqual(2, kallisto_bus.call_count)
            kallisto_bus.assert_has_calls([
                call(
                    self.fastqs,
                    'index.0',
                    self.technology,
                    os.path.join(temp_dir, 'bus_part0'),
                    threads=threads,
                    n=True
                ),
                call(
                    self.fastqs,
                    'index.1',
                    self.technology,
                    os.path.join(temp_dir, 'bus_part1'),
                    threads=threads,
                    n=True
                )
            ])
            self.assertEqual(2, move_file.call_count)
            move_file.assert_has_calls([
                call(
                    'temp1', os.path.join(temp_dir, 'bus_part0', 'output.bus')
                ),
                call(
                    'temp2', os.path.join(temp_dir, 'bus_part1', 'output.bus')
                )
            ])
            bustools_merge.assert_called_once_with([
                os.path.join(temp_dir, 'bus_part0'),
                os.path.join(temp_dir, 'bus_part1')
            ],
                                                   out_dir,
                                                   threads=threads)
            self.assertEqual(bustools_sort.call_count, 6)
            bustools_sort.assert_has_calls([
                call(
                    os.path.join(temp_dir, 'bus_part0', 'output.bus'),
                    get_temporary_filename.return_value,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                    flags=True
                ),
                call(
                    os.path.join(temp_dir, 'bus_part1', 'output.bus'),
                    get_temporary_filename.return_value,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                    flags=True
                ),
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    cdna_capture_path,
                    cdna_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    intron_capture_path,
                    intron_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            self.assertEqual(3, bustools_inspect.call_count)
            bustools_inspect.assert_has_calls([
                call(bus_s_path, inspect_path, self.whitelist_path, ecmap_path),
                call(
                    cdna_s_path, inspect_cdna_path, self.whitelist_path,
                    ecmap_path
                ),
                call(
                    intron_s_path, inspect_intron_path, self.whitelist_path,
                    ecmap_path
                )
            ])
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            self.assertEqual(2, bustools_count.call_count)
            bustools_count.assert_has_calls([
                call(
                    cdna_s_path,
                    os.path.join(counts_dir, BUS_CDNA_PREFIX),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                ),
                call(
                    intron_s_path,
                    os.path.join(counts_dir, BUS_INTRON_PREFIX),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                )
            ])
            filter_with_bustools.assert_not_called()
            convert_matrices.assert_not_called()

    def test_count_velocity_convert(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrices') as convert_matrices,\
            mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report'),\
            mock.patch('kb_python.count.import_matrix_as_anndata'):
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            inspect_cdna_path = os.path.join(
                out_dir, f'inspect.{BUS_CDNA_PREFIX}.json'
            )
            inspect_intron_path = os.path.join(
                out_dir, f'inspect.{BUS_INTRON_PREFIX}.json'
            )
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            cdna_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
            )
            intron_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
            )
            cdna_s_path = os.path.join(
                out_dir, '{}{}'.format(BUS_CDNA_PREFIX, BUS_UNFILTERED_SUFFIX)
            )
            intron_s_path = os.path.join(
                out_dir,
                '{}{}'.format(BUS_INTRON_PREFIX, BUS_UNFILTERED_SUFFIX)
            )
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            adata = mock.MagicMock()
            loom_path = os.path.join(counts_dir, '{}.loom'.format(ADATA_PREFIX))
            adata.write_loom.return_value = loom_path
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path
            }
            bustools_sort.side_effect = [{
                'bus': bus_s_path
            }, {
                'bus': bus_scs_path
            }, {
                'bus': cdna_s_path
            }, {
                'bus': intron_s_path
            }]
            bustools_inspect.side_effect = [{
                'inspect': inspect_path
            }, {
                'inspect': inspect_cdna_path
            }, {
                'inspect': inspect_intron_path
            }]
            bustools_capture.side_effect = [{
                'bus': cdna_capture_path
            }, {
                'bus': intron_capture_path
            }]
            bustools_correct.return_value = {'bus': bus_sc_path}
            bustools_count.side_effect = [{
                'mtx':
                    '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ),
            }, {
                'mtx':
                    '{}.mtx'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
            }]
            convert_matrices.return_value = {'loom': loom_path}
            STATS.save.return_value = 'stats'

            self.assertEqual({
                'stats': 'stats',
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'info': info_path,
                    'inspect': inspect_path,
                    'bus_scs': bus_scs_path,
                    'loom': loom_path,
                    BUS_CDNA_PREFIX: {
                        'bus':
                            cdna_s_path,
                        'mtx':
                            '{}.mtx'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'genes':
                            '{}.genes.txt'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'barcodes':
                            '{}.barcodes.txt'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'inspect':
                            inspect_cdna_path
                    },
                    BUS_INTRON_PREFIX: {
                        'bus':
                            intron_s_path,
                        'mtx':
                            '{}.mtx'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'genes':
                            '{}.genes.txt'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'barcodes':
                            '{}.barcodes.txt'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'inspect':
                            inspect_intron_path
                    }
                }
            },
                             count.count_velocity(
                                 self.index_path,
                                 self.t2g_path,
                                 cdna_t2c_path,
                                 intron_t2c_path,
                                 self.technology,
                                 out_dir,
                                 self.fastqs,
                                 whitelist_path=self.whitelist_path,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                                 loom=True
                             ))
            stream_fastqs.assert_called_once_with(
                self.fastqs, temp_dir=temp_dir
            )
            kallisto_bus.assert_called_once_with(
                self.fastqs,
                self.index_path,
                self.technology,
                out_dir,
                threads=threads
            )
            self.assertEqual(bustools_sort.call_count, 4)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    cdna_capture_path,
                    cdna_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    intron_capture_path,
                    intron_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            self.assertEqual(3, bustools_inspect.call_count)
            bustools_inspect.assert_has_calls([
                call(bus_s_path, inspect_path, self.whitelist_path, ecmap_path),
                call(
                    cdna_s_path, inspect_cdna_path, self.whitelist_path,
                    ecmap_path
                ),
                call(
                    intron_s_path, inspect_intron_path, self.whitelist_path,
                    ecmap_path
                )
            ])
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            self.assertEqual(2, bustools_count.call_count)
            bustools_count.assert_has_calls([
                call(
                    cdna_s_path,
                    os.path.join(counts_dir, BUS_CDNA_PREFIX),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                ),
                call(
                    intron_s_path,
                    os.path.join(counts_dir, BUS_INTRON_PREFIX),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                )
            ])
            filter_with_bustools.assert_not_called()
            convert_matrices.assert_called_once_with(
                counts_dir,
                [
                    '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
                    '{}.mtx'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    )
                ],
                [
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ), '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    )
                ],
                genes_paths=[
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ), '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    )
                ],
                t2g_path=self.t2g_path,
                ec_paths=[None, None],
                txnames_path=txnames_path,
                loom=True,
                h5ad=False,
                tcc=False,
                nucleus=False,
            )

    def test_count_velocity_without_whitelist(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrices') as convert_matrices,\
            mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report'),\
            mock.patch('kb_python.count.import_matrix_as_anndata'):
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            inspect_cdna_path = os.path.join(
                out_dir, f'inspect.{BUS_CDNA_PREFIX}.json'
            )
            inspect_intron_path = os.path.join(
                out_dir, f'inspect.{BUS_INTRON_PREFIX}.json'
            )
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            cdna_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
            )
            intron_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
            )
            cdna_s_path = os.path.join(
                out_dir, '{}{}'.format(BUS_CDNA_PREFIX, BUS_UNFILTERED_SUFFIX)
            )
            intron_s_path = os.path.join(
                out_dir,
                '{}{}'.format(BUS_INTRON_PREFIX, BUS_UNFILTERED_SUFFIX)
            )
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path
            }
            bustools_sort.side_effect = [{
                'bus': bus_s_path
            }, {
                'bus': bus_scs_path
            }, {
                'bus': cdna_s_path
            }, {
                'bus': intron_s_path
            }]
            copy_or_create_whitelist.return_value = self.whitelist_path
            bustools_inspect.side_effect = [{
                'inspect': inspect_path
            }, {
                'inspect': inspect_cdna_path
            }, {
                'inspect': inspect_intron_path
            }]
            bustools_capture.side_effect = [{
                'bus': cdna_capture_path
            }, {
                'bus': intron_capture_path
            }]
            bustools_correct.return_value = {'bus': bus_sc_path}
            bustools_count.side_effect = [{
                'mtx':
                    '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ),
            }, {
                'mtx':
                    '{}.mtx'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
            }]
            STATS.save.return_value = 'stats'

            self.assertEqual({
                'stats': 'stats',
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'info': info_path,
                    'inspect': inspect_path,
                    'bus_scs': bus_scs_path,
                    'whitelist': self.whitelist_path,
                    BUS_CDNA_PREFIX: {
                        'bus':
                            cdna_s_path,
                        'mtx':
                            '{}.mtx'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'genes':
                            '{}.genes.txt'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'barcodes':
                            '{}.barcodes.txt'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'inspect':
                            inspect_cdna_path
                    },
                    BUS_INTRON_PREFIX: {
                        'bus':
                            intron_s_path,
                        'mtx':
                            '{}.mtx'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'genes':
                            '{}.genes.txt'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'barcodes':
                            '{}.barcodes.txt'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'inspect':
                            inspect_intron_path
                    }
                }
            },
                             count.count_velocity(
                                 self.index_path,
                                 self.t2g_path,
                                 cdna_t2c_path,
                                 intron_t2c_path,
                                 self.technology,
                                 out_dir,
                                 self.fastqs,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory
                             ))
            stream_fastqs.assert_called_once_with(
                self.fastqs, temp_dir=temp_dir
            )
            kallisto_bus.assert_called_once_with(
                self.fastqs,
                self.index_path,
                self.technology,
                out_dir,
                threads=threads
            )
            self.assertEqual(bustools_sort.call_count, 4)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    cdna_capture_path,
                    cdna_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    intron_capture_path,
                    intron_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            self.assertEqual(3, bustools_inspect.call_count)
            bustools_inspect.assert_has_calls([
                call(bus_s_path, inspect_path, self.whitelist_path, ecmap_path),
                call(
                    cdna_s_path, inspect_cdna_path, self.whitelist_path,
                    ecmap_path
                ),
                call(
                    intron_s_path, inspect_intron_path, self.whitelist_path,
                    ecmap_path
                )
            ])
            copy_or_create_whitelist.assert_called_once_with(
                self.technology, bus_s_path, out_dir
            )
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            self.assertEqual(2, bustools_count.call_count)
            bustools_count.assert_has_calls([
                call(
                    cdna_s_path,
                    os.path.join(counts_dir, BUS_CDNA_PREFIX),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                ),
                call(
                    intron_s_path,
                    os.path.join(counts_dir, BUS_INTRON_PREFIX),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                )
            ])
            filter_with_bustools.assert_not_called()
            convert_matrices.assert_not_called()

    def test_count_velocity_filter(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrices') as convert_matrices,\
            mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report'),\
            mock.patch('kb_python.count.import_matrix_as_anndata'):
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            inspect_cdna_path = os.path.join(
                out_dir, f'inspect.{BUS_CDNA_PREFIX}.json'
            )
            inspect_intron_path = os.path.join(
                out_dir, f'inspect.{BUS_INTRON_PREFIX}.json'
            )
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            cdna_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
            )
            intron_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
            )
            cdna_s_path = os.path.join(
                out_dir, '{}{}'.format(BUS_CDNA_PREFIX, BUS_UNFILTERED_SUFFIX)
            )
            intron_s_path = os.path.join(
                out_dir,
                '{}{}'.format(BUS_INTRON_PREFIX, BUS_UNFILTERED_SUFFIX)
            )
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path
            }
            cdna_filtered_path = mock.MagicMock()
            intron_filtered_path = mock.MagicMock()
            bustools_sort.side_effect = [{
                'bus': bus_s_path
            }, {
                'bus': bus_scs_path
            }, {
                'bus': cdna_s_path
            }, {
                'bus': intron_s_path
            }, {
                'bus': cdna_filtered_path
            }, {
                'bus': intron_filtered_path
            }]
            bustools_inspect.side_effect = [{
                'inspect': inspect_path
            }, {
                'inspect': inspect_cdna_path
            }, {
                'inspect': inspect_intron_path
            }]
            cdna_filtered_capture_path = mock.MagicMock()
            intron_filtered_capture_path = mock.MagicMock()
            bustools_capture.side_effect = [{
                'bus': cdna_capture_path
            }, {
                'bus': intron_capture_path
            }, {
                'bus': cdna_filtered_capture_path
            }, {
                'bus': intron_filtered_capture_path
            }]
            bustools_correct.return_value = {'bus': bus_sc_path}
            bustools_count.side_effect = [{
                'mtx':
                    '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ),
            }, {
                'mtx':
                    '{}.mtx'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
            }, {
                'mtx':
                    '{}.mtx'.format(
                        os.path.join(
                            out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX
                        )
                    ),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(
                            out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX
                        )
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(
                            out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX
                        )
                    ),
            }, {
                'mtx':
                    '{}.mtx'.format(
                        os.path.join(
                            out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
                        )
                    ),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(
                            out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
                        )
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(
                            out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
                        )
                    ),
            }]
            filtered_whitelist_path = os.path.join(
                out_dir, FILTER_WHITELIST_FILENAME
            )
            filtered_bus_path = os.path.join(out_dir, BUS_FILTERED_FILENAME)

            filter_result = {
                'whitelist': filtered_whitelist_path,
                'bus_scs': filtered_bus_path,
            }
            filter_with_bustools.return_value = filter_result
            STATS.save.return_value = 'stats'

            self.assertEqual({
                'stats': 'stats',
                'filtered': {
                    'whitelist': filtered_whitelist_path,
                    'bus_scs': filtered_bus_path,
                    BUS_CDNA_PREFIX: {
                        'bus':
                            cdna_filtered_path,
                        'mtx':
                            '{}.mtx'.format(
                                os.path.join(
                                    out_dir, FILTERED_COUNTS_DIR,
                                    BUS_CDNA_PREFIX
                                )
                            ),
                        'genes':
                            '{}.genes.txt'.format(
                                os.path.join(
                                    out_dir, FILTERED_COUNTS_DIR,
                                    BUS_CDNA_PREFIX
                                )
                            ),
                        'barcodes':
                            '{}.barcodes.txt'.format(
                                os.path.join(
                                    out_dir, FILTERED_COUNTS_DIR,
                                    BUS_CDNA_PREFIX
                                )
                            ),
                    },
                    BUS_INTRON_PREFIX: {
                        'bus':
                            intron_filtered_path,
                        'mtx':
                            '{}.mtx'.format(
                                os.path.join(
                                    out_dir, FILTERED_COUNTS_DIR,
                                    BUS_INTRON_PREFIX
                                )
                            ),
                        'genes':
                            '{}.genes.txt'.format(
                                os.path.join(
                                    out_dir, FILTERED_COUNTS_DIR,
                                    BUS_INTRON_PREFIX
                                )
                            ),
                        'barcodes':
                            '{}.barcodes.txt'.format(
                                os.path.join(
                                    out_dir, FILTERED_COUNTS_DIR,
                                    BUS_INTRON_PREFIX
                                )
                            ),
                    }
                },
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'info': info_path,
                    'inspect': inspect_path,
                    'bus_scs': bus_scs_path,
                    BUS_CDNA_PREFIX: {
                        'bus':
                            cdna_s_path,
                        'mtx':
                            '{}.mtx'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'genes':
                            '{}.genes.txt'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'barcodes':
                            '{}.barcodes.txt'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'inspect':
                            inspect_cdna_path
                    },
                    BUS_INTRON_PREFIX: {
                        'bus':
                            intron_s_path,
                        'mtx':
                            '{}.mtx'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'genes':
                            '{}.genes.txt'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'barcodes':
                            '{}.barcodes.txt'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'inspect':
                            inspect_intron_path
                    }
                }
            },
                             count.count_velocity(
                                 self.index_path,
                                 self.t2g_path,
                                 cdna_t2c_path,
                                 intron_t2c_path,
                                 self.technology,
                                 out_dir,
                                 self.fastqs,
                                 filter='bustools',
                                 whitelist_path=self.whitelist_path,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory
                             ))
            stream_fastqs.assert_called_once_with(
                self.fastqs, temp_dir=temp_dir
            )
            kallisto_bus.assert_called_once_with(
                self.fastqs,
                self.index_path,
                self.technology,
                out_dir,
                threads=threads
            )
            self.assertEqual(bustools_sort.call_count, 6)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    cdna_capture_path,
                    cdna_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    intron_capture_path,
                    intron_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    cdna_filtered_capture_path,
                    os.path.join(
                        out_dir,
                        '{}{}'.format(BUS_CDNA_PREFIX, BUS_FILTERED_SUFFIX)
                    ),
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    intron_filtered_capture_path,
                    os.path.join(
                        out_dir,
                        '{}{}'.format(BUS_INTRON_PREFIX, BUS_FILTERED_SUFFIX)
                    ),
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            self.assertEqual(3, bustools_inspect.call_count)
            bustools_inspect.assert_has_calls([
                call(bus_s_path, inspect_path, self.whitelist_path, ecmap_path),
                call(
                    cdna_s_path, inspect_cdna_path, self.whitelist_path,
                    ecmap_path
                ),
                call(
                    intron_s_path, inspect_intron_path, self.whitelist_path,
                    ecmap_path
                )
            ])
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            self.assertEqual(4, bustools_count.call_count)
            bustools_count.assert_has_calls([
                call(
                    cdna_s_path,
                    os.path.join(counts_dir, BUS_CDNA_PREFIX),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                ),
                call(
                    intron_s_path,
                    os.path.join(counts_dir, BUS_INTRON_PREFIX),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                ),
                call(
                    cdna_filtered_path,
                    os.path.join(out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                ),
                call(
                    intron_filtered_path,
                    os.path.join(
                        out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
                    ),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                )
            ])
            filter_with_bustools.assert_called_once_with(
                bus_scs_path,
                ecmap_path,
                txnames_path,
                self.t2g_path,
                filtered_whitelist_path,
                filtered_bus_path,
                count=False,
            )
            convert_matrices.assert_not_called()

    def test_count_velocity_filter_convert(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrices') as convert_matrices,\
            mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report'),\
            mock.patch('kb_python.count.import_matrix_as_anndata'):
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            inspect_cdna_path = os.path.join(
                out_dir, f'inspect.{BUS_CDNA_PREFIX}.json'
            )
            inspect_intron_path = os.path.join(
                out_dir, f'inspect.{BUS_INTRON_PREFIX}.json'
            )
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            cdna_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
            )
            intron_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
            )
            cdna_s_path = os.path.join(
                out_dir, '{}{}'.format(BUS_CDNA_PREFIX, BUS_UNFILTERED_SUFFIX)
            )
            intron_s_path = os.path.join(
                out_dir,
                '{}{}'.format(BUS_INTRON_PREFIX, BUS_UNFILTERED_SUFFIX)
            )
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path
            }
            cdna_filtered_path = mock.MagicMock()
            intron_filtered_path = mock.MagicMock()
            bustools_sort.side_effect = [{
                'bus': bus_s_path
            }, {
                'bus': bus_scs_path
            }, {
                'bus': cdna_s_path
            }, {
                'bus': intron_s_path
            }, {
                'bus': cdna_filtered_path
            }, {
                'bus': intron_filtered_path
            }]
            bustools_inspect.side_effect = [{
                'inspect': inspect_path
            }, {
                'inspect': inspect_cdna_path
            }, {
                'inspect': inspect_intron_path
            }]
            cdna_filtered_capture_path = mock.MagicMock()
            intron_filtered_capture_path = mock.MagicMock()
            bustools_capture.side_effect = [{
                'bus': cdna_capture_path
            }, {
                'bus': intron_capture_path
            }, {
                'bus': cdna_filtered_capture_path
            }, {
                'bus': intron_filtered_capture_path
            }]
            bustools_correct.return_value = {'bus': bus_sc_path}
            bustools_count.side_effect = [{
                'mtx':
                    '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ),
            }, {
                'mtx':
                    '{}.mtx'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ),
            }, {
                'mtx':
                    '{}.mtx'.format(
                        os.path.join(
                            out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX
                        )
                    ),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(
                            out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX
                        )
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(
                            out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX
                        )
                    ),
            }, {
                'mtx':
                    '{}.mtx'.format(
                        os.path.join(
                            out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
                        )
                    ),
                'genes':
                    '{}.genes.txt'.format(
                        os.path.join(
                            out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
                        )
                    ),
                'barcodes':
                    '{}.barcodes.txt'.format(
                        os.path.join(
                            out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
                        )
                    ),
            }]
            filtered_whitelist_path = os.path.join(
                out_dir, FILTER_WHITELIST_FILENAME
            )
            filtered_bus_path = os.path.join(out_dir, BUS_FILTERED_FILENAME)

            filter_result = {
                'whitelist': filtered_whitelist_path,
                'bus_scs': filtered_bus_path,
            }
            filter_with_bustools.return_value = filter_result
            STATS.save.return_value = 'stats'

            self.assertEqual({
                'stats': 'stats',
                'filtered': {
                    'whitelist': filtered_whitelist_path,
                    'bus_scs': filtered_bus_path,
                    BUS_CDNA_PREFIX: {
                        'bus':
                            cdna_filtered_path,
                        'mtx':
                            '{}.mtx'.format(
                                os.path.join(
                                    out_dir, FILTERED_COUNTS_DIR,
                                    BUS_CDNA_PREFIX
                                )
                            ),
                        'genes':
                            '{}.genes.txt'.format(
                                os.path.join(
                                    out_dir, FILTERED_COUNTS_DIR,
                                    BUS_CDNA_PREFIX
                                )
                            ),
                        'barcodes':
                            '{}.barcodes.txt'.format(
                                os.path.join(
                                    out_dir, FILTERED_COUNTS_DIR,
                                    BUS_CDNA_PREFIX
                                )
                            ),
                    },
                    BUS_INTRON_PREFIX: {
                        'bus':
                            intron_filtered_path,
                        'mtx':
                            '{}.mtx'.format(
                                os.path.join(
                                    out_dir, FILTERED_COUNTS_DIR,
                                    BUS_INTRON_PREFIX
                                )
                            ),
                        'genes':
                            '{}.genes.txt'.format(
                                os.path.join(
                                    out_dir, FILTERED_COUNTS_DIR,
                                    BUS_INTRON_PREFIX
                                )
                            ),
                        'barcodes':
                            '{}.barcodes.txt'.format(
                                os.path.join(
                                    out_dir, FILTERED_COUNTS_DIR,
                                    BUS_INTRON_PREFIX
                                )
                            ),
                    }
                },
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'info': info_path,
                    'inspect': inspect_path,
                    'bus_scs': bus_scs_path,
                    BUS_CDNA_PREFIX: {
                        'bus':
                            cdna_s_path,
                        'mtx':
                            '{}.mtx'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'genes':
                            '{}.genes.txt'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'barcodes':
                            '{}.barcodes.txt'.format(
                                os.path.join(counts_dir, BUS_CDNA_PREFIX)
                            ),
                        'inspect':
                            inspect_cdna_path
                    },
                    BUS_INTRON_PREFIX: {
                        'bus':
                            intron_s_path,
                        'mtx':
                            '{}.mtx'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'genes':
                            '{}.genes.txt'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'barcodes':
                            '{}.barcodes.txt'.format(
                                os.path.join(counts_dir, BUS_INTRON_PREFIX)
                            ),
                        'inspect':
                            inspect_intron_path
                    }
                }
            },
                             count.count_velocity(
                                 self.index_path,
                                 self.t2g_path,
                                 cdna_t2c_path,
                                 intron_t2c_path,
                                 self.technology,
                                 out_dir,
                                 self.fastqs,
                                 filter='bustools',
                                 whitelist_path=self.whitelist_path,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                                 loom=True,
                             ))
            stream_fastqs.assert_called_once_with(
                self.fastqs, temp_dir=temp_dir
            )
            kallisto_bus.assert_called_once_with(
                self.fastqs,
                self.index_path,
                self.technology,
                out_dir,
                threads=threads
            )
            self.assertEqual(bustools_sort.call_count, 6)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    cdna_capture_path,
                    cdna_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    intron_capture_path,
                    intron_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    cdna_filtered_capture_path,
                    os.path.join(
                        out_dir,
                        '{}{}'.format(BUS_CDNA_PREFIX, BUS_FILTERED_SUFFIX)
                    ),
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                ),
                call(
                    intron_filtered_capture_path,
                    os.path.join(
                        out_dir,
                        '{}{}'.format(BUS_INTRON_PREFIX, BUS_FILTERED_SUFFIX)
                    ),
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            self.assertEqual(3, bustools_inspect.call_count)
            bustools_inspect.assert_has_calls([
                call(bus_s_path, inspect_path, self.whitelist_path, ecmap_path),
                call(
                    cdna_s_path, inspect_cdna_path, self.whitelist_path,
                    ecmap_path
                ),
                call(
                    intron_s_path, inspect_intron_path, self.whitelist_path,
                    ecmap_path
                )
            ])
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            self.assertEqual(4, bustools_count.call_count)
            bustools_count.assert_has_calls([
                call(
                    cdna_s_path,
                    os.path.join(counts_dir, BUS_CDNA_PREFIX),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                ),
                call(
                    intron_s_path,
                    os.path.join(counts_dir, BUS_INTRON_PREFIX),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                ),
                call(
                    cdna_filtered_path,
                    os.path.join(out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                ),
                call(
                    intron_filtered_path,
                    os.path.join(
                        out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
                    ),
                    self.t2g_path,
                    ecmap_path,
                    txnames_path,
                    tcc=False,
                    mm=False
                )
            ])
            filter_with_bustools.assert_called_once_with(
                bus_scs_path,
                ecmap_path,
                txnames_path,
                self.t2g_path,
                filtered_whitelist_path,
                filtered_bus_path,
                count=False,
            )
            self.assertEqual(2, convert_matrices.call_count)
            args = [
                call(
                    counts_dir, [
                        '{}.mtx'.format(
                            os.path.join(counts_dir, BUS_CDNA_PREFIX)
                        ), '{}.mtx'.format(
                            os.path.join(counts_dir, BUS_INTRON_PREFIX)
                        )
                    ], [
                        '{}.barcodes.txt'.format(
                            os.path.join(counts_dir, BUS_CDNA_PREFIX)
                        ), '{}.barcodes.txt'.format(
                            os.path.join(counts_dir, BUS_INTRON_PREFIX)
                        )
                    ],
                    genes_paths=[
                        '{}.genes.txt'.format(
                            os.path.join(counts_dir, BUS_CDNA_PREFIX)
                        ), '{}.genes.txt'.format(
                            os.path.join(counts_dir, BUS_INTRON_PREFIX)
                        )
                    ],
                    ec_paths=[None, None],
                    t2g_path=self.t2g_path,
                    txnames_path=txnames_path,
                    loom=True,
                    h5ad=False,
                    tcc=False,
                    nucleus=False
                ),
                call(
                    os.path.join(out_dir, FILTERED_COUNTS_DIR), [
                        '{}.mtx'.format(
                            os.path.join(
                                out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX
                            )
                        ), '{}.mtx'.format(
                            os.path.join(
                                out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
                            )
                        )
                    ], [
                        '{}.barcodes.txt'.format(
                            os.path.join(
                                out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX
                            )
                        ), '{}.barcodes.txt'.format(
                            os.path.join(
                                out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
                            )
                        )
                    ],
                    genes_paths=[
                        '{}.genes.txt'.format(
                            os.path.join(
                                out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX
                            )
                        ), '{}.genes.txt'.format(
                            os.path.join(
                                out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
                            )
                        )
                    ],
                    ec_paths=[None, None],
                    t2g_path=self.t2g_path,
                    txnames_path=txnames_path,
                    loom=True,
                    h5ad=False,
                    tcc=False,
                    nucleus=False
                )
            ]
            self.assertEqual(args[0], convert_matrices.call_args_list[0])
            self.assertEqual(args[1], convert_matrices.call_args_list[1])
