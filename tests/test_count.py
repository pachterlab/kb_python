import os
from unittest import mock, TestCase
from unittest.mock import ANY, call

import kb_python.count as count
from kb_python.constants import (
    ABUNDANCE_FILENAME,
    ABUNDANCE_GENE_FILENAME,
    ABUNDANCE_GENE_TPM_FILENAME,
    ABUNDANCE_TPM_FILENAME,
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
    CAPTURE_FILENAME,
    CELLRANGER_DIR,
    CELLRANGER_BARCODES,
    CELLRANGER_GENES,
    CELLRANGER_MATRIX,
    COUNTS_PREFIX,
    ECMAP_FILENAME,
    FEATURE_PREFIX,
    FILTER_WHITELIST_FILENAME,
    FILTERED_COUNTS_DIR,
    FLD_FILENAME,
    FLENS_FILENAME,
    GENES_FILENAME,
    GENOMEBAM_FILENAME,
    GENOMEBAM_INDEX_FILENAME,
    INSPECT_FILENAME,
    INSPECT_INTERNAL_FILENAME,
    INSPECT_UMI_FILENAME,
    INTERNAL_SUFFIX,
    KALLISTO_INFO_FILENAME,
    KB_INFO_FILENAME,
    REPORT_HTML_FILENAME,
    REPORT_NOTEBOOK_FILENAME,
    SAVED_INDEX_FILENAME,
    TCC_PREFIX,
    TXNAMES_FILENAME,
    UMI_SUFFIX,
    UNFILTERED_COUNTS_DIR,
    UNFILTERED_QUANT_DIR,
    WHITELIST_FILENAME,
)
from tests.mixins import TestMixin


class TestCount(TestMixin, TestCase):

    def setUp(self):
        super(TestCount, self).setUp()
        makedirs_mock = mock.patch('kb_python.count.make_directory')
        makedirs_mock.start()
        self.addCleanup(makedirs_mock.stop)

    # def test_kallisto_bus(self):
    #     out_dir = self.temp_dir
    #     result = count.kallisto_bus(
    #         self.fastqs, self.index_path, self.technology, out_dir, threads=1
    #     )
    #     self.assertEqual({
    #         'bus': os.path.join(out_dir, BUS_FILENAME),
    #         'ecmap': os.path.join(out_dir, ECMAP_FILENAME),
    #         'txnames': os.path.join(out_dir, TXNAMES_FILENAME),
    #         'info': os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    #     }, result)
    #     for key, path in result.items():
    #         self.assertTrue(os.path.exists(path))
    # 
    # def test_kallisto_bus_batch(self):
    #     out_dir = self.temp_dir
    #     result = count.kallisto_bus(
    #         self.smartseq3_single_batch_path,
    #         self.ref_index_path,
    #         'BULK',
    #         out_dir,
    #         threads=1
    #     )
    #     self.assertEqual({
    #         'bus': os.path.join(out_dir, BUS_FILENAME),
    #         'ecmap': os.path.join(out_dir, ECMAP_FILENAME),
    #         'txnames': os.path.join(out_dir, TXNAMES_FILENAME),
    #         'info': os.path.join(out_dir, KALLISTO_INFO_FILENAME),
    #         'saved_index': os.path.join(out_dir, SAVED_INDEX_FILENAME)
    #     }, result)
    # 
    # def test_kallisto_bus_paired(self):
    #     out_dir = self.temp_dir
    #     result = count.kallisto_bus(
    #         self.smartseq3_paired_batch_path,
    #         self.ref_index_path,
    #         'BULK',
    #         out_dir,
    #         threads=1,
    #         paired=True
    #     )
    #     self.assertEqual({
    #         'bus': os.path.join(out_dir, BUS_FILENAME),
    #         'ecmap': os.path.join(out_dir, ECMAP_FILENAME),
    #         'txnames': os.path.join(out_dir, TXNAMES_FILENAME),
    #         'info': os.path.join(out_dir, KALLISTO_INFO_FILENAME),
    #         'saved_index': os.path.join(out_dir, SAVED_INDEX_FILENAME),
    #         'flens': os.path.join(out_dir, FLENS_FILENAME)
    #     }, result)

    # def test_kallisto_quant_tcc_flens(self):
    #     out_dir = self.temp_dir
    #     result = count.kallisto_quant_tcc(
    #         self.quant_mtx_path,
    #         self.saved_index_path,
    #         self.quant_ecmap_path,
    #         self.quant_t2g_path,
    #         out_dir,
    #         flens_path=self.flens_path,
    #         threads=1
    #     )
    #     self.assertEqual({
    #         'genes': os.path.join(out_dir, GENES_FILENAME),
    #         'gene_mtx': os.path.join(out_dir, ABUNDANCE_GENE_FILENAME),
    #         'gene_tpm_mtx': os.path.join(out_dir, ABUNDANCE_GENE_TPM_FILENAME),
    #         'mtx': os.path.join(out_dir, ABUNDANCE_FILENAME),
    #         'tpm_mtx': os.path.join(out_dir, ABUNDANCE_TPM_FILENAME),
    #         'fld': os.path.join(out_dir, FLD_FILENAME),
    #         'txnames': os.path.join(out_dir, TXNAMES_FILENAME),
    #     }, result)
    #     for key, path in result.items():
    #         self.assertTrue(os.path.exists(path))
    # 
    # def test_kallisto_quant_tcc_l_s(self):
    #     out_dir = self.temp_dir
    #     result = count.kallisto_quant_tcc(
    #         self.quant_mtx_path,
    #         self.saved_index_path,
    #         self.quant_ecmap_path,
    #         self.quant_t2g_path,
    #         out_dir,
    #         l=200,
    #         s=20,
    #         threads=1
    #     )
    #     self.assertEqual({
    #         'genes': os.path.join(out_dir, GENES_FILENAME),
    #         'gene_mtx': os.path.join(out_dir, ABUNDANCE_GENE_FILENAME),
    #         'gene_tpm_mtx': os.path.join(out_dir, ABUNDANCE_GENE_TPM_FILENAME),
    #         'mtx': os.path.join(out_dir, ABUNDANCE_FILENAME),
    #         'tpm_mtx': os.path.join(out_dir, ABUNDANCE_TPM_FILENAME),
    #         'fld': os.path.join(out_dir, FLD_FILENAME),
    #         'txnames': os.path.join(out_dir, TXNAMES_FILENAME),
    #     }, result)
    #     for key, path in result.items():
    #         self.assertTrue(os.path.exists(path))

    def test_bustools_project(self):
        out_dir = self.temp_dir
        out_path = os.path.join(out_dir, 'projected.bus')
        result = count.bustools_project(
            self.bus_path, out_path, self.kite_map_path, self.ecmap_path,
            self.txnames_path
        )
        self.assertEqual({'bus': out_path}, result)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_sort(self):
        out_dir = self.temp_dir
        out_path = os.path.join(out_dir, BUS_S_FILENAME)
        result = count.bustools_sort(
            self.bus_path, out_path, threads=1, memory='1G'
        )
        self.assertEqual({'bus': out_path}, result)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_inspect(self):
        out_dir = self.temp_dir
        out_path = os.path.join(out_dir, INSPECT_FILENAME)
        result = count.bustools_inspect(
            self.bus_s_path, out_path, self.whitelist_path, self.ecmap_path
        )
        self.assertEqual({'inspect': out_path}, result)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_correct(self):
        out_dir = self.temp_dir
        out_path = os.path.join(out_dir, BUS_SC_FILENAME)
        result = count.bustools_correct(
            self.bus_s_path, out_path, self.whitelist_path
        )
        self.assertEqual({'bus': out_path}, result)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_count(self):
        out_dir = self.temp_dir
        counts_path = os.path.join(out_dir, COUNTS_PREFIX)
        result = count.bustools_count(
            self.bus_scs_path, counts_path, self.t2g_path, self.ecmap_path,
            self.txnames_path, umi_gene=False
        )
        self.assertEqual({
            'mtx': '{}.mtx'.format(counts_path),
            'genes': '{}.genes.txt'.format(counts_path),
            'barcodes': '{}.barcodes.txt'.format(counts_path),
        }, result)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_count_removes_existing_dir(self):
        out_dir = self.temp_dir
        counts_path = os.path.join(out_dir, COUNTS_PREFIX)
        os.makedirs(counts_path, exist_ok=True)
        result = count.bustools_count(
            self.bus_scs_path, counts_path, self.t2g_path, self.ecmap_path,
            self.txnames_path, umi_gene=False
        )
        self.assertEqual({
            'mtx': '{}.mtx'.format(counts_path),
            'genes': '{}.genes.txt'.format(counts_path),
            'barcodes': '{}.barcodes.txt'.format(counts_path),
        }, result)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_capture(self):
        out_dir = self.temp_dir
        out_path = os.path.join(out_dir, 'capture.bus')
        result = count.bustools_capture(
            self.lamanno_bus_scs_path, out_path, self.lamanno_cdna_t2c_path,
            self.lamanno_txnames_path, self.lamanno_txnames_path
        )
        self.assertEqual({'bus': out_path}, result)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_whitelist(self):
        out_dir = self.temp_dir
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
                name='gene',
                by_name=False,
                loom=True,
                loom_names=['barcode', 'target_name'],
                batch_barcodes_path=None
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
                name='gene',
                by_name=False,
                loom=False,
                loom_names=['barcode', 'target_name'],
                batch_barcodes_path=None
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
                matrix_path,
                barcodes_path,
                ec_path,
                txnames_path,
                threads=8,
                loom=True,
                loom_names=['barcode', 'target_name'],
                batch_barcodes_path=None
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
                    name='gene',
                    by_name=False,
                    loom=True,
                    loom_names=['barcode', 'target_name'],
                    batch_barcodes_path=None
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
                    batch_barcodes_path=None,
                    loom=False,
                    loom_names=['barcode', 'target_name'],
                    t2g_path=t2g_path,
                    name='gene',
                    by_name=False
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
                    matrix_path,
                    barcode_path,
                    ec_path,
                    txnames_path,
                    threads=8,
                    loom=True,
                    loom_names=['barcode', 'target_name'],
                    batch_barcodes_path=None
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
                    name='gene',
                    by_name=False,
                    loom=True,
                    loom_names=['barcode', 'target_name'],
                    batch_barcodes_path=None
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
            temp_dir = self.temp_dir
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
                                 counts_prefix=counts_prefix,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory
                             ))

            bustools_whitelist.assert_called_once_with(
                bus_path, whitelist_path, threshold=None
            )
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
                memory=memory,
            )
            bustools_count.assert_called_once_with(
                sort_path,
                counts_prefix,
                t2g_path,
                ecmap_path,
                txnames_path,
                tcc=False,
                mm=False,
                umi_gene=True,
                em=False,
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
            temp_dir = self.temp_dir
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
                                 counts_prefix=counts_prefix,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                                 loom=True,
                                 loom_names=['barcode', 'target_name'],
                             ))

            bustools_whitelist.assert_called_once_with(
                bus_path, whitelist_path, threshold=None
            )
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
                memory=memory,
            )
            bustools_count.assert_called_once_with(
                sort_path,
                counts_prefix,
                t2g_path,
                ecmap_path,
                txnames_path,
                tcc=False,
                mm=False,
                umi_gene=True,
                em=False,
            )
            convert_matrix.assert_called_once_with(
                counts_dir,
                '{}.mtx'.format(counts_prefix),
                '{}.barcodes.txt'.format(counts_prefix),
                batch_barcodes_path=None,
                genes_path='{}.genes.txt'.format(counts_prefix),
                t2g_path=t2g_path,
                ec_path=None,
                txnames_path=txnames_path,
                name='gene',
                loom=True,
                loom_names=['barcode', 'target_name'],
                h5ad=False,
                by_name=False,
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
            temp_dir = self.temp_dir
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
                                 counts_prefix=counts_prefix,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                                 count=False,
                             ))

            bustools_whitelist.assert_called_once_with(
                bus_path,
                whitelist_path,
                threshold=None,
            )
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
                memory=memory,
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
            temp_dir = self.temp_dir
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
                                 counts_prefix=counts_prefix,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                                 tcc=True,
                             ))

            bustools_whitelist.assert_called_once_with(
                bus_path, whitelist_path, threshold=None
            )
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
                memory=memory,
            )
            bustools_count.assert_called_once_with(
                sort_path,
                counts_prefix,
                t2g_path,
                ecmap_path,
                txnames_path,
                tcc=True,
                mm=False,
                umi_gene=True,
                em=False,
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
            temp_dir = self.temp_dir
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
                                 counts_prefix=counts_prefix,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                                 cellranger=True
                             ))

            bustools_whitelist.assert_called_once_with(
                bus_path, whitelist_path, threshold=None
            )
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
                memory=memory,
            )
            bustools_count.assert_called_once_with(
                sort_path,
                counts_prefix,
                t2g_path,
                ecmap_path,
                txnames_path,
                tcc=False,
                mm=False,
                umi_gene=True,
                em=False,
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
            temp_dir = self.temp_dir
            fastqs = ['path/to/file1.gz', 'path/to/file2.gz']
            stream_file.side_effect = ['FILE 1', 'FILE 2']
            self.assertEqual(
                fastqs, count.stream_fastqs(fastqs, temp_dir=temp_dir)
            )
            stream_file.assert_not_called()

    def test_stream_fastqs_remote(self):
        with mock.patch('kb_python.count.stream_file') as stream_file:
            temp_dir = self.temp_dir
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

    def test_stream_batch(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs:
            stream_fastqs.return_value = ['r1', 'r2']
            path = count.stream_batch(self.smartseq_batch_path, self.temp_dir)
            stream_fastqs.assert_called_once_with([
                'R1.fastq.gz', 'R2.fastq.gz'
            ])
            with open(path, 'r') as f:
                self.assertEqual('1\tr1\tr2\n', f.read())

    def test_copy_or_create_whitelist_provided(self):
        with mock.patch('kb_python.count.copy_whitelist') as copy_whitelist,\
            mock.patch('kb_python.count.bustools_whitelist') as bustools_whitelist:
            out_dir = self.temp_dir
            count.copy_or_create_whitelist(
                self.technology, self.bus_s_path, out_dir
            )
            copy_whitelist.assert_called_once_with(self.technology, out_dir)
            bustools_whitelist.assert_not_called()

    def test_copy_or_create_whitelist_not_provided(self):
        with mock.patch('kb_python.count.copy_whitelist') as copy_whitelist,\
            mock.patch('kb_python.count.bustools_whitelist') as bustools_whitelist:
            out_dir = self.temp_dir
            count.copy_or_create_whitelist(
                'UNSUPPORTED', self.bus_s_path, out_dir
            )
            copy_whitelist.assert_not_called()
            bustools_whitelist.assert_called_once_with(
                self.bus_s_path,
                os.path.join(out_dir, WHITELIST_FILENAME),
            )

    def test_convert_transcripts_to_genes(self):
        with mock.patch('kb_python.count.read_t2g') as read_t2g:
            read_t2g.return_value = {
                'ENST00000003583.12': ('ENSG00000001460.18', 'STPG1')
            }
            genes_path = os.path.join(self.temp_dir, 'genes.txt')
            self.assertEqual(
                genes_path,
                count.convert_transcripts_to_genes(
                    self.smartseq_txnames_path, self.smartseq_t2g_path,
                    genes_path
                )
            )
            with open(genes_path, 'r') as f:
                self.assertEqual(['ENSG00000001460.18', 'ENST00000003912.7'], [
                    line.strip() for line in f if not line.isspace()
                ])

    def test_matrix_to_cellranger(self):
        out_dir = self.temp_dir
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

    def test_write_smartseq3_capture(self):
        capture_path = os.path.join(self.temp_dir, 'capture.txt')
        self.assertEqual(
            capture_path, count.write_smartseq3_capture(capture_path)
        )
        with open(capture_path, 'r') as f:
            self.assertEqual(('T' * 32) + '\n', f.read())

    def test_count_with_whitelist(self):
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
            mock.patch('kb_python.count.render_report') as render_report,\
            mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata:
            out_dir = self.temp_dir
            temp_dir = self.temp_dir
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
            kallisto_bus.assert_called_once()
            self.assertEqual(bustools_sort.call_count, 2)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                    store_num=False
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path,
                inspect_path,
                whitelist_path=self.whitelist_path,
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
                mm=False,
                cm=False,
                umi_gene=True,
                em=False,
                batch_barcodes=False,
            )
            convert_matrix.assert_not_called()
            filter_with_bustools.assert_not_called()

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
            out_dir = self.temp_dir
            temp_dir = self.temp_dir
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
            kallisto_bus.assert_called_once()
            self.assertEqual(bustools_sort.call_count, 2)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                    store_num=False
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path,
                inspect_path,
                whitelist_path=self.whitelist_path,
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
                mm=False,
                cm=False,
                umi_gene=True,
                em=False,
                batch_barcodes=False,
            )
            convert_matrix.assert_not_called()
            filter_with_bustools.assert_not_called()

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
            out_dir = self.temp_dir
            temp_dir = self.temp_dir
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
            kallisto_bus.assert_called_once()
            self.assertEqual(bustools_sort.call_count, 2)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                    store_num=False
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path,
                inspect_path,
                whitelist_path=self.whitelist_path,
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
                mm=False,
                cm=False,
                umi_gene=True,
                em=False,
                batch_barcodes=False
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
                by_name=False,
                tcc=False,
                threads=threads,
                loom_names=['barcode', 'target_name'],
                batch_barcodes_path=None
            )
            filter_with_bustools.assert_not_called()

    def test_count_cellranger(self):
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
            mock.patch('kb_python.count.render_report') as render_report,\
            mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata,\
            mock.patch('kb_python.count.matrix_to_cellranger') as matrix_to_cellranger:
            out_dir = self.temp_dir
            temp_dir = self.temp_dir
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
            kallisto_bus.assert_called_once()
            self.assertEqual(bustools_sort.call_count, 2)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                    store_num=False
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path,
                inspect_path,
                whitelist_path=self.whitelist_path,
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
                mm=False,
                cm=False,
                umi_gene=True,
                em=False,
                batch_barcodes=False,
            )
            convert_matrix.assert_not_called()
            filter_with_bustools.assert_not_called()

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
            out_dir = self.temp_dir
            temp_dir = self.temp_dir
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
            kallisto_bus.assert_called_once()
            self.assertEqual(2, bustools_sort.call_count)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                    store_num=False
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path,
                inspect_path,
                whitelist_path=self.whitelist_path,
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
                mm=False,
                cm=False,
                umi_gene=True,
                em=False,
                batch_barcodes=False
            )
            filter_with_bustools.assert_called_once_with(
                bus_scs_path,
                ecmap_path,
                txnames_path,
                self.t2g_path,
                os.path.join(out_dir, FILTER_WHITELIST_FILENAME),
                os.path.join(out_dir, BUS_FILTERED_FILENAME),
                filter_threshold=None,
                counts_prefix=os.path.join(
                    out_dir, FILTERED_COUNTS_DIR, COUNTS_PREFIX
                ),
                kite=False,
                temp_dir=temp_dir,
                threads=threads,
                memory=memory,
                loom=False,
                h5ad=False,
                by_name=False,
                tcc=False,
                umi_gene=True,
                em=False,
                loom_names=['barcode', 'target_name'],
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
            out_dir = self.temp_dir
            temp_dir = self.temp_dir
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
            kallisto_bus.assert_called_once()
            self.assertEqual(bustools_sort.call_count, 2)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                    store_num=False
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path,
                inspect_path,
                whitelist_path=self.whitelist_path,
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
                mm=False,
                cm=False,
                umi_gene=True,
                em=False,
                batch_barcodes=False,
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
            out_dir = self.temp_dir
            temp_dir = self.temp_dir
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
            kallisto_bus.assert_called_once()
            self.assertEqual(bustools_sort.call_count, 2)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                    store_num=False
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path,
                inspect_path,
                whitelist_path=self.whitelist_path,
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
                mm=False,
                cm=False,
                umi_gene=True,
                em=False,
                batch_barcodes=False,
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
                by_name=False,
                tcc=False,
                threads=threads,
                loom_names=['barcode', 'target_name'],
                batch_barcodes_path=None
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
            out_dir = self.temp_dir
            temp_dir = self.temp_dir
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
            kallisto_bus.assert_called_once()
            self.assertEqual(2, bustools_sort.call_count)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                    store_num=False
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path,
                inspect_path,
                whitelist_path=self.whitelist_path,
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
                mm=False,
                cm=False,
                umi_gene=True,
                em=False,
                batch_barcodes=False,
            )
            filter_with_bustools.assert_called_once_with(
                bus_scs_path,
                ecmap_path,
                txnames_path,
                self.t2g_path,
                os.path.join(out_dir, FILTER_WHITELIST_FILENAME),
                os.path.join(out_dir, BUS_FILTERED_FILENAME),
                filter_threshold=None,
                counts_prefix=os.path.join(
                    out_dir, FILTERED_COUNTS_DIR, FEATURE_PREFIX
                ),
                kite=True,
                temp_dir=temp_dir,
                threads=threads,
                memory=memory,
                loom=False,
                h5ad=False,
                by_name=False,
                tcc=False,
                umi_gene=True,
                em=False,
                loom_names=['barcode', 'target_name'],
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
            mock.patch('kb_python.count.create_10x_feature_barcode_map') as create_10x_feature_barcode_map,\
            mock.patch('kb_python.count.STATS') as STATS,\
            mock.patch('kb_python.count.render_report'),\
            mock.patch('kb_python.count.import_matrix_as_anndata'):
            out_dir = self.temp_dir
            temp_dir = self.temp_dir
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
            map_path = os.path.join(out_dir, '10x_feature_barcode_map.txt')
            bus_s_path = os.path.join(temp_dir, 'output.s.bus')
            bus_sc_path = os.path.join(temp_dir, 'output.s.c.bus')
            bus_scs_path = os.path.join(temp_dir, 'output.s.c.s.bus')
            bus_scsp_path = os.path.join(temp_dir, 'output.s.c.s.p.bus')
            bus_scsps_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
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
                'bus': bus_scsps_path
            }]
            create_10x_feature_barcode_map.return_value = map_path
            bustools_project.return_value = {'bus': bus_scsp_path}
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
                    'bus_scs': bus_scsps_path,
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
            kallisto_bus.assert_called_once()
            self.assertEqual(3, bustools_sort.call_count)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                    store_num=False
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                ),
                call(
                    bus_scsp_path,
                    bus_scsps_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                )
            ])
            create_10x_feature_barcode_map.assert_called_once_with(map_path)
            bustools_project.assert_called_once_with(
                bus_scs_path, bus_scsp_path, map_path, ecmap_path, txnames_path
            )
            bustools_inspect.assert_called_once_with(
                bus_s_path,
                inspect_path,
                whitelist_path=self.whitelist_path,
            )
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            bustools_count.assert_called_once_with(
                bus_scsps_path,
                counts_prefix,
                self.t2g_path,
                ecmap_path,
                txnames_path,
                tcc=False,
                mm=False,
                cm=False,
                umi_gene=True,
                em=False,
                batch_barcodes=False,
            )
            convert_matrix.assert_not_called()
            filter_with_bustools.assert_not_called()

    def test_count_bulk_multi_paired(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.stream_batch') as stream_batch,\
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
            out_dir = self.temp_dir
            temp_dir = self.temp_dir
            counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
            counts_prefix = os.path.join(counts_dir, COUNTS_PREFIX)
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            flens_path = os.path.join(out_dir, FLENS_FILENAME)
            saved_index_path = os.path.join(out_dir, SAVED_INDEX_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            fastqs = [
                self.smartseq3_1_i1_fastq_path, self.smartseq3_1_i2_fastq_path,
                self.smartseq3_1_R1_fastq_path, self.smartseq3_1_R2_fastq_path,
                self.smartseq3_2_i1_fastq_path, self.smartseq3_2_i2_fastq_path,
                self.smartseq3_2_R1_fastq_path, self.smartseq3_2_R2_fastq_path
            ]
            stream_fastqs.return_value = fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path,
                'flens': flens_path,
                'saved_index': saved_index_path
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
                    'flens': flens_path,
                    'saved_index': saved_index_path,
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
                                 'SMARTSEQ2',
                                 out_dir,
                                 fastqs,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                                 paired=True,
                                 h5ad=True
                             ))
            stream_fastqs.assert_called_once_with(fastqs, temp_dir=temp_dir)
            kallisto_bus.assert_called_once()
            self.assertEqual(bustools_sort.call_count, 2)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                    store_num=False
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path,
                inspect_path,
                whitelist_path=self.whitelist_path,
            )
            copy_or_create_whitelist.assert_called_once_with(
                'SMARTSEQ2', bus_s_path, out_dir
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
                mm=False,
                cm=True,
                umi_gene=True,
                em=False,
                batch_barcodes=False,
            )
            convert_matrix.assert_called_once_with(
                counts_dir,
                f'{counts_prefix}.mtx',
                f'{counts_prefix}.barcodes.txt',
                genes_path=f'{counts_prefix}.genes.txt',
                t2g_path=self.t2g_path,
                ec_path=None,
                txnames_path=txnames_path,
                name='gene',
                loom=False,
                h5ad=True,
                by_name=False,
                tcc=False,
                threads=threads,
                loom_names=['barcode', 'target_name'],
                batch_barcodes_path=None
            )
            filter_with_bustools.assert_not_called()
            stream_batch.assert_not_called()

    def test_count_bulk_multi_single(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.stream_batch') as stream_batch,\
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
            out_dir = self.temp_dir
            temp_dir = self.temp_dir
            counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
            counts_prefix = os.path.join(counts_dir, COUNTS_PREFIX)
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
            flens_path = os.path.join(out_dir, FLENS_FILENAME)
            saved_index_path = os.path.join(out_dir, SAVED_INDEX_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
            fastqs = [
                self.smartseq3_1_i1_fastq_path,
                self.smartseq3_1_i2_fastq_path,
                self.smartseq3_1_R1_fastq_path,
                self.smartseq3_2_i1_fastq_path,
                self.smartseq3_2_i2_fastq_path,
                self.smartseq3_2_R1_fastq_path,
            ]
            stream_fastqs.return_value = fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'info': info_path,
                'flens': flens_path,
                'saved_index': saved_index_path,
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
                    'flens': flens_path,
                    'saved_index': saved_index_path,
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
                                 'SMARTSEQ2',
                                 out_dir,
                                 fastqs,
                                 temp_dir=temp_dir,
                                 threads=threads,
                                 memory=memory,
                                 paired=False,
                                 h5ad=True
                             ))
            stream_fastqs.assert_called_once_with(fastqs, temp_dir=temp_dir)
            kallisto_bus.assert_called_once()
            self.assertEqual(bustools_sort.call_count, 2)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                    store_num=False
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path,
                inspect_path,
                whitelist_path=self.whitelist_path,
            )
            copy_or_create_whitelist.assert_called_once_with(
                'SMARTSEQ2', bus_s_path, out_dir
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
                mm=False,
                cm=True,
                umi_gene=True,
                em=False,
                batch_barcodes=False,
            )
            convert_matrix.assert_called_once_with(
                counts_dir,
                f'{counts_prefix}.mtx',
                f'{counts_prefix}.barcodes.txt',
                genes_path=f'{counts_prefix}.genes.txt',
                t2g_path=self.t2g_path,
                ec_path=None,
                txnames_path=txnames_path,
                name='gene',
                loom=False,
                h5ad=True,
                by_name=False,
                tcc=False,
                threads=threads,
                loom_names=['barcode', 'target_name'],
                batch_barcodes_path=None
            )
            filter_with_bustools.assert_not_called()
            stream_batch.assert_not_called()

    # def test_count_bulk_demux_paired(self):
    #     with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
    #         mock.patch('kb_python.count.stream_batch') as stream_batch,\
    #         mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
    #         mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
    #         mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
    #         mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
    #         mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
    #         mock.patch('kb_python.count.bustools_count') as bustools_count,\
    #         mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
    #         mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
    #         mock.patch('kb_python.count.STATS') as STATS,\
    #         mock.patch('kb_python.count.render_report'),\
    #         mock.patch('kb_python.count.import_matrix_as_anndata'):
    #         out_dir = self.temp_dir
    #         temp_dir = self.temp_dir
    #         counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
    #         counts_prefix = os.path.join(counts_dir, COUNTS_PREFIX)
    #         threads = 99999
    #         memory = 'TEST'
    #         bus_path = os.path.join(out_dir, BUS_FILENAME)
    #         ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
    #         txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
    #         inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
    #         info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    #         flens_path = os.path.join(out_dir, FLENS_FILENAME)
    #         saved_index_path = os.path.join(out_dir, SAVED_INDEX_FILENAME)
    #         bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
    #         batch_path = self.smartseq3_paired_batch_path
    #         stream_batch.return_value = batch_path
    #         kallisto_bus.return_value = {
    #             'bus': bus_path,
    #             'ecmap': ecmap_path,
    #             'txnames': txnames_path,
    #             'info': info_path,
    #             'flens': flens_path,
    #             'saved_index': saved_index_path
    #         }
    #         bustools_sort.return_value = {'bus': bus_s_path}
    #         bustools_inspect.return_value = {'inspect': inspect_path}
    #         bustools_count.return_value = {
    #             'mtx': '{}.mtx'.format(counts_prefix),
    #             'genes': '{}.genes.txt'.format(counts_prefix),
    #             'barcodes': '{}.barcodes.txt'.format(counts_prefix),
    #         }
    #         STATS.save.return_value = 'stats'
    #         self.assertEqual({
    #             'stats': 'stats',
    #             'unfiltered': {
    #                 'bus': bus_path,
    #                 'ecmap': ecmap_path,
    #                 'txnames': txnames_path,
    #                 'info': info_path,
    #                 'flens': flens_path,
    #                 'saved_index': saved_index_path,
    #                 'inspect': inspect_path,
    #                 'mtx': '{}.mtx'.format(counts_prefix),
    #                 'genes': '{}.genes.txt'.format(counts_prefix),
    #                 'barcodes': '{}.barcodes.txt'.format(counts_prefix),
    #             }
    #         },
    #                          count.count(
    #                              self.index_path,
    #                              self.t2g_path,
    #                              'SMARTSEQ2',
    #                              out_dir,
    #                              batch_path,
    #                              temp_dir=temp_dir,
    #                              threads=threads,
    #                              memory=memory,
    #                              paired=True,
    #                              h5ad=True
    #                          ))
    #         stream_batch.assert_called_once_with(batch_path, temp_dir=temp_dir)
    #         stream_fastqs.assert_not_called()
    #         kallisto_bus.assert_called_once_with(
    #             batch_path,
    #             self.index_path,
    #             'BULK',
    #             out_dir,
    #             threads=threads,
    #             paired=True,
    #             genomebam=False,
    #             strand=None,
    #             gtf_path=None,
    #             chromosomes_path=None,
    #         )
    #         bustools_sort.assert_called_once_with(
    #             bus_path,
    #             bus_s_path,
    #             temp_dir=temp_dir,
    #             threads=threads,
    #             memory=memory,
    #             store_num=False
    #         )
    #         bustools_inspect.assert_called_once_with(
    #             bus_s_path,
    #             inspect_path,
    #             whitelist_path=None,
    #         )
    #         copy_or_create_whitelist.assert_not_called()
    #         bustools_correct.assert_not_called()
    #         bustools_count.assert_called_once_with(
    #             bus_s_path,
    #             counts_prefix,
    #             self.t2g_path,
    #             ecmap_path,
    #             txnames_path,
    #             tcc=False,
    #             mm=False,
    #             cm=True,
    #             umi_gene=True,
    #             em=False,
    #         )
    #         convert_matrix.assert_called_once_with(
    #             counts_dir,
    #             f'{counts_prefix}.mtx',
    #             f'{counts_prefix}.barcodes.txt',
    #             genes_path=f'{counts_prefix}.genes.txt',
    #             t2g_path=self.t2g_path,
    #             ec_path=None,
    #             txnames_path=txnames_path,
    #             name='gene',
    #             loom=False,
    #             h5ad=True,
    #             by_name=False,
    #             tcc=False,
    #             threads=threads
    #         )
    #         filter_with_bustools.assert_not_called()

    # def test_count_bulk_demux_single(self):
    #     with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
    #         mock.patch('kb_python.count.stream_batch') as stream_batch,\
    #         mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
    #         mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
    #         mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
    #         mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
    #         mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
    #         mock.patch('kb_python.count.bustools_count') as bustools_count,\
    #         mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
    #         mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
    #         mock.patch('kb_python.count.STATS') as STATS,\
    #         mock.patch('kb_python.count.render_report'),\
    #         mock.patch('kb_python.count.import_matrix_as_anndata'):
    #         out_dir = self.temp_dir
    #         temp_dir = self.temp_dir
    #         counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
    #         counts_prefix = os.path.join(counts_dir, COUNTS_PREFIX)
    #         threads = 99999
    #         memory = 'TEST'
    #         bus_path = os.path.join(out_dir, BUS_FILENAME)
    #         ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
    #         txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
    #         inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
    #         info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    #         saved_index_path = os.path.join(out_dir, SAVED_INDEX_FILENAME)
    #         bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
    #         batch_path = self.smartseq3_paired_batch_path
    #         stream_batch.return_value = batch_path
    #         kallisto_bus.return_value = {
    #             'bus': bus_path,
    #             'ecmap': ecmap_path,
    #             'txnames': txnames_path,
    #             'info': info_path,
    #             'saved_index': saved_index_path
    #         }
    #         bustools_sort.return_value = {'bus': bus_s_path}
    #         bustools_inspect.return_value = {'inspect': inspect_path}
    #         bustools_count.return_value = {
    #             'mtx': '{}.mtx'.format(counts_prefix),
    #             'genes': '{}.genes.txt'.format(counts_prefix),
    #             'barcodes': '{}.barcodes.txt'.format(counts_prefix),
    #         }
    #         STATS.save.return_value = 'stats'
    #         self.maxDiff = None
    #         self.assertEqual({
    #             'stats': 'stats',
    #             'unfiltered': {
    #                 'bus': bus_path,
    #                 'ecmap': ecmap_path,
    #                 'txnames': txnames_path,
    #                 'info': info_path,
    #                 'saved_index': saved_index_path,
    #                 'inspect': inspect_path,
    #                 'mtx': '{}.mtx'.format(counts_prefix),
    #                 'genes': '{}.genes.txt'.format(counts_prefix),
    #                 'barcodes': '{}.barcodes.txt'.format(counts_prefix),
    #             }
    #         },
    #                          count.count(
    #                              self.index_path,
    #                              self.t2g_path,
    #                              'SMARTSEQ2',
    #                              out_dir,
    #                              batch_path,
    #                              temp_dir=temp_dir,
    #                              threads=threads,
    #                              memory=memory,
    #                              paired=False,
    #                              h5ad=True
    #                          ))
    #         stream_batch.assert_called_once_with(batch_path, temp_dir=temp_dir)
    #         stream_fastqs.assert_not_called()
    #         kallisto_bus.assert_called_once_with(
    #             batch_path,
    #             self.index_path,
    #             'BULK',
    #             out_dir,
    #             threads=threads,
    #             paired=False,
    #             genomebam=False,
    #             strand=None,
    #             gtf_path=None,
    #             chromosomes_path=None,
    #         )
    #         bustools_sort.assert_called_once_with(
    #             bus_path,
    #             bus_s_path,
    #             temp_dir=temp_dir,
    #             threads=threads,
    #             memory=memory,
    #             store_num=False
    #         )
    #         bustools_inspect.assert_called_once_with(
    #             bus_s_path,
    #             inspect_path,
    #             whitelist_path=None,
    #         )
    #         copy_or_create_whitelist.assert_not_called()
    #         bustools_correct.assert_not_called()
    #         bustools_count.assert_called_once_with(
    #             bus_s_path,
    #             counts_prefix,
    #             self.t2g_path,
    #             ecmap_path,
    #             txnames_path,
    #             tcc=False,
    #             mm=False,
    #             cm=True,
    #             umi_gene=True,
    #             em=False,
    #         )
    #         convert_matrix.assert_called_once_with(
    #             counts_dir,
    #             f'{counts_prefix}.mtx',
    #             f'{counts_prefix}.barcodes.txt',
    #             genes_path=f'{counts_prefix}.genes.txt',
    #             t2g_path=self.t2g_path,
    #             ec_path=None,
    #             txnames_path=txnames_path,
    #             name='gene',
    #             loom=False,
    #             h5ad=True,
    #             by_name=False,
    #             tcc=False,
    #             threads=threads
    #         )
    #         filter_with_bustools.assert_not_called()
    # 
    # def test_count_bulk_demux_paired_tcc(self):
    #     with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
    #         mock.patch('kb_python.count.stream_batch') as stream_batch,\
    #         mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
    #         mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
    #         mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
    #         mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
    #         mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
    #         mock.patch('kb_python.count.bustools_count') as bustools_count,\
    #         mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
    #         mock.patch('kb_python.count.kallisto_quant_tcc') as kallisto_quant_tcc,\
    #         mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
    #         mock.patch('kb_python.count.STATS') as STATS,\
    #         mock.patch('kb_python.count.render_report'),\
    #         mock.patch('kb_python.count.import_matrix_as_anndata'):
    #         out_dir = self.temp_dir
    #         temp_dir = self.temp_dir
    #         counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
    #         counts_prefix = os.path.join(counts_dir, TCC_PREFIX)
    #         quant_dir = os.path.join(out_dir, UNFILTERED_QUANT_DIR)
    #         threads = 99999
    #         memory = 'TEST'
    #         bus_path = os.path.join(out_dir, BUS_FILENAME)
    #         ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
    #         txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
    #         inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
    #         info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    #         flens_path = os.path.join(out_dir, FLENS_FILENAME)
    #         saved_index_path = os.path.join(out_dir, SAVED_INDEX_FILENAME)
    #         bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
    #         batch_path = self.smartseq3_paired_batch_path
    #         stream_batch.return_value = batch_path
    #         kallisto_bus.return_value = {
    #             'bus': bus_path,
    #             'ecmap': ecmap_path,
    #             'txnames': txnames_path,
    #             'info': info_path,
    #             'flens': flens_path,
    #             'saved_index': saved_index_path
    #         }
    #         bustools_sort.return_value = {'bus': bus_s_path}
    #         bustools_inspect.return_value = {'inspect': inspect_path}
    #         bustools_count.return_value = {
    #             'mtx': '{}.mtx'.format(counts_prefix),
    #             'ec': '{}.ec.txt'.format(counts_prefix),
    #             'barcodes': '{}.barcodes.txt'.format(counts_prefix),
    #         }
    #         kallisto_quant_tcc.return_value = {
    #             'genes':
    #                 os.path.join(quant_dir, GENES_FILENAME),
    #             'gene_mtx':
    #                 os.path.join(quant_dir, ABUNDANCE_GENE_FILENAME),
    #             'gene_tpm_mtx':
    #                 os.path.join(quant_dir, ABUNDANCE_GENE_TPM_FILENAME),
    #             'mtx':
    #                 os.path.join(quant_dir, ABUNDANCE_FILENAME),
    #             'tpm_mtx':
    #                 os.path.join(quant_dir, ABUNDANCE_TPM_FILENAME),
    #             'fld':
    #                 os.path.join(quant_dir, FLD_FILENAME),
    #             'txnames':
    #                 os.path.join(quant_dir, TXNAMES_FILENAME),
    #         }
    #         STATS.save.return_value = 'stats'
    #         self.assertEqual({
    #             'stats': 'stats',
    #             'unfiltered': {
    #                 'bus':
    #                     bus_path,
    #                 'ecmap':
    #                     ecmap_path,
    #                 'ec':
    #                     f'{counts_prefix}.ec.txt',
    #                 'info':
    #                     info_path,
    #                 'flens':
    #                     flens_path,
    #                 'saved_index':
    #                     saved_index_path,
    #                 'inspect':
    #                     inspect_path,
    #                 'genes':
    #                     os.path.join(quant_dir, GENES_FILENAME),
    #                 'gene_mtx':
    #                     os.path.join(quant_dir, ABUNDANCE_GENE_FILENAME),
    #                 'gene_tpm_mtx':
    #                     os.path.join(quant_dir, ABUNDANCE_GENE_TPM_FILENAME),
    #                 'mtx':
    #                     os.path.join(quant_dir, ABUNDANCE_FILENAME),
    #                 'tpm_mtx':
    #                     os.path.join(quant_dir, ABUNDANCE_TPM_FILENAME),
    #                 'fld':
    #                     os.path.join(quant_dir, FLD_FILENAME),
    #                 'txnames':
    #                     os.path.join(quant_dir, TXNAMES_FILENAME),
    #                 'barcodes':
    #                     '{}.barcodes.txt'.format(counts_prefix),
    #             }
    #         },
    #                          count.count(
    #                              self.index_path,
    #                              self.t2g_path,
    #                              'SMARTSEQ2',
    #                              out_dir,
    #                              batch_path,
    #                              temp_dir=temp_dir,
    #                              threads=threads,
    #                              memory=memory,
    #                              paired=True,
    #                              h5ad=True,
    #                              tcc=True
    #                          ))
    #         stream_batch.assert_called_once_with(batch_path, temp_dir=temp_dir)
    #         stream_fastqs.assert_not_called()
    #         kallisto_bus.assert_called_once_with(
    #             batch_path,
    #             self.index_path,
    #             'BULK',
    #             out_dir,
    #             threads=threads,
    #             paired=True,
    #             genomebam=False,
    #             strand=None,
    #             gtf_path=None,
    #             chromosomes_path=None,
    #         )
    #         bustools_sort.assert_called_once_with(
    #             bus_path,
    #             bus_s_path,
    #             temp_dir=temp_dir,
    #             threads=threads,
    #             memory=memory,
    #             store_num=False
    #         )
    #         bustools_inspect.assert_called_once_with(
    #             bus_s_path,
    #             inspect_path,
    #             whitelist_path=None,
    #         )
    #         copy_or_create_whitelist.assert_not_called()
    #         bustools_correct.assert_not_called()
    #         bustools_count.assert_called_once_with(
    #             bus_s_path,
    #             counts_prefix,
    #             self.t2g_path,
    #             ecmap_path,
    #             txnames_path,
    #             tcc=True,
    #             mm=True,
    #             cm=True,
    #             umi_gene=True,
    #             em=False,
    #         )
    #         kallisto_quant_tcc.assert_called_once_with(
    #             f'{counts_prefix}.mtx',
    #             saved_index_path,
    #             ecmap_path,
    #             self.t2g_path,
    #             quant_dir,
    #             flens_path=flens_path,
    #             l=None,
    #             s=None,
    #             threads=threads
    #         )
    #         convert_matrix.assert_called_once_with(
    #             quant_dir,
    #             os.path.join(quant_dir, ABUNDANCE_FILENAME),
    #             f'{counts_prefix}.barcodes.txt',
    #             genes_path=os.path.join(quant_dir, TXNAMES_FILENAME),
    #             t2g_path=self.t2g_path,
    #             ec_path=f'{counts_prefix}.ec.txt',
    #             txnames_path=os.path.join(out_dir, TXNAMES_FILENAME),
    #             name='transcript',
    #             loom=False,
    #             h5ad=True,
    #             by_name=False,
    #             tcc=False,
    #             threads=threads
    #         )
    #         filter_with_bustools.assert_not_called()

    # def test_count_bulk_demux_single_tcc(self):
    #     with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
    #         mock.patch('kb_python.count.stream_batch') as stream_batch,\
    #         mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
    #         mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
    #         mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
    #         mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
    #         mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
    #         mock.patch('kb_python.count.bustools_count') as bustools_count,\
    #         mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
    #         mock.patch('kb_python.count.kallisto_quant_tcc') as kallisto_quant_tcc,\
    #         mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
    #         mock.patch('kb_python.count.STATS') as STATS,\
    #         mock.patch('kb_python.count.render_report'),\
    #         mock.patch('kb_python.count.import_matrix_as_anndata'):
    #         out_dir = self.temp_dir
    #         temp_dir = self.temp_dir
    #         counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
    #         counts_prefix = os.path.join(counts_dir, TCC_PREFIX)
    #         quant_dir = os.path.join(out_dir, UNFILTERED_QUANT_DIR)
    #         threads = 99999
    #         memory = 'TEST'
    #         bus_path = os.path.join(out_dir, BUS_FILENAME)
    #         ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
    #         txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
    #         inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
    #         info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    #         flens_path = os.path.join(out_dir, FLENS_FILENAME)
    #         saved_index_path = os.path.join(out_dir, SAVED_INDEX_FILENAME)
    #         bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
    #         batch_path = self.smartseq3_paired_batch_path
    #         stream_batch.return_value = batch_path
    #         kallisto_bus.return_value = {
    #             'bus': bus_path,
    #             'ecmap': ecmap_path,
    #             'txnames': txnames_path,
    #             'info': info_path,
    #             'flens': flens_path,
    #             'saved_index': saved_index_path
    #         }
    #         bustools_sort.return_value = {'bus': bus_s_path}
    #         bustools_inspect.return_value = {'inspect': inspect_path}
    #         bustools_count.return_value = {
    #             'mtx': '{}.mtx'.format(counts_prefix),
    #             'ec': '{}.ec.txt'.format(counts_prefix),
    #             'barcodes': '{}.barcodes.txt'.format(counts_prefix),
    #         }
    #         kallisto_quant_tcc.return_value = {
    #             'genes':
    #                 os.path.join(quant_dir, GENES_FILENAME),
    #             'gene_mtx':
    #                 os.path.join(quant_dir, ABUNDANCE_GENE_FILENAME),
    #             'gene_tpm_mtx':
    #                 os.path.join(quant_dir, ABUNDANCE_GENE_TPM_FILENAME),
    #             'mtx':
    #                 os.path.join(quant_dir, ABUNDANCE_FILENAME),
    #             'tpm_mtx':
    #                 os.path.join(quant_dir, ABUNDANCE_TPM_FILENAME),
    #             'fld':
    #                 os.path.join(quant_dir, FLD_FILENAME),
    #             'txnames':
    #                 os.path.join(quant_dir, TXNAMES_FILENAME),
    #         }
    #         STATS.save.return_value = 'stats'
    #         self.assertEqual({
    #             'stats': 'stats',
    #             'unfiltered': {
    #                 'bus':
    #                     bus_path,
    #                 'ecmap':
    #                     ecmap_path,
    #                 'ec':
    #                     f'{counts_prefix}.ec.txt',
    #                 'info':
    #                     info_path,
    #                 'flens':
    #                     flens_path,
    #                 'saved_index':
    #                     saved_index_path,
    #                 'inspect':
    #                     inspect_path,
    #                 'genes':
    #                     os.path.join(quant_dir, GENES_FILENAME),
    #                 'gene_mtx':
    #                     os.path.join(quant_dir, ABUNDANCE_GENE_FILENAME),
    #                 'gene_tpm_mtx':
    #                     os.path.join(quant_dir, ABUNDANCE_GENE_TPM_FILENAME),
    #                 'mtx':
    #                     os.path.join(quant_dir, ABUNDANCE_FILENAME),
    #                 'tpm_mtx':
    #                     os.path.join(quant_dir, ABUNDANCE_TPM_FILENAME),
    #                 'fld':
    #                     os.path.join(quant_dir, FLD_FILENAME),
    #                 'txnames':
    #                     os.path.join(quant_dir, TXNAMES_FILENAME),
    #                 'barcodes':
    #                     '{}.barcodes.txt'.format(counts_prefix),
    #             }
    #         },
    #                          count.count(
    #                              self.index_path,
    #                              self.t2g_path,
    #                              'SMARTSEQ2',
    #                              out_dir,
    #                              batch_path,
    #                              temp_dir=temp_dir,
    #                              threads=threads,
    #                              memory=memory,
    #                              paired=False,
    #                              h5ad=True,
    #                              tcc=True
    #                          ))
    #         stream_batch.assert_called_once_with(batch_path, temp_dir=temp_dir)
    #         stream_fastqs.assert_not_called()
    #         kallisto_bus.assert_called_once_with(
    #             batch_path,
    #             self.index_path,
    #             'BULK',
    #             out_dir,
    #             threads=threads,
    #             paired=False,
    #             genomebam=False,
    #             strand=None,
    #             gtf_path=None,
    #             chromosomes_path=None,
    #         )
    #         bustools_sort.assert_called_once_with(
    #             bus_path,
    #             bus_s_path,
    #             temp_dir=temp_dir,
    #             threads=threads,
    #             memory=memory,
    #             store_num=False
    #         )
    #         bustools_inspect.assert_called_once_with(
    #             bus_s_path,
    #             inspect_path,
    #             whitelist_path=None,
    #         )
    #         copy_or_create_whitelist.assert_not_called()
    #         bustools_correct.assert_not_called()
    #         bustools_count.assert_called_once_with(
    #             bus_s_path,
    #             counts_prefix,
    #             self.t2g_path,
    #             ecmap_path,
    #             txnames_path,
    #             tcc=True,
    #             mm=True,
    #             cm=True,
    #             umi_gene=True,
    #             em=False,
    #         )
    #         kallisto_quant_tcc.assert_called_once_with(
    #             f'{counts_prefix}.mtx',
    #             saved_index_path,
    #             ecmap_path,
    #             self.t2g_path,
    #             quant_dir,
    #             flens_path=flens_path,
    #             l=None,
    #             s=None,
    #             threads=threads
    #         )
    #         convert_matrix.assert_called_once_with(
    #             quant_dir,
    #             os.path.join(quant_dir, ABUNDANCE_FILENAME),
    #             f'{counts_prefix}.barcodes.txt',
    #             genes_path=os.path.join(quant_dir, TXNAMES_FILENAME),
    #             t2g_path=self.t2g_path,
    #             ec_path=f'{counts_prefix}.ec.txt',
    #             txnames_path=os.path.join(out_dir, TXNAMES_FILENAME),
    #             name='transcript',
    #             loom=False,
    #             h5ad=True,
    #             by_name=False,
    #             tcc=False,
    #             threads=threads
    #         )
    #         filter_with_bustools.assert_not_called()
    # 
    # def test_count_smartseq3(self):
    #     with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
    #         mock.patch('kb_python.count.stream_batch') as stream_batch,\
    #         mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
    #         mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
    #         mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
    #         mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
    #         mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
    #         mock.patch('kb_python.count.write_smartseq3_capture') as write_smartseq3_capture,\
    #         mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
    #         mock.patch('kb_python.count.bustools_count') as bustools_count,\
    #         mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
    #         mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
    #         mock.patch('kb_python.count.STATS') as STATS,\
    #         mock.patch('kb_python.count.render_report'),\
    #         mock.patch('kb_python.count.import_matrix_as_anndata'):
    #         out_dir = self.temp_dir
    #         temp_dir = self.temp_dir
    #         counts_internal_dir = os.path.join(
    #             out_dir, f'{UNFILTERED_COUNTS_DIR}{INTERNAL_SUFFIX}'
    #         )
    #         counts_umi_dir = os.path.join(
    #             out_dir, f'{UNFILTERED_COUNTS_DIR}{UMI_SUFFIX}'
    #         )
    #         counts_internal_prefix = os.path.join(
    #             counts_internal_dir, COUNTS_PREFIX
    #         )
    #         counts_umi_prefix = os.path.join(counts_umi_dir, COUNTS_PREFIX)
    #         threads = 99999
    #         memory = 'TEST'
    #         bus_path = os.path.join(out_dir, BUS_FILENAME)
    #         ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
    #         txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
    #         inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
    #         inspect_internal_path = os.path.join(
    #             out_dir, INSPECT_INTERNAL_FILENAME
    #         )
    #         inspect_umi_path = os.path.join(out_dir, INSPECT_UMI_FILENAME)
    #         info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    #         flens_path = os.path.join(out_dir, FLENS_FILENAME)
    #         saved_index_path = os.path.join(out_dir, SAVED_INDEX_FILENAME)
    #         bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
    #         bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
    #         bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
    #         capture_path = os.path.join(out_dir, CAPTURE_FILENAME)
    #         bus_internal_path = os.path.join(
    #             out_dir, f'output{INTERNAL_SUFFIX}.bus'
    #         )
    #         bus_umi_path = os.path.join(out_dir, f'output{UMI_SUFFIX}.bus')
    #         fastqs = [
    #             self.smartseq3_1_i1_fastq_path, self.smartseq3_1_i2_fastq_path,
    #             self.smartseq3_1_R1_fastq_path, self.smartseq3_1_R2_fastq_path,
    #             self.smartseq3_2_i1_fastq_path, self.smartseq3_2_i2_fastq_path,
    #             self.smartseq3_2_R1_fastq_path, self.smartseq3_2_R2_fastq_path
    #         ]
    #         stream_fastqs.return_value = fastqs
    #         kallisto_bus.return_value = {
    #             'bus': bus_path,
    #             'ecmap': ecmap_path,
    #             'txnames': txnames_path,
    #             'info': info_path,
    #             'flens': flens_path,
    #             'saved_index': saved_index_path
    #         }
    #         bustools_sort.side_effect = [{
    #             'bus': bus_s_path
    #         }, {
    #             'bus': bus_scs_path
    #         }]
    #         bustools_inspect.side_effect = [{
    #             'inspect': inspect_path
    #         }, {
    #             'inspect': inspect_internal_path
    #         }, {
    #             'inspect': inspect_umi_path
    #         }]
    #         copy_or_create_whitelist.return_value = self.whitelist_path
    #         bustools_correct.return_value = {'bus': bus_sc_path}
    #         write_smartseq3_capture.return_value = capture_path
    #         bustools_capture.side_effect = [{
    #             'bus': bus_internal_path
    #         }, {
    #             'bus': bus_umi_path
    #         }]
    #         bustools_count.side_effect = [{
    #             'mtx': f'{counts_internal_prefix}.mtx',
    #             'genes': f'{counts_internal_prefix}.genes.txt',
    #             'barcodes': f'{counts_internal_prefix}.barcodes.txt',
    #         }, {
    #             'mtx': f'{counts_umi_prefix}.mtx',
    #             'genes': f'{counts_umi_prefix}.genes.txt',
    #             'barcodes': f'{counts_umi_prefix}.barcodes.txt',
    #         }]
    #         convert_matrix.side_effect = [{
    #             'h5ad': os.path.join(counts_internal_dir, 'adata.h5ad')
    #         }, {
    #             'h5ad': os.path.join(counts_umi_dir, 'adata.h5ad')
    #         }]
    #         STATS.save.return_value = 'stats'
    # 
    #         self.assertEqual({
    #             'stats': 'stats',
    #             'unfiltered': {
    #                 'bus':
    #                     bus_path,
    #                 'ecmap':
    #                     ecmap_path,
    #                 'txnames':
    #                     txnames_path,
    #                 'info':
    #                     info_path,
    #                 'flens':
    #                     flens_path,
    #                 'saved_index':
    #                     saved_index_path,
    #                 'whitelist':
    #                     self.whitelist_path,
    #                 'inspect':
    #                     inspect_path,
    #                 'inspect_umi':
    #                     inspect_umi_path,
    #                 'inspect_internal':
    #                     inspect_internal_path,
    #                 'bus_scs':
    #                     bus_scs_path,
    #                 'bus_internal':
    #                     bus_internal_path,
    #                 'bus_umi':
    #                     bus_umi_path,
    #                 'mtx_internal':
    #                     f'{counts_internal_prefix}.mtx',
    #                 'genes_internal':
    #                     f'{counts_internal_prefix}.genes.txt',
    #                 'barcodes_internal':
    #                     f'{counts_internal_prefix}.barcodes.txt',
    #                 'mtx_umi':
    #                     f'{counts_umi_prefix}.mtx',
    #                 'genes_umi':
    #                     f'{counts_umi_prefix}.genes.txt',
    #                 'barcodes_umi':
    #                     f'{counts_umi_prefix}.barcodes.txt',
    #                 'h5ad_internal':
    #                     os.path.join(counts_internal_dir, 'adata.h5ad'),
    #                 'h5ad_umi':
    #                     os.path.join(counts_umi_dir, 'adata.h5ad'),
    #             }
    #         },
    #                          count.count(
    #                              self.index_path,
    #                              self.t2g_path,
    #                              "SMARTSEQ3",
    #                              out_dir,
    #                              fastqs,
    #                              temp_dir=temp_dir,
    #                              threads=threads,
    #                              memory=memory,
    #                              h5ad=True
    #                          ))
    #         stream_fastqs.assert_called_once_with(fastqs, temp_dir=temp_dir)
    #         kallisto_bus.assert_called_once_with(
    #             fastqs,
    #             self.index_path,
    #             'SMARTSEQ3',
    #             out_dir,
    #             threads=threads,
    #             paired=True,
    #             genomebam=False,
    #             strand=None,
    #             gtf_path=None,
    #             chromosomes_path=None,
    #         )
    #         self.assertEqual(bustools_sort.call_count, 2)
    #         bustools_sort.assert_has_calls([
    #             call(
    #                 bus_path,
    #                 bus_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 bus_sc_path,
    #                 bus_scs_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             )
    #         ])
    #         bustools_inspect.assert_has_calls([
    #             call(
    #                 bus_s_path,
    #                 inspect_path,
    #                 whitelist_path=self.whitelist_path,
    #             ),
    #             call(
    #                 bus_internal_path,
    #                 inspect_internal_path,
    #                 whitelist_path=self.whitelist_path
    #             ),
    #             call(
    #                 bus_umi_path,
    #                 inspect_umi_path,
    #                 whitelist_path=self.whitelist_path
    #             ),
    #         ])
    #         copy_or_create_whitelist.assert_called_once_with(
    #             'SMARTSEQ3', bus_s_path, out_dir
    #         )
    #         bustools_correct.assert_called_once_with(
    #             bus_s_path, bus_sc_path, self.whitelist_path
    #         )
    #         self.assertEqual(2, bustools_capture.call_count)
    #         bustools_capture.assert_has_calls([
    #             call(
    #                 bus_scs_path,
    #                 bus_internal_path,
    #                 capture_path,
    #                 capture_type='umis',
    #                 complement=False
    #             ),
    #             call(
    #                 bus_scs_path,
    #                 bus_umi_path,
    #                 capture_path,
    #                 capture_type='umis',
    #                 complement=True
    #             )
    #         ])
    #         self.assertEqual(2, bustools_count.call_count)
    #         bustools_count.assert_has_calls([
    #             call(
    #                 bus_internal_path,
    #                 counts_internal_prefix,
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 cm=True,
    #                 umi_gene=False
    #             ),
    #             call(
    #                 bus_umi_path,
    #                 counts_umi_prefix,
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 cm=False,
    #                 umi_gene=True
    #             ),
    #         ])
    #         self.assertEqual(2, convert_matrix.call_count)
    #         convert_matrix.assert_has_calls([
    #             call(
    #                 counts_internal_dir,
    #                 f'{counts_internal_prefix}.mtx',
    #                 f'{counts_internal_prefix}.barcodes.txt',
    #                 genes_path=f'{counts_internal_prefix}.genes.txt',
    #                 t2g_path=self.t2g_path,
    #                 ec_path=None,
    #                 txnames_path=txnames_path,
    #                 name='gene',
    #                 loom=False,
    #                 h5ad=True,
    #                 by_name=False,
    #                 tcc=False,
    #                 threads=threads
    #             ),
    #             call(
    #                 counts_umi_dir,
    #                 f'{counts_umi_prefix}.mtx',
    #                 f'{counts_umi_prefix}.barcodes.txt',
    #                 genes_path=f'{counts_umi_prefix}.genes.txt',
    #                 t2g_path=self.t2g_path,
    #                 ec_path=None,
    #                 txnames_path=txnames_path,
    #                 name='gene',
    #                 loom=False,
    #                 h5ad=True,
    #                 by_name=False,
    #                 tcc=False,
    #                 threads=threads
    #             ),
    #         ])
    #         filter_with_bustools.assert_not_called()
    #         stream_batch.assert_not_called()
    # 
    # def test_count_smartseq3_tcc(self):
    #     with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
    #         mock.patch('kb_python.count.stream_batch') as stream_batch,\
    #         mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
    #         mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
    #         mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
    #         mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
    #         mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
    #         mock.patch('kb_python.count.write_smartseq3_capture') as write_smartseq3_capture,\
    #         mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
    #         mock.patch('kb_python.count.bustools_count') as bustools_count,\
    #         mock.patch('kb_python.count.kallisto_quant_tcc') as kallisto_quant_tcc,\
    #         mock.patch('kb_python.count.convert_matrix') as convert_matrix,\
    #         mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
    #         mock.patch('kb_python.count.STATS') as STATS,\
    #         mock.patch('kb_python.count.render_report'),\
    #         mock.patch('kb_python.count.import_matrix_as_anndata'):
    #         out_dir = self.temp_dir
    #         temp_dir = self.temp_dir
    #         counts_internal_dir = os.path.join(
    #             out_dir, f'{UNFILTERED_COUNTS_DIR}{INTERNAL_SUFFIX}'
    #         )
    #         counts_umi_dir = os.path.join(
    #             out_dir, f'{UNFILTERED_COUNTS_DIR}{UMI_SUFFIX}'
    #         )
    #         counts_internal_prefix = os.path.join(
    #             counts_internal_dir, TCC_PREFIX
    #         )
    #         counts_umi_prefix = os.path.join(counts_umi_dir, TCC_PREFIX)
    #         quant_internal_dir = os.path.join(
    #             out_dir, f'{UNFILTERED_QUANT_DIR}{INTERNAL_SUFFIX}'
    #         )
    #         quant_umi_dir = os.path.join(
    #             out_dir, f'{UNFILTERED_QUANT_DIR}{UMI_SUFFIX}'
    #         )
    #         threads = 99999
    #         memory = 'TEST'
    #         bus_path = os.path.join(out_dir, BUS_FILENAME)
    #         ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
    #         txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
    #         inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
    #         inspect_internal_path = os.path.join(
    #             out_dir, INSPECT_INTERNAL_FILENAME
    #         )
    #         inspect_umi_path = os.path.join(out_dir, INSPECT_UMI_FILENAME)
    #         info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    #         flens_path = os.path.join(out_dir, FLENS_FILENAME)
    #         saved_index_path = os.path.join(out_dir, SAVED_INDEX_FILENAME)
    #         bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
    #         bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
    #         bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
    #         capture_path = os.path.join(out_dir, CAPTURE_FILENAME)
    #         bus_internal_path = os.path.join(
    #             out_dir, f'output{INTERNAL_SUFFIX}.bus'
    #         )
    #         bus_umi_path = os.path.join(out_dir, f'output{UMI_SUFFIX}.bus')
    #         fastqs = [
    #             self.smartseq3_1_i1_fastq_path, self.smartseq3_1_i2_fastq_path,
    #             self.smartseq3_1_R1_fastq_path, self.smartseq3_1_R2_fastq_path,
    #             self.smartseq3_2_i1_fastq_path, self.smartseq3_2_i2_fastq_path,
    #             self.smartseq3_2_R1_fastq_path, self.smartseq3_2_R2_fastq_path
    #         ]
    #         stream_fastqs.return_value = fastqs
    #         kallisto_bus.return_value = {
    #             'bus': bus_path,
    #             'ecmap': ecmap_path,
    #             'txnames': txnames_path,
    #             'info': info_path,
    #             'flens': flens_path,
    #             'saved_index': saved_index_path
    #         }
    #         bustools_sort.side_effect = [{
    #             'bus': bus_s_path
    #         }, {
    #             'bus': bus_scs_path
    #         }]
    #         bustools_inspect.side_effect = [{
    #             'inspect': inspect_path
    #         }, {
    #             'inspect': inspect_internal_path
    #         }, {
    #             'inspect': inspect_umi_path
    #         }]
    #         copy_or_create_whitelist.return_value = self.whitelist_path
    #         bustools_correct.return_value = {'bus': bus_sc_path}
    #         write_smartseq3_capture.return_value = capture_path
    #         bustools_capture.side_effect = [{
    #             'bus': bus_internal_path
    #         }, {
    #             'bus': bus_umi_path
    #         }]
    #         bustools_count.side_effect = [{
    #             'mtx': f'{counts_internal_prefix}.mtx',
    #             'ec': f'{counts_internal_prefix}.ec.txt',
    #             'barcodes': f'{counts_internal_prefix}.barcodes.txt',
    #         }, {
    #             'mtx': f'{counts_umi_prefix}.mtx',
    #             'ec': f'{counts_umi_prefix}.ec.txt',
    #             'barcodes': f'{counts_umi_prefix}.barcodes.txt',
    #         }]
    #         kallisto_quant_tcc.side_effect = [{
    #             'genes':
    #                 os.path.join(quant_internal_dir, GENES_FILENAME),
    #             'gene_mtx':
    #                 os.path.join(quant_internal_dir, ABUNDANCE_GENE_FILENAME),
    #             'gene_tpm_mtx':
    #                 os.path.join(
    #                     quant_internal_dir, ABUNDANCE_GENE_TPM_FILENAME
    #                 ),
    #             'mtx':
    #                 os.path.join(quant_internal_dir, ABUNDANCE_FILENAME),
    #             'tpm_mtx':
    #                 os.path.join(quant_internal_dir, ABUNDANCE_TPM_FILENAME),
    #             'fld':
    #                 os.path.join(quant_internal_dir, FLD_FILENAME),
    #             'txnames':
    #                 os.path.join(quant_internal_dir, TXNAMES_FILENAME),
    #         }, {
    #             'genes':
    #                 os.path.join(quant_umi_dir, GENES_FILENAME),
    #             'gene_mtx':
    #                 os.path.join(quant_umi_dir, ABUNDANCE_GENE_FILENAME),
    #             'gene_tpm_mtx':
    #                 os.path.join(quant_umi_dir, ABUNDANCE_GENE_TPM_FILENAME),
    #             'mtx':
    #                 os.path.join(quant_umi_dir, ABUNDANCE_FILENAME),
    #             'tpm_mtx':
    #                 os.path.join(quant_umi_dir, ABUNDANCE_TPM_FILENAME),
    #             'fld':
    #                 os.path.join(quant_umi_dir, FLD_FILENAME),
    #             'txnames':
    #                 os.path.join(quant_umi_dir, TXNAMES_FILENAME),
    #         }]
    #         convert_matrix.side_effect = [{
    #             'h5ad': os.path.join(counts_internal_dir, 'adata.h5ad')
    #         }, {
    #             'h5ad': os.path.join(counts_umi_dir, 'adata.h5ad')
    #         }]
    #         STATS.save.return_value = 'stats'
    # 
    #         self.assertEqual({
    #             'stats': 'stats',
    #             'unfiltered': {
    #                 'bus':
    #                     bus_path,
    #                 'ecmap':
    #                     ecmap_path,
    #                 'txnames':
    #                     txnames_path,
    #                 'info':
    #                     info_path,
    #                 'flens':
    #                     flens_path,
    #                 'saved_index':
    #                     saved_index_path,
    #                 'whitelist':
    #                     self.whitelist_path,
    #                 'inspect':
    #                     inspect_path,
    #                 'inspect_umi':
    #                     inspect_umi_path,
    #                 'inspect_internal':
    #                     inspect_internal_path,
    #                 'bus_scs':
    #                     bus_scs_path,
    #                 'bus_internal':
    #                     bus_internal_path,
    #                 'bus_umi':
    #                     bus_umi_path,
    #                 'ec_internal':
    #                     f'{counts_internal_prefix}.ec.txt',
    #                 'barcodes_internal':
    #                     f'{counts_internal_prefix}.barcodes.txt',
    #                 'ec_umi':
    #                     f'{counts_umi_prefix}.ec.txt',
    #                 'barcodes_umi':
    #                     f'{counts_umi_prefix}.barcodes.txt',
    #                 'h5ad_internal':
    #                     os.path.join(counts_internal_dir, 'adata.h5ad'),
    #                 'h5ad_umi':
    #                     os.path.join(counts_umi_dir, 'adata.h5ad'),
    #                 'genes_internal':
    #                     os.path.join(quant_internal_dir, GENES_FILENAME),
    #                 'gene_mtx_internal':
    #                     os.path.join(
    #                         quant_internal_dir, ABUNDANCE_GENE_FILENAME
    #                     ),
    #                 'gene_tpm_mtx_internal':
    #                     os.path.join(
    #                         quant_internal_dir, ABUNDANCE_GENE_TPM_FILENAME
    #                     ),
    #                 'mtx_internal':
    #                     os.path.join(quant_internal_dir, ABUNDANCE_FILENAME),
    #                 'tpm_mtx_internal':
    #                     os.path.join(
    #                         quant_internal_dir, ABUNDANCE_TPM_FILENAME
    #                     ),
    #                 'fld_internal':
    #                     os.path.join(quant_internal_dir, FLD_FILENAME),
    #                 'txnames_internal':
    #                     os.path.join(quant_internal_dir, TXNAMES_FILENAME),
    #                 'genes_umi':
    #                     os.path.join(quant_umi_dir, GENES_FILENAME),
    #                 'gene_mtx_umi':
    #                     os.path.join(quant_umi_dir, ABUNDANCE_GENE_FILENAME),
    #                 'gene_tpm_mtx_umi':
    #                     os.path.join(
    #                         quant_umi_dir, ABUNDANCE_GENE_TPM_FILENAME
    #                     ),
    #                 'mtx_umi':
    #                     os.path.join(quant_umi_dir, ABUNDANCE_FILENAME),
    #                 'tpm_mtx_umi':
    #                     os.path.join(quant_umi_dir, ABUNDANCE_TPM_FILENAME),
    #                 'fld_umi':
    #                     os.path.join(quant_umi_dir, FLD_FILENAME),
    #                 'txnames_umi':
    #                     os.path.join(quant_umi_dir, TXNAMES_FILENAME),
    #             }
    #         },
    #                          count.count_smartseq3(
    #                              self.index_path,
    #                              self.t2g_path,
    #                              "SMARTSEQ3",
    #                              out_dir,
    #                              fastqs,
    #                              temp_dir=temp_dir,
    #                              threads=threads,
    #                              memory=memory,
    #                              h5ad=True,
    #                              tcc=True
    #                          ))
    #         stream_fastqs.assert_called_once_with(fastqs, temp_dir=temp_dir)
    #         kallisto_bus.assert_called_once_with(
    #             fastqs,
    #             self.index_path,
    #             'SMARTSEQ3',
    #             out_dir,
    #             threads=threads,
    #             paired=True,
    #             genomebam=False,
    #             strand=None,
    #             gtf_path=None,
    #             chromosomes_path=None,
    #         )
    #         self.assertEqual(bustools_sort.call_count, 2)
    #         bustools_sort.assert_has_calls([
    #             call(
    #                 bus_path,
    #                 bus_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 bus_sc_path,
    #                 bus_scs_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             )
    #         ])
    #         bustools_inspect.assert_has_calls([
    #             call(
    #                 bus_s_path,
    #                 inspect_path,
    #                 whitelist_path=self.whitelist_path,
    #             ),
    #             call(
    #                 bus_internal_path,
    #                 inspect_internal_path,
    #                 whitelist_path=self.whitelist_path
    #             ),
    #             call(
    #                 bus_umi_path,
    #                 inspect_umi_path,
    #                 whitelist_path=self.whitelist_path
    #             ),
    #         ])
    #         copy_or_create_whitelist.assert_called_once_with(
    #             'SMARTSEQ3', bus_s_path, out_dir
    #         )
    #         bustools_correct.assert_called_once_with(
    #             bus_s_path, bus_sc_path, self.whitelist_path
    #         )
    #         self.assertEqual(2, bustools_capture.call_count)
    #         bustools_capture.assert_has_calls([
    #             call(
    #                 bus_scs_path,
    #                 bus_internal_path,
    #                 capture_path,
    #                 capture_type='umis',
    #                 complement=False
    #             ),
    #             call(
    #                 bus_scs_path,
    #                 bus_umi_path,
    #                 capture_path,
    #                 capture_type='umis',
    #                 complement=True
    #             )
    #         ])
    #         self.assertEqual(2, bustools_count.call_count)
    #         bustools_count.assert_has_calls([
    #             call(
    #                 bus_internal_path,
    #                 counts_internal_prefix,
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=True,
    #                 mm=True,
    #                 cm=True,
    #                 umi_gene=False
    #             ),
    #             call(
    #                 bus_umi_path,
    #                 counts_umi_prefix,
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=True,
    #                 mm=True,
    #                 cm=False,
    #                 umi_gene=True
    #             ),
    #         ])
    #         self.assertEqual(2, convert_matrix.call_count)
    #         convert_matrix.assert_has_calls([
    #             call(
    #                 quant_internal_dir,
    #                 os.path.join(quant_internal_dir, ABUNDANCE_FILENAME),
    #                 f'{counts_internal_prefix}.barcodes.txt',
    #                 genes_path=os.path.join(
    #                     quant_internal_dir, TXNAMES_FILENAME
    #                 ),
    #                 t2g_path=self.t2g_path,
    #                 ec_path=f'{counts_internal_prefix}.ec.txt',
    #                 txnames_path=txnames_path,
    #                 name='transcript',
    #                 loom=False,
    #                 h5ad=True,
    #                 by_name=False,
    #                 tcc=False,
    #                 threads=threads
    #             ),
    #             call(
    #                 quant_umi_dir,
    #                 os.path.join(quant_umi_dir, ABUNDANCE_FILENAME),
    #                 f'{counts_umi_prefix}.barcodes.txt',
    #                 genes_path=os.path.join(quant_umi_dir, TXNAMES_FILENAME),
    #                 t2g_path=self.t2g_path,
    #                 ec_path=f'{counts_umi_prefix}.ec.txt',
    #                 txnames_path=txnames_path,
    #                 name='transcript',
    #                 loom=False,
    #                 h5ad=True,
    #                 by_name=False,
    #                 tcc=False,
    #                 threads=threads
    #             ),
    #         ])
    #         filter_with_bustools.assert_not_called()
    #         stream_batch.assert_not_called()

    def test_count_strand(self):
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
            mock.patch('kb_python.count.render_report') as render_report,\
            mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata:
            out_dir = self.temp_dir
            temp_dir = self.temp_dir
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
                                 strand='unstranded'
                             ))

            stream_fastqs.assert_called_once_with(
                self.fastqs, temp_dir=temp_dir
            )
            kallisto_bus.assert_called_once()
            self.assertEqual(bustools_sort.call_count, 2)
            bustools_sort.assert_has_calls([
                call(
                    bus_path,
                    bus_s_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                    store_num=False
                ),
                call(
                    bus_sc_path,
                    bus_scs_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory,
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path,
                inspect_path,
                whitelist_path=self.whitelist_path,
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
                mm=False,
                cm=False,
                umi_gene=True,
                em=False,
                batch_barcodes=False,
            )
            convert_matrix.assert_not_called()
            filter_with_bustools.assert_not_called()

            STATS.start.assert_called_once()
            STATS.end.assert_called_once()
            STATS.save.assert_called_once_with(
                os.path.join(out_dir, KB_INFO_FILENAME)
            )
            import_matrix_as_anndata.assert_not_called()
            render_report.assert_not_called()


    # def test_count_velocity_with_whitelist(self):
    #     with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
    #         mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
    #         mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
    #         mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
    #         mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
    #         mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
    #         mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
    #         mock.patch('kb_python.count.bustools_count') as bustools_count,\
    #         mock.patch('kb_python.count.convert_matrices') as convert_matrices,\
    #         mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
    #         mock.patch('kb_python.count.STATS') as STATS,\
    #         mock.patch('kb_python.count.render_report') as render_report,\
    #         mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata:
    #         out_dir = self.temp_dir
    #         temp_dir = self.temp_dir
    #         counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
    #         threads = 99999
    #         memory = 'TEST'
    #         bus_path = os.path.join(out_dir, BUS_FILENAME)
    #         ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
    #         txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
    #         info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    #         inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
    #         inspect_cdna_path = os.path.join(
    #             out_dir, f'inspect.{BUS_CDNA_PREFIX}.json'
    #         )
    #         inspect_intron_path = os.path.join(
    #             out_dir, f'inspect.{BUS_INTRON_PREFIX}.json'
    #         )
    #         bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
    #         bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
    #         bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
    #         cdna_capture_path = os.path.join(
    #             temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
    #         )
    #         intron_capture_path = os.path.join(
    #             temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
    #         )
    #         cdna_s_path = os.path.join(
    #             out_dir, '{}{}'.format(BUS_CDNA_PREFIX, BUS_UNFILTERED_SUFFIX)
    #         )
    #         intron_s_path = os.path.join(
    #             out_dir,
    #             '{}{}'.format(BUS_INTRON_PREFIX, BUS_UNFILTERED_SUFFIX)
    #         )
    #         cdna_t2c_path = mock.MagicMock()
    #         intron_t2c_path = mock.MagicMock()
    #         stream_fastqs.return_value = self.fastqs
    #         kallisto_bus.return_value = {
    #             'bus': bus_path,
    #             'ecmap': ecmap_path,
    #             'txnames': txnames_path,
    #             'info': info_path
    #         }
    #         bustools_sort.side_effect = [{
    #             'bus': bus_s_path
    #         }, {
    #             'bus': bus_scs_path
    #         }, {
    #             'bus': cdna_s_path
    #         }, {
    #             'bus': intron_s_path
    #         }]
    #         bustools_inspect.side_effect = [{
    #             'inspect': inspect_path
    #         }, {
    #             'inspect': inspect_cdna_path
    #         }, {
    #             'inspect': inspect_intron_path
    #         }]
    #         bustools_capture.side_effect = [{
    #             'bus': cdna_capture_path
    #         }, {
    #             'bus': intron_capture_path
    #         }]
    #         bustools_correct.return_value = {'bus': bus_sc_path}
    #         bustools_count.side_effect = [{
    #             'mtx':
    #                 '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ),
    #         }, {
    #             'mtx':
    #                 '{}.mtx'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #         }]
    #         STATS.save.return_value = 'stats'
    # 
    #         self.assertEqual({
    #             'stats': 'stats',
    #             'unfiltered': {
    #                 'bus': bus_path,
    #                 'ecmap': ecmap_path,
    #                 'txnames': txnames_path,
    #                 'info': info_path,
    #                 'inspect': inspect_path,
    #                 'bus_scs': bus_scs_path,
    #                 BUS_CDNA_PREFIX: {
    #                     'bus':
    #                         cdna_s_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'inspect':
    #                         inspect_cdna_path
    #                 },
    #                 BUS_INTRON_PREFIX: {
    #                     'bus':
    #                         intron_s_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'inspect':
    #                         inspect_intron_path
    #                 }
    #             }
    #         },
    #                          count.count_nac(
    #                              self.index_path,
    #                              self.t2g_path,
    #                              cdna_t2c_path,
    #                              intron_t2c_path,
    #                              self.technology,
    #                              out_dir,
    #                              self.fastqs,
    #                              whitelist_path=self.whitelist_path,
    #                              temp_dir=temp_dir,
    #                              threads=threads,
    #                              memory=memory
    #                          ))
    #         stream_fastqs.assert_called_once_with(
    #             self.fastqs, temp_dir=temp_dir
    #         )
    #         kallisto_bus.assert_called_once_with(
    #             self.fastqs,
    #             self.index_path,
    #             self.technology,
    #             out_dir,
    #             threads=threads,
    #             paired=False,
    #             genomebam=False,
    #             strand=None,
    #             gtf_path=None,
    #             chromosomes_path=None,
    #         )
    #         self.assertEqual(bustools_sort.call_count, 4)
    #         bustools_sort.assert_has_calls([
    #             call(
    #                 bus_path,
    #                 bus_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 bus_sc_path,
    #                 bus_scs_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 cdna_capture_path,
    #                 cdna_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 intron_capture_path,
    #                 intron_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             )
    #         ])
    #         self.assertEqual(3, bustools_inspect.call_count)
    #         bustools_inspect.assert_has_calls([
    #             call(
    #                 bus_s_path,
    #                 inspect_path,
    #                 whitelist_path=self.whitelist_path,
    #             ),
    #             call(
    #                 cdna_s_path,
    #                 inspect_cdna_path,
    #                 whitelist_path=self.whitelist_path,
    #             ),
    #             call(
    #                 intron_s_path,
    #                 inspect_intron_path,
    #                 whitelist_path=self.whitelist_path,
    #             )
    #         ])
    #         copy_or_create_whitelist.assert_not_called()
    #         bustools_correct.assert_called_once_with(
    #             bus_s_path, bus_sc_path, self.whitelist_path
    #         )
    #         self.assertEqual(2, bustools_count.call_count)
    #         bustools_count.assert_has_calls([
    #             call(
    #                 cdna_s_path,
    #                 os.path.join(counts_dir, BUS_CDNA_PREFIX),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 cm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             ),
    #             call(
    #                 intron_s_path,
    #                 os.path.join(counts_dir, BUS_INTRON_PREFIX),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 cm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             )
    #         ])
    #         filter_with_bustools.assert_not_called()
    #         convert_matrices.assert_not_called()
    # 
    #         STATS.start.assert_called_once()
    #         STATS.end.assert_called_once()
    #         STATS.to_dict.assert_not_called()
    #         import_matrix_as_anndata.assert_not_called()
    #         render_report.assert_not_called()
    # 
    # def test_count_velocity_cellranger(self):
    #     with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
    #         mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
    #         mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
    #         mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
    #         mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
    #         mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
    #         mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
    #         mock.patch('kb_python.count.bustools_count') as bustools_count,\
    #         mock.patch('kb_python.count.convert_matrices') as convert_matrices,\
    #         mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
    #         mock.patch('kb_python.count.STATS') as STATS,\
    #         mock.patch('kb_python.count.render_report') as render_report,\
    #         mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata,\
    #         mock.patch('kb_python.count.matrix_to_cellranger') as matrix_to_cellranger:
    #         out_dir = self.temp_dir
    #         temp_dir = self.temp_dir
    #         counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
    #         cellranger_cdna_dir = os.path.join(
    #             counts_dir, f'{CELLRANGER_DIR}_{BUS_CDNA_PREFIX}'
    #         )
    #         cellranger_intron_dir = os.path.join(
    #             counts_dir, f'{CELLRANGER_DIR}_{BUS_INTRON_PREFIX}'
    #         )
    #         threads = 99999
    #         memory = 'TEST'
    #         bus_path = os.path.join(out_dir, BUS_FILENAME)
    #         ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
    #         txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
    #         info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    #         inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
    #         inspect_cdna_path = os.path.join(
    #             out_dir, f'inspect.{BUS_CDNA_PREFIX}.json'
    #         )
    #         inspect_intron_path = os.path.join(
    #             out_dir, f'inspect.{BUS_INTRON_PREFIX}.json'
    #         )
    #         bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
    #         bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
    #         bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
    #         cdna_capture_path = os.path.join(
    #             temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
    #         )
    #         intron_capture_path = os.path.join(
    #             temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
    #         )
    #         cdna_s_path = os.path.join(
    #             out_dir, '{}{}'.format(BUS_CDNA_PREFIX, BUS_UNFILTERED_SUFFIX)
    #         )
    #         intron_s_path = os.path.join(
    #             out_dir,
    #             '{}{}'.format(BUS_INTRON_PREFIX, BUS_UNFILTERED_SUFFIX)
    #         )
    #         cdna_t2c_path = mock.MagicMock()
    #         intron_t2c_path = mock.MagicMock()
    #         stream_fastqs.return_value = self.fastqs
    #         kallisto_bus.return_value = {
    #             'bus': bus_path,
    #             'ecmap': ecmap_path,
    #             'txnames': txnames_path,
    #             'info': info_path
    #         }
    #         bustools_sort.side_effect = [{
    #             'bus': bus_s_path
    #         }, {
    #             'bus': bus_scs_path
    #         }, {
    #             'bus': cdna_s_path
    #         }, {
    #             'bus': intron_s_path
    #         }]
    #         bustools_inspect.side_effect = [{
    #             'inspect': inspect_path
    #         }, {
    #             'inspect': inspect_cdna_path
    #         }, {
    #             'inspect': inspect_intron_path
    #         }]
    #         bustools_capture.side_effect = [{
    #             'bus': cdna_capture_path
    #         }, {
    #             'bus': intron_capture_path
    #         }]
    #         bustools_correct.return_value = {'bus': bus_sc_path}
    #         bustools_count.side_effect = [{
    #             'mtx':
    #                 '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ),
    #         }, {
    #             'mtx':
    #                 '{}.mtx'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #         }]
    #         matrix_to_cellranger.side_effect = [{
    #             'mtx':
    #                 os.path.join(cellranger_cdna_dir, CELLRANGER_MATRIX),
    #             'genes':
    #                 os.path.join(cellranger_cdna_dir, CELLRANGER_GENES),
    #             'barcodes':
    #                 os.path.join(cellranger_cdna_dir, CELLRANGER_BARCODES),
    #         }, {
    #             'mtx':
    #                 os.path.join(cellranger_intron_dir, CELLRANGER_MATRIX),
    #             'genes':
    #                 os.path.join(cellranger_intron_dir, CELLRANGER_GENES),
    #             'barcodes':
    #                 os.path.join(cellranger_intron_dir, CELLRANGER_BARCODES),
    #         }]
    #         STATS.save.return_value = 'stats'
    # 
    #         self.assertEqual({
    #             'stats': 'stats',
    #             'unfiltered': {
    #                 'bus': bus_path,
    #                 'ecmap': ecmap_path,
    #                 'txnames': txnames_path,
    #                 'info': info_path,
    #                 'inspect': inspect_path,
    #                 'bus_scs': bus_scs_path,
    #                 BUS_CDNA_PREFIX: {
    #                     'bus':
    #                         cdna_s_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'inspect':
    #                         inspect_cdna_path,
    #                     'cellranger': {
    #                         'mtx':
    #                             os.path.join(
    #                                 cellranger_cdna_dir, CELLRANGER_MATRIX
    #                             ),
    #                         'genes':
    #                             os.path.join(
    #                                 cellranger_cdna_dir, CELLRANGER_GENES
    #                             ),
    #                         'barcodes':
    #                             os.path.join(
    #                                 cellranger_cdna_dir, CELLRANGER_BARCODES
    #                             ),
    #                     }
    #                 },
    #                 BUS_INTRON_PREFIX: {
    #                     'bus':
    #                         intron_s_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'inspect':
    #                         inspect_intron_path,
    #                     'cellranger': {
    #                         'mtx':
    #                             os.path.join(
    #                                 cellranger_intron_dir, CELLRANGER_MATRIX
    #                             ),
    #                         'genes':
    #                             os.path.join(
    #                                 cellranger_intron_dir, CELLRANGER_GENES
    #                             ),
    #                         'barcodes':
    #                             os.path.join(
    #                                 cellranger_intron_dir, CELLRANGER_BARCODES
    #                             ),
    #                     }
    #                 }
    #             }
    #         },
    #                          count.count_nac(
    #                              self.index_path,
    #                              self.t2g_path,
    #                              cdna_t2c_path,
    #                              intron_t2c_path,
    #                              self.technology,
    #                              out_dir,
    #                              self.fastqs,
    #                              whitelist_path=self.whitelist_path,
    #                              temp_dir=temp_dir,
    #                              threads=threads,
    #                              memory=memory,
    #                              cellranger=True
    #                          ))
    #         stream_fastqs.assert_called_once_with(
    #             self.fastqs, temp_dir=temp_dir
    #         )
    #         kallisto_bus.assert_called_once_with(
    #             self.fastqs,
    #             self.index_path,
    #             self.technology,
    #             out_dir,
    #             threads=threads,
    #             paired=False,
    #             genomebam=False,
    #             strand=None,
    #             gtf_path=None,
    #             chromosomes_path=None,
    #         )
    #         self.assertEqual(bustools_sort.call_count, 4)
    #         bustools_sort.assert_has_calls([
    #             call(
    #                 bus_path,
    #                 bus_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 bus_sc_path,
    #                 bus_scs_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 cdna_capture_path,
    #                 cdna_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 intron_capture_path,
    #                 intron_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             )
    #         ])
    #         self.assertEqual(3, bustools_inspect.call_count)
    #         bustools_inspect.assert_has_calls([
    #             call(
    #                 bus_s_path,
    #                 inspect_path,
    #                 whitelist_path=self.whitelist_path,
    #             ),
    #             call(
    #                 cdna_s_path,
    #                 inspect_cdna_path,
    #                 whitelist_path=self.whitelist_path,
    #             ),
    #             call(
    #                 intron_s_path,
    #                 inspect_intron_path,
    #                 whitelist_path=self.whitelist_path,
    #             )
    #         ])
    #         copy_or_create_whitelist.assert_not_called()
    #         bustools_correct.assert_called_once_with(
    #             bus_s_path, bus_sc_path, self.whitelist_path
    #         )
    #         self.assertEqual(2, bustools_count.call_count)
    #         bustools_count.assert_has_calls([
    #             call(
    #                 cdna_s_path,
    #                 os.path.join(counts_dir, BUS_CDNA_PREFIX),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 cm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             ),
    #             call(
    #                 intron_s_path,
    #                 os.path.join(counts_dir, BUS_INTRON_PREFIX),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 cm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             )
    #         ])
    #         filter_with_bustools.assert_not_called()
    #         convert_matrices.assert_not_called()
    # 
    #         STATS.start.assert_called_once()
    #         STATS.end.assert_called_once()
    #         STATS.to_dict.assert_not_called()
    #         import_matrix_as_anndata.assert_not_called()
    #         render_report.assert_not_called()
    #         self.assertEqual(2, matrix_to_cellranger.call_count)
    #         matrix_to_cellranger.assert_has_calls([
    #             call(
    #                 '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ), '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ), self.t2g_path, cellranger_cdna_dir
    #             ),
    #             call(
    #                 '{}.mtx'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ), '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ), '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ), self.t2g_path, cellranger_intron_dir
    #             ),
    #         ])
    # 
    # def test_count_velocity_report(self):
    #     with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
    #         mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
    #         mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
    #         mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
    #         mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
    #         mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
    #         mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
    #         mock.patch('kb_python.count.bustools_count') as bustools_count,\
    #         mock.patch('kb_python.count.convert_matrices') as convert_matrices,\
    #         mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
    #         mock.patch('kb_python.count.STATS') as STATS,\
    #         mock.patch('kb_python.count.render_report') as render_report,\
    #         mock.patch('kb_python.count.import_matrix_as_anndata'):
    #         out_dir = self.temp_dir
    #         temp_dir = self.temp_dir
    #         counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
    #         threads = 99999
    #         memory = 'TEST'
    #         bus_path = os.path.join(out_dir, BUS_FILENAME)
    #         ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
    #         txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
    #         info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    #         inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
    #         inspect_cdna_path = os.path.join(
    #             out_dir, f'inspect.{BUS_CDNA_PREFIX}.json'
    #         )
    #         inspect_intron_path = os.path.join(
    #             out_dir, f'inspect.{BUS_INTRON_PREFIX}.json'
    #         )
    #         report_path = os.path.join(out_dir, REPORT_NOTEBOOK_FILENAME)
    #         report_cdna_path = os.path.join(
    #             out_dir, f'report.{BUS_CDNA_PREFIX}.html'
    #         )
    #         report_intron_path = os.path.join(
    #             out_dir, f'report.{BUS_INTRON_PREFIX}.html'
    #         )
    #         bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
    #         bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
    #         bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
    #         cdna_capture_path = os.path.join(
    #             temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
    #         )
    #         intron_capture_path = os.path.join(
    #             temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
    #         )
    #         cdna_s_path = os.path.join(
    #             out_dir, '{}{}'.format(BUS_CDNA_PREFIX, BUS_UNFILTERED_SUFFIX)
    #         )
    #         intron_s_path = os.path.join(
    #             out_dir,
    #             '{}{}'.format(BUS_INTRON_PREFIX, BUS_UNFILTERED_SUFFIX)
    #         )
    #         cdna_t2c_path = mock.MagicMock()
    #         intron_t2c_path = mock.MagicMock()
    #         stream_fastqs.return_value = self.fastqs
    #         kallisto_bus.return_value = {
    #             'bus': bus_path,
    #             'ecmap': ecmap_path,
    #             'txnames': txnames_path,
    #             'info': info_path
    #         }
    #         bustools_sort.side_effect = [{
    #             'bus': bus_s_path
    #         }, {
    #             'bus': bus_scs_path
    #         }, {
    #             'bus': cdna_s_path
    #         }, {
    #             'bus': intron_s_path
    #         }]
    #         bustools_inspect.side_effect = [{
    #             'inspect': inspect_path
    #         }, {
    #             'inspect': inspect_cdna_path
    #         }, {
    #             'inspect': inspect_intron_path
    #         }]
    #         bustools_capture.side_effect = [{
    #             'bus': cdna_capture_path
    #         }, {
    #             'bus': intron_capture_path
    #         }]
    #         bustools_correct.return_value = {'bus': bus_sc_path}
    #         bustools_count.side_effect = [{
    #             'mtx':
    #                 '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ),
    #         }, {
    #             'mtx':
    #                 '{}.mtx'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #         }]
    #         render_report.side_effect = [{
    #             'report': report_path
    #         }, {
    #             'report': report_cdna_path
    #         }, {
    #             'report': report_intron_path
    #         }]
    #         STATS.save.return_value = 'stats'
    # 
    #         self.assertEqual({
    #             'stats': 'stats',
    #             'unfiltered': {
    #                 'bus': bus_path,
    #                 'ecmap': ecmap_path,
    #                 'txnames': txnames_path,
    #                 'info': info_path,
    #                 'inspect': inspect_path,
    #                 'report': report_path,
    #                 'bus_scs': bus_scs_path,
    #                 BUS_CDNA_PREFIX: {
    #                     'bus':
    #                         cdna_s_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'report':
    #                         report_cdna_path,
    #                     'inspect':
    #                         inspect_cdna_path
    #                 },
    #                 BUS_INTRON_PREFIX: {
    #                     'bus':
    #                         intron_s_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'report':
    #                         report_intron_path,
    #                     'inspect':
    #                         inspect_intron_path
    #                 }
    #             }
    #         },
    #                          count.count_nac(
    #                              self.index_path,
    #                              self.t2g_path,
    #                              cdna_t2c_path,
    #                              intron_t2c_path,
    #                              self.technology,
    #                              out_dir,
    #                              self.fastqs,
    #                              whitelist_path=self.whitelist_path,
    #                              temp_dir=temp_dir,
    #                              threads=threads,
    #                              memory=memory,
    #                              report=True
    #                          ))
    #         stream_fastqs.assert_called_once_with(
    #             self.fastqs, temp_dir=temp_dir
    #         )
    #         kallisto_bus.assert_called_once_with(
    #             self.fastqs,
    #             self.index_path,
    #             self.technology,
    #             out_dir,
    #             threads=threads,
    #             paired=False,
    #             genomebam=False,
    #             strand=None,
    #             gtf_path=None,
    #             chromosomes_path=None,
    #         )
    #         self.assertEqual(bustools_sort.call_count, 4)
    #         bustools_sort.assert_has_calls([
    #             call(
    #                 bus_path,
    #                 bus_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 bus_sc_path,
    #                 bus_scs_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 cdna_capture_path,
    #                 cdna_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 intron_capture_path,
    #                 intron_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             )
    #         ])
    #         self.assertEqual(3, bustools_inspect.call_count)
    #         bustools_inspect.assert_has_calls([
    #             call(
    #                 bus_s_path,
    #                 inspect_path,
    #                 whitelist_path=self.whitelist_path,
    #             ),
    #             call(
    #                 cdna_s_path,
    #                 inspect_cdna_path,
    #                 whitelist_path=self.whitelist_path,
    #             ),
    #             call(
    #                 intron_s_path,
    #                 inspect_intron_path,
    #                 whitelist_path=self.whitelist_path,
    #             )
    #         ])
    #         copy_or_create_whitelist.assert_not_called()
    #         bustools_correct.assert_called_once_with(
    #             bus_s_path, bus_sc_path, self.whitelist_path
    #         )
    #         self.assertEqual(2, bustools_count.call_count)
    #         bustools_count.assert_has_calls([
    #             call(
    #                 cdna_s_path,
    #                 os.path.join(counts_dir, BUS_CDNA_PREFIX),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 cm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             ),
    #             call(
    #                 intron_s_path,
    #                 os.path.join(counts_dir, BUS_INTRON_PREFIX),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 cm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             )
    #         ])
    #         filter_with_bustools.assert_not_called()
    #         convert_matrices.assert_not_called()
    # 
    #         STATS.start.assert_called_once()
    #         STATS.end.assert_called_once()
    #         self.assertEqual(3, render_report.call_count)
    #         render_report.assert_has_calls([
    #             call(
    #                 'stats',
    #                 info_path,
    #                 inspect_path,
    #                 ANY,
    #                 ANY,
    #                 temp_dir=temp_dir
    #             ),
    #             call(
    #                 'stats',
    #                 info_path,
    #                 inspect_cdna_path,
    #                 ANY,
    #                 ANY,
    #                 '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ),
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ),
    #                 self.t2g_path,
    #                 temp_dir=temp_dir
    #             ),
    #             call(
    #                 'stats',
    #                 info_path,
    #                 inspect_intron_path,
    #                 ANY,
    #                 ANY,
    #                 '{}.mtx'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #                 self.t2g_path,
    #                 temp_dir=temp_dir
    #             )
    #         ])
    # 
    # def test_count_velocity_convert(self):
    #     with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
    #         mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
    #         mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
    #         mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
    #         mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
    #         mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
    #         mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
    #         mock.patch('kb_python.count.bustools_count') as bustools_count,\
    #         mock.patch('kb_python.count.convert_matrices') as convert_matrices,\
    #         mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
    #         mock.patch('kb_python.count.STATS') as STATS,\
    #         mock.patch('kb_python.count.render_report'),\
    #         mock.patch('kb_python.count.import_matrix_as_anndata'):
    #         out_dir = self.temp_dir
    #         temp_dir = self.temp_dir
    #         counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
    #         threads = 99999
    #         memory = 'TEST'
    #         bus_path = os.path.join(out_dir, BUS_FILENAME)
    #         ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
    #         txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
    #         info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    #         inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
    #         inspect_cdna_path = os.path.join(
    #             out_dir, f'inspect.{BUS_CDNA_PREFIX}.json'
    #         )
    #         inspect_intron_path = os.path.join(
    #             out_dir, f'inspect.{BUS_INTRON_PREFIX}.json'
    #         )
    #         bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
    #         bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
    #         bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
    #         cdna_capture_path = os.path.join(
    #             temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
    #         )
    #         intron_capture_path = os.path.join(
    #             temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
    #         )
    #         cdna_s_path = os.path.join(
    #             out_dir, '{}{}'.format(BUS_CDNA_PREFIX, BUS_UNFILTERED_SUFFIX)
    #         )
    #         intron_s_path = os.path.join(
    #             out_dir,
    #             '{}{}'.format(BUS_INTRON_PREFIX, BUS_UNFILTERED_SUFFIX)
    #         )
    #         cdna_t2c_path = mock.MagicMock()
    #         intron_t2c_path = mock.MagicMock()
    #         adata = mock.MagicMock()
    #         loom_path = os.path.join(counts_dir, '{}.loom'.format(ADATA_PREFIX))
    #         adata.write_loom.return_value = loom_path
    #         stream_fastqs.return_value = self.fastqs
    #         kallisto_bus.return_value = {
    #             'bus': bus_path,
    #             'ecmap': ecmap_path,
    #             'txnames': txnames_path,
    #             'info': info_path
    #         }
    #         bustools_sort.side_effect = [{
    #             'bus': bus_s_path
    #         }, {
    #             'bus': bus_scs_path
    #         }, {
    #             'bus': cdna_s_path
    #         }, {
    #             'bus': intron_s_path
    #         }]
    #         bustools_inspect.side_effect = [{
    #             'inspect': inspect_path
    #         }, {
    #             'inspect': inspect_cdna_path
    #         }, {
    #             'inspect': inspect_intron_path
    #         }]
    #         bustools_capture.side_effect = [{
    #             'bus': cdna_capture_path
    #         }, {
    #             'bus': intron_capture_path
    #         }]
    #         bustools_correct.return_value = {'bus': bus_sc_path}
    #         bustools_count.side_effect = [{
    #             'mtx':
    #                 '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ),
    #         }, {
    #             'mtx':
    #                 '{}.mtx'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #         }]
    #         convert_matrices.return_value = {'loom': loom_path}
    #         STATS.save.return_value = 'stats'
    # 
    #         self.assertEqual({
    #             'stats': 'stats',
    #             'unfiltered': {
    #                 'bus': bus_path,
    #                 'ecmap': ecmap_path,
    #                 'txnames': txnames_path,
    #                 'info': info_path,
    #                 'inspect': inspect_path,
    #                 'bus_scs': bus_scs_path,
    #                 'loom': loom_path,
    #                 BUS_CDNA_PREFIX: {
    #                     'bus':
    #                         cdna_s_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'inspect':
    #                         inspect_cdna_path
    #                 },
    #                 BUS_INTRON_PREFIX: {
    #                     'bus':
    #                         intron_s_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'inspect':
    #                         inspect_intron_path
    #                 }
    #             }
    #         },
    #                          count.count_nac(
    #                              self.index_path,
    #                              self.t2g_path,
    #                              cdna_t2c_path,
    #                              intron_t2c_path,
    #                              self.technology,
    #                              out_dir,
    #                              self.fastqs,
    #                              whitelist_path=self.whitelist_path,
    #                              temp_dir=temp_dir,
    #                              threads=threads,
    #                              memory=memory,
    #                              loom=True
    #                          ))
    #         stream_fastqs.assert_called_once_with(
    #             self.fastqs, temp_dir=temp_dir
    #         )
    #         kallisto_bus.assert_called_once_with(
    #             self.fastqs,
    #             self.index_path,
    #             self.technology,
    #             out_dir,
    #             threads=threads,
    #             paired=False,
    #             genomebam=False,
    #             strand=None,
    #             gtf_path=None,
    #             chromosomes_path=None,
    #         )
    #         self.assertEqual(bustools_sort.call_count, 4)
    #         bustools_sort.assert_has_calls([
    #             call(
    #                 bus_path,
    #                 bus_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 bus_sc_path,
    #                 bus_scs_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 cdna_capture_path,
    #                 cdna_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 intron_capture_path,
    #                 intron_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             )
    #         ])
    #         self.assertEqual(3, bustools_inspect.call_count)
    #         bustools_inspect.assert_has_calls([
    #             call(
    #                 bus_s_path,
    #                 inspect_path,
    #                 whitelist_path=self.whitelist_path,
    #             ),
    #             call(
    #                 cdna_s_path,
    #                 inspect_cdna_path,
    #                 whitelist_path=self.whitelist_path,
    #             ),
    #             call(
    #                 intron_s_path,
    #                 inspect_intron_path,
    #                 whitelist_path=self.whitelist_path,
    #             )
    #         ])
    #         copy_or_create_whitelist.assert_not_called()
    #         bustools_correct.assert_called_once_with(
    #             bus_s_path, bus_sc_path, self.whitelist_path
    #         )
    #         self.assertEqual(2, bustools_count.call_count)
    #         bustools_count.assert_has_calls([
    #             call(
    #                 cdna_s_path,
    #                 os.path.join(counts_dir, BUS_CDNA_PREFIX),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 cm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             ),
    #             call(
    #                 intron_s_path,
    #                 os.path.join(counts_dir, BUS_INTRON_PREFIX),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 cm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             )
    #         ])
    #         filter_with_bustools.assert_not_called()
    #         convert_matrices.assert_called_once_with(
    #             counts_dir,
    #             [
    #                 '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
    #                 '{}.mtx'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 )
    #             ],
    #             [
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ), '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 )
    #             ],
    #             genes_paths=[
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ), '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 )
    #             ],
    #             t2g_path=self.t2g_path,
    #             ec_paths=[None, None],
    #             txnames_path=txnames_path,
    #             name='gene',
    #             loom=True,
    #             h5ad=False,
    #             by_name=False,
    #             tcc=False,
    #             nucleus=False,
    #             threads=threads,
    #         )
    # 
    # def test_count_velocity_without_whitelist(self):
    #     with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
    #         mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
    #         mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
    #         mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
    #         mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
    #         mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
    #         mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
    #         mock.patch('kb_python.count.bustools_count') as bustools_count,\
    #         mock.patch('kb_python.count.convert_matrices') as convert_matrices,\
    #         mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
    #         mock.patch('kb_python.count.STATS') as STATS,\
    #         mock.patch('kb_python.count.render_report'),\
    #         mock.patch('kb_python.count.import_matrix_as_anndata'):
    #         out_dir = self.temp_dir
    #         temp_dir = self.temp_dir
    #         counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
    #         threads = 99999
    #         memory = 'TEST'
    #         bus_path = os.path.join(out_dir, BUS_FILENAME)
    #         ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
    #         txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
    #         info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    #         inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
    #         inspect_cdna_path = os.path.join(
    #             out_dir, f'inspect.{BUS_CDNA_PREFIX}.json'
    #         )
    #         inspect_intron_path = os.path.join(
    #             out_dir, f'inspect.{BUS_INTRON_PREFIX}.json'
    #         )
    #         bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
    #         bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
    #         bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
    #         cdna_capture_path = os.path.join(
    #             temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
    #         )
    #         intron_capture_path = os.path.join(
    #             temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
    #         )
    #         cdna_s_path = os.path.join(
    #             out_dir, '{}{}'.format(BUS_CDNA_PREFIX, BUS_UNFILTERED_SUFFIX)
    #         )
    #         intron_s_path = os.path.join(
    #             out_dir,
    #             '{}{}'.format(BUS_INTRON_PREFIX, BUS_UNFILTERED_SUFFIX)
    #         )
    #         cdna_t2c_path = mock.MagicMock()
    #         intron_t2c_path = mock.MagicMock()
    #         stream_fastqs.return_value = self.fastqs
    #         kallisto_bus.return_value = {
    #             'bus': bus_path,
    #             'ecmap': ecmap_path,
    #             'txnames': txnames_path,
    #             'info': info_path
    #         }
    #         bustools_sort.side_effect = [{
    #             'bus': bus_s_path
    #         }, {
    #             'bus': bus_scs_path
    #         }, {
    #             'bus': cdna_s_path
    #         }, {
    #             'bus': intron_s_path
    #         }]
    #         copy_or_create_whitelist.return_value = self.whitelist_path
    #         bustools_inspect.side_effect = [{
    #             'inspect': inspect_path
    #         }, {
    #             'inspect': inspect_cdna_path
    #         }, {
    #             'inspect': inspect_intron_path
    #         }]
    #         bustools_capture.side_effect = [{
    #             'bus': cdna_capture_path
    #         }, {
    #             'bus': intron_capture_path
    #         }]
    #         bustools_correct.return_value = {'bus': bus_sc_path}
    #         bustools_count.side_effect = [{
    #             'mtx':
    #                 '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ),
    #         }, {
    #             'mtx':
    #                 '{}.mtx'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #         }]
    #         STATS.save.return_value = 'stats'
    # 
    #         self.assertEqual({
    #             'stats': 'stats',
    #             'unfiltered': {
    #                 'bus': bus_path,
    #                 'ecmap': ecmap_path,
    #                 'txnames': txnames_path,
    #                 'info': info_path,
    #                 'inspect': inspect_path,
    #                 'bus_scs': bus_scs_path,
    #                 'whitelist': self.whitelist_path,
    #                 BUS_CDNA_PREFIX: {
    #                     'bus':
    #                         cdna_s_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'inspect':
    #                         inspect_cdna_path
    #                 },
    #                 BUS_INTRON_PREFIX: {
    #                     'bus':
    #                         intron_s_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'inspect':
    #                         inspect_intron_path
    #                 }
    #             }
    #         },
    #                          count.count_nac(
    #                              self.index_path,
    #                              self.t2g_path,
    #                              cdna_t2c_path,
    #                              intron_t2c_path,
    #                              self.technology,
    #                              out_dir,
    #                              self.fastqs,
    #                              temp_dir=temp_dir,
    #                              threads=threads,
    #                              memory=memory
    #                          ))
    #         stream_fastqs.assert_called_once_with(
    #             self.fastqs, temp_dir=temp_dir
    #         )
    #         kallisto_bus.assert_called_once_with(
    #             self.fastqs,
    #             self.index_path,
    #             self.technology,
    #             out_dir,
    #             threads=threads,
    #             paired=False,
    #             genomebam=False,
    #             strand=None,
    #             gtf_path=None,
    #             chromosomes_path=None,
    #         )
    #         self.assertEqual(bustools_sort.call_count, 4)
    #         bustools_sort.assert_has_calls([
    #             call(
    #                 bus_path,
    #                 bus_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 bus_sc_path,
    #                 bus_scs_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 cdna_capture_path,
    #                 cdna_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 intron_capture_path,
    #                 intron_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             )
    #         ])
    #         self.assertEqual(3, bustools_inspect.call_count)
    #         bustools_inspect.assert_has_calls([
    #             call(
    #                 bus_s_path,
    #                 inspect_path,
    #                 whitelist_path=self.whitelist_path,
    #             ),
    #             call(
    #                 cdna_s_path,
    #                 inspect_cdna_path,
    #                 whitelist_path=self.whitelist_path,
    #             ),
    #             call(
    #                 intron_s_path,
    #                 inspect_intron_path,
    #                 whitelist_path=self.whitelist_path,
    #             )
    #         ])
    #         copy_or_create_whitelist.assert_called_once_with(
    #             self.technology, bus_s_path, out_dir
    #         )
    #         bustools_correct.assert_called_once_with(
    #             bus_s_path, bus_sc_path, self.whitelist_path
    #         )
    #         self.assertEqual(2, bustools_count.call_count)
    #         bustools_count.assert_has_calls([
    #             call(
    #                 cdna_s_path,
    #                 os.path.join(counts_dir, BUS_CDNA_PREFIX),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 cm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             ),
    #             call(
    #                 intron_s_path,
    #                 os.path.join(counts_dir, BUS_INTRON_PREFIX),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 cm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             )
    #         ])
    #         filter_with_bustools.assert_not_called()
    #         convert_matrices.assert_not_called()
    # 
    # def test_count_velocity_filter(self):
    #     with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
    #         mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
    #         mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
    #         mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
    #         mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
    #         mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
    #         mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
    #         mock.patch('kb_python.count.bustools_count') as bustools_count,\
    #         mock.patch('kb_python.count.convert_matrices') as convert_matrices,\
    #         mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
    #         mock.patch('kb_python.count.STATS') as STATS,\
    #         mock.patch('kb_python.count.render_report'),\
    #         mock.patch('kb_python.count.import_matrix_as_anndata'):
    #         out_dir = self.temp_dir
    #         temp_dir = self.temp_dir
    #         counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
    #         threads = 99999
    #         memory = 'TEST'
    #         bus_path = os.path.join(out_dir, BUS_FILENAME)
    #         ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
    #         txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
    #         info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    #         inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
    #         inspect_cdna_path = os.path.join(
    #             out_dir, f'inspect.{BUS_CDNA_PREFIX}.json'
    #         )
    #         inspect_intron_path = os.path.join(
    #             out_dir, f'inspect.{BUS_INTRON_PREFIX}.json'
    #         )
    #         bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
    #         bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
    #         bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
    #         cdna_capture_path = os.path.join(
    #             temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
    #         )
    #         intron_capture_path = os.path.join(
    #             temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
    #         )
    #         cdna_s_path = os.path.join(
    #             out_dir, '{}{}'.format(BUS_CDNA_PREFIX, BUS_UNFILTERED_SUFFIX)
    #         )
    #         intron_s_path = os.path.join(
    #             out_dir,
    #             '{}{}'.format(BUS_INTRON_PREFIX, BUS_UNFILTERED_SUFFIX)
    #         )
    #         cdna_t2c_path = mock.MagicMock()
    #         intron_t2c_path = mock.MagicMock()
    #         stream_fastqs.return_value = self.fastqs
    #         kallisto_bus.return_value = {
    #             'bus': bus_path,
    #             'ecmap': ecmap_path,
    #             'txnames': txnames_path,
    #             'info': info_path
    #         }
    #         cdna_filtered_path = mock.MagicMock()
    #         intron_filtered_path = mock.MagicMock()
    #         bustools_sort.side_effect = [{
    #             'bus': bus_s_path
    #         }, {
    #             'bus': bus_scs_path
    #         }, {
    #             'bus': cdna_s_path
    #         }, {
    #             'bus': intron_s_path
    #         }, {
    #             'bus': cdna_filtered_path
    #         }, {
    #             'bus': intron_filtered_path
    #         }]
    #         bustools_inspect.side_effect = [{
    #             'inspect': inspect_path
    #         }, {
    #             'inspect': inspect_cdna_path
    #         }, {
    #             'inspect': inspect_intron_path
    #         }]
    #         cdna_filtered_capture_path = mock.MagicMock()
    #         intron_filtered_capture_path = mock.MagicMock()
    #         bustools_capture.side_effect = [{
    #             'bus': cdna_capture_path
    #         }, {
    #             'bus': intron_capture_path
    #         }, {
    #             'bus': cdna_filtered_capture_path
    #         }, {
    #             'bus': intron_filtered_capture_path
    #         }]
    #         bustools_correct.return_value = {'bus': bus_sc_path}
    #         bustools_count.side_effect = [{
    #             'mtx':
    #                 '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ),
    #         }, {
    #             'mtx':
    #                 '{}.mtx'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #         }, {
    #             'mtx':
    #                 '{}.mtx'.format(
    #                     os.path.join(
    #                         out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX
    #                     )
    #                 ),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(
    #                         out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX
    #                     )
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(
    #                         out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX
    #                     )
    #                 ),
    #         }, {
    #             'mtx':
    #                 '{}.mtx'.format(
    #                     os.path.join(
    #                         out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
    #                     )
    #                 ),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(
    #                         out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
    #                     )
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(
    #                         out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
    #                     )
    #                 ),
    #         }]
    #         filtered_whitelist_path = os.path.join(
    #             out_dir, FILTER_WHITELIST_FILENAME
    #         )
    #         filtered_bus_path = os.path.join(out_dir, BUS_FILTERED_FILENAME)
    # 
    #         filter_result = {
    #             'whitelist': filtered_whitelist_path,
    #             'bus_scs': filtered_bus_path,
    #         }
    #         filter_with_bustools.return_value = filter_result
    #         STATS.save.return_value = 'stats'
    # 
    #         self.assertEqual({
    #             'stats': 'stats',
    #             'filtered': {
    #                 'whitelist': filtered_whitelist_path,
    #                 'bus_scs': filtered_bus_path,
    #                 BUS_CDNA_PREFIX: {
    #                     'bus':
    #                         cdna_filtered_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(
    #                                 out_dir, FILTERED_COUNTS_DIR,
    #                                 BUS_CDNA_PREFIX
    #                             )
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(
    #                                 out_dir, FILTERED_COUNTS_DIR,
    #                                 BUS_CDNA_PREFIX
    #                             )
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(
    #                                 out_dir, FILTERED_COUNTS_DIR,
    #                                 BUS_CDNA_PREFIX
    #                             )
    #                         ),
    #                 },
    #                 BUS_INTRON_PREFIX: {
    #                     'bus':
    #                         intron_filtered_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(
    #                                 out_dir, FILTERED_COUNTS_DIR,
    #                                 BUS_INTRON_PREFIX
    #                             )
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(
    #                                 out_dir, FILTERED_COUNTS_DIR,
    #                                 BUS_INTRON_PREFIX
    #                             )
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(
    #                                 out_dir, FILTERED_COUNTS_DIR,
    #                                 BUS_INTRON_PREFIX
    #                             )
    #                         ),
    #                 }
    #             },
    #             'unfiltered': {
    #                 'bus': bus_path,
    #                 'ecmap': ecmap_path,
    #                 'txnames': txnames_path,
    #                 'info': info_path,
    #                 'inspect': inspect_path,
    #                 'bus_scs': bus_scs_path,
    #                 BUS_CDNA_PREFIX: {
    #                     'bus':
    #                         cdna_s_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'inspect':
    #                         inspect_cdna_path
    #                 },
    #                 BUS_INTRON_PREFIX: {
    #                     'bus':
    #                         intron_s_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'inspect':
    #                         inspect_intron_path
    #                 }
    #             }
    #         },
    #                          count.count_nac(
    #                              self.index_path,
    #                              self.t2g_path,
    #                              cdna_t2c_path,
    #                              intron_t2c_path,
    #                              self.technology,
    #                              out_dir,
    #                              self.fastqs,
    #                              filter='bustools',
    #                              whitelist_path=self.whitelist_path,
    #                              temp_dir=temp_dir,
    #                              threads=threads,
    #                              memory=memory
    #                          ))
    #         stream_fastqs.assert_called_once_with(
    #             self.fastqs, temp_dir=temp_dir
    #         )
    #         kallisto_bus.assert_called_once_with(
    #             self.fastqs,
    #             self.index_path,
    #             self.technology,
    #             out_dir,
    #             threads=threads,
    #             paired=False,
    #             genomebam=False,
    #             strand=None,
    #             gtf_path=None,
    #             chromosomes_path=None,
    #         )
    #         self.assertEqual(bustools_sort.call_count, 6)
    #         bustools_sort.assert_has_calls([
    #             call(
    #                 bus_path,
    #                 bus_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 bus_sc_path,
    #                 bus_scs_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 cdna_capture_path,
    #                 cdna_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 intron_capture_path,
    #                 intron_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 cdna_filtered_capture_path,
    #                 os.path.join(
    #                     out_dir,
    #                     '{}{}'.format(BUS_CDNA_PREFIX, BUS_FILTERED_SUFFIX)
    #                 ),
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 intron_filtered_capture_path,
    #                 os.path.join(
    #                     out_dir,
    #                     '{}{}'.format(BUS_INTRON_PREFIX, BUS_FILTERED_SUFFIX)
    #                 ),
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             )
    #         ])
    #         self.assertEqual(3, bustools_inspect.call_count)
    #         bustools_inspect.assert_has_calls([
    #             call(
    #                 bus_s_path,
    #                 inspect_path,
    #                 whitelist_path=self.whitelist_path,
    #             ),
    #             call(
    #                 cdna_s_path,
    #                 inspect_cdna_path,
    #                 whitelist_path=self.whitelist_path,
    #             ),
    #             call(
    #                 intron_s_path,
    #                 inspect_intron_path,
    #                 whitelist_path=self.whitelist_path,
    #             )
    #         ])
    #         copy_or_create_whitelist.assert_not_called()
    #         bustools_correct.assert_called_once_with(
    #             bus_s_path, bus_sc_path, self.whitelist_path
    #         )
    #         self.assertEqual(4, bustools_count.call_count)
    #         bustools_count.assert_has_calls([
    #             call(
    #                 cdna_s_path,
    #                 os.path.join(counts_dir, BUS_CDNA_PREFIX),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 cm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             ),
    #             call(
    #                 intron_s_path,
    #                 os.path.join(counts_dir, BUS_INTRON_PREFIX),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 cm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             ),
    #             call(
    #                 cdna_filtered_path,
    #                 os.path.join(out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             ),
    #             call(
    #                 intron_filtered_path,
    #                 os.path.join(
    #                     out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
    #                 ),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             )
    #         ])
    #         filter_with_bustools.assert_called_once_with(
    #             bus_scs_path,
    #             ecmap_path,
    #             txnames_path,
    #             self.t2g_path,
    #             filtered_whitelist_path,
    #             filtered_bus_path,
    #             filter_threshold=None,
    #             temp_dir=temp_dir,
    #             memory=memory,
    #             count=False,
    #             umi_gene=False,
    #             em=False,
    #         )
    #         convert_matrices.assert_not_called()
    # 
    # def test_count_velocity_filter_convert(self):
    #     with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
    #         mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
    #         mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
    #         mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
    #         mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
    #         mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
    #         mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
    #         mock.patch('kb_python.count.bustools_count') as bustools_count,\
    #         mock.patch('kb_python.count.convert_matrices') as convert_matrices,\
    #         mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
    #         mock.patch('kb_python.count.STATS') as STATS,\
    #         mock.patch('kb_python.count.render_report'),\
    #         mock.patch('kb_python.count.import_matrix_as_anndata'):
    #         out_dir = self.temp_dir
    #         temp_dir = self.temp_dir
    #         counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
    #         threads = 99999
    #         memory = 'TEST'
    #         bus_path = os.path.join(out_dir, BUS_FILENAME)
    #         ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
    #         txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
    #         info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    #         inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
    #         inspect_cdna_path = os.path.join(
    #             out_dir, f'inspect.{BUS_CDNA_PREFIX}.json'
    #         )
    #         inspect_intron_path = os.path.join(
    #             out_dir, f'inspect.{BUS_INTRON_PREFIX}.json'
    #         )
    #         bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
    #         bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
    #         bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
    #         cdna_capture_path = os.path.join(
    #             temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
    #         )
    #         intron_capture_path = os.path.join(
    #             temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
    #         )
    #         cdna_s_path = os.path.join(
    #             out_dir, '{}{}'.format(BUS_CDNA_PREFIX, BUS_UNFILTERED_SUFFIX)
    #         )
    #         intron_s_path = os.path.join(
    #             out_dir,
    #             '{}{}'.format(BUS_INTRON_PREFIX, BUS_UNFILTERED_SUFFIX)
    #         )
    #         cdna_t2c_path = mock.MagicMock()
    #         intron_t2c_path = mock.MagicMock()
    #         stream_fastqs.return_value = self.fastqs
    #         kallisto_bus.return_value = {
    #             'bus': bus_path,
    #             'ecmap': ecmap_path,
    #             'txnames': txnames_path,
    #             'info': info_path
    #         }
    #         cdna_filtered_path = mock.MagicMock()
    #         intron_filtered_path = mock.MagicMock()
    #         bustools_sort.side_effect = [{
    #             'bus': bus_s_path
    #         }, {
    #             'bus': bus_scs_path
    #         }, {
    #             'bus': cdna_s_path
    #         }, {
    #             'bus': intron_s_path
    #         }, {
    #             'bus': cdna_filtered_path
    #         }, {
    #             'bus': intron_filtered_path
    #         }]
    #         bustools_inspect.side_effect = [{
    #             'inspect': inspect_path
    #         }, {
    #             'inspect': inspect_cdna_path
    #         }, {
    #             'inspect': inspect_intron_path
    #         }]
    #         cdna_filtered_capture_path = mock.MagicMock()
    #         intron_filtered_capture_path = mock.MagicMock()
    #         bustools_capture.side_effect = [{
    #             'bus': cdna_capture_path
    #         }, {
    #             'bus': intron_capture_path
    #         }, {
    #             'bus': cdna_filtered_capture_path
    #         }, {
    #             'bus': intron_filtered_capture_path
    #         }]
    #         bustools_correct.return_value = {'bus': bus_sc_path}
    #         bustools_count.side_effect = [{
    #             'mtx':
    #                 '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ),
    #         }, {
    #             'mtx':
    #                 '{}.mtx'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #         }, {
    #             'mtx':
    #                 '{}.mtx'.format(
    #                     os.path.join(
    #                         out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX
    #                     )
    #                 ),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(
    #                         out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX
    #                     )
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(
    #                         out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX
    #                     )
    #                 ),
    #         }, {
    #             'mtx':
    #                 '{}.mtx'.format(
    #                     os.path.join(
    #                         out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
    #                     )
    #                 ),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(
    #                         out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
    #                     )
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(
    #                         out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
    #                     )
    #                 ),
    #         }]
    #         filtered_whitelist_path = os.path.join(
    #             out_dir, FILTER_WHITELIST_FILENAME
    #         )
    #         filtered_bus_path = os.path.join(out_dir, BUS_FILTERED_FILENAME)
    # 
    #         filter_result = {
    #             'whitelist': filtered_whitelist_path,
    #             'bus_scs': filtered_bus_path,
    #         }
    #         filter_with_bustools.return_value = filter_result
    #         STATS.save.return_value = 'stats'
    # 
    #         self.assertEqual({
    #             'stats': 'stats',
    #             'filtered': {
    #                 'whitelist': filtered_whitelist_path,
    #                 'bus_scs': filtered_bus_path,
    #                 BUS_CDNA_PREFIX: {
    #                     'bus':
    #                         cdna_filtered_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(
    #                                 out_dir, FILTERED_COUNTS_DIR,
    #                                 BUS_CDNA_PREFIX
    #                             )
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(
    #                                 out_dir, FILTERED_COUNTS_DIR,
    #                                 BUS_CDNA_PREFIX
    #                             )
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(
    #                                 out_dir, FILTERED_COUNTS_DIR,
    #                                 BUS_CDNA_PREFIX
    #                             )
    #                         ),
    #                 },
    #                 BUS_INTRON_PREFIX: {
    #                     'bus':
    #                         intron_filtered_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(
    #                                 out_dir, FILTERED_COUNTS_DIR,
    #                                 BUS_INTRON_PREFIX
    #                             )
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(
    #                                 out_dir, FILTERED_COUNTS_DIR,
    #                                 BUS_INTRON_PREFIX
    #                             )
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(
    #                                 out_dir, FILTERED_COUNTS_DIR,
    #                                 BUS_INTRON_PREFIX
    #                             )
    #                         ),
    #                 }
    #             },
    #             'unfiltered': {
    #                 'bus': bus_path,
    #                 'ecmap': ecmap_path,
    #                 'txnames': txnames_path,
    #                 'info': info_path,
    #                 'inspect': inspect_path,
    #                 'bus_scs': bus_scs_path,
    #                 BUS_CDNA_PREFIX: {
    #                     'bus':
    #                         cdna_s_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'inspect':
    #                         inspect_cdna_path
    #                 },
    #                 BUS_INTRON_PREFIX: {
    #                     'bus':
    #                         intron_s_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'inspect':
    #                         inspect_intron_path
    #                 }
    #             }
    #         },
    #                          count.count_nac(
    #                              self.index_path,
    #                              self.t2g_path,
    #                              cdna_t2c_path,
    #                              intron_t2c_path,
    #                              self.technology,
    #                              out_dir,
    #                              self.fastqs,
    #                              filter='bustools',
    #                              whitelist_path=self.whitelist_path,
    #                              temp_dir=temp_dir,
    #                              threads=threads,
    #                              memory=memory,
    #                              loom=True,
    #                          ))
    #         stream_fastqs.assert_called_once_with(
    #             self.fastqs, temp_dir=temp_dir
    #         )
    #         kallisto_bus.assert_called_once_with(
    #             self.fastqs,
    #             self.index_path,
    #             self.technology,
    #             out_dir,
    #             threads=threads,
    #             paired=False,
    #             genomebam=False,
    #             strand=None,
    #             gtf_path=None,
    #             chromosomes_path=None,
    #         )
    #         self.assertEqual(bustools_sort.call_count, 6)
    #         bustools_sort.assert_has_calls([
    #             call(
    #                 bus_path,
    #                 bus_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 bus_sc_path,
    #                 bus_scs_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 cdna_capture_path,
    #                 cdna_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 intron_capture_path,
    #                 intron_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 cdna_filtered_capture_path,
    #                 os.path.join(
    #                     out_dir,
    #                     '{}{}'.format(BUS_CDNA_PREFIX, BUS_FILTERED_SUFFIX)
    #                 ),
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 intron_filtered_capture_path,
    #                 os.path.join(
    #                     out_dir,
    #                     '{}{}'.format(BUS_INTRON_PREFIX, BUS_FILTERED_SUFFIX)
    #                 ),
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             )
    #         ])
    #         self.assertEqual(3, bustools_inspect.call_count)
    #         bustools_inspect.assert_has_calls([
    #             call(
    #                 bus_s_path,
    #                 inspect_path,
    #                 whitelist_path=self.whitelist_path,
    #             ),
    #             call(
    #                 cdna_s_path,
    #                 inspect_cdna_path,
    #                 whitelist_path=self.whitelist_path,
    #             ),
    #             call(
    #                 intron_s_path,
    #                 inspect_intron_path,
    #                 whitelist_path=self.whitelist_path,
    #             )
    #         ])
    #         copy_or_create_whitelist.assert_not_called()
    #         bustools_correct.assert_called_once_with(
    #             bus_s_path, bus_sc_path, self.whitelist_path
    #         )
    #         self.assertEqual(4, bustools_count.call_count)
    #         bustools_count.assert_has_calls([
    #             call(
    #                 cdna_s_path,
    #                 os.path.join(counts_dir, BUS_CDNA_PREFIX),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 cm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             ),
    #             call(
    #                 intron_s_path,
    #                 os.path.join(counts_dir, BUS_INTRON_PREFIX),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 cm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             ),
    #             call(
    #                 cdna_filtered_path,
    #                 os.path.join(out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             ),
    #             call(
    #                 intron_filtered_path,
    #                 os.path.join(
    #                     out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
    #                 ),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             )
    #         ])
    #         filter_with_bustools.assert_called_once_with(
    #             bus_scs_path,
    #             ecmap_path,
    #             txnames_path,
    #             self.t2g_path,
    #             filtered_whitelist_path,
    #             filtered_bus_path,
    #             filter_threshold=None,
    #             temp_dir=temp_dir,
    #             memory=memory,
    #             count=False,
    #             umi_gene=False,
    #             em=False,
    #         )
    #         self.assertEqual(2, convert_matrices.call_count)
    #         args = [
    #             call(
    #                 counts_dir,
    #                 [
    #                     '{}.mtx'.format(
    #                         os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                     ), '{}.mtx'.format(
    #                         os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                     )
    #                 ],
    #                 [
    #                     '{}.barcodes.txt'.format(
    #                         os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                     ), '{}.barcodes.txt'.format(
    #                         os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                     )
    #                 ],
    #                 genes_paths=[
    #                     '{}.genes.txt'.format(
    #                         os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                     ), '{}.genes.txt'.format(
    #                         os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                     )
    #                 ],
    #                 ec_paths=[None, None],
    #                 t2g_path=self.t2g_path,
    #                 txnames_path=txnames_path,
    #                 loom=True,
    #                 h5ad=False,
    #                 name='gene',
    #                 by_name=False,
    #                 tcc=False,
    #                 nucleus=False,
    #                 threads=threads,
    #             ),
    #             call(
    #                 os.path.join(out_dir, FILTERED_COUNTS_DIR),
    #                 [
    #                     '{}.mtx'.format(
    #                         os.path.join(
    #                             out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX
    #                         )
    #                     ), '{}.mtx'.format(
    #                         os.path.join(
    #                             out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
    #                         )
    #                     )
    #                 ],
    #                 [
    #                     '{}.barcodes.txt'.format(
    #                         os.path.join(
    #                             out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX
    #                         )
    #                     ), '{}.barcodes.txt'.format(
    #                         os.path.join(
    #                             out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
    #                         )
    #                     )
    #                 ],
    #                 genes_paths=[
    #                     '{}.genes.txt'.format(
    #                         os.path.join(
    #                             out_dir, FILTERED_COUNTS_DIR, BUS_CDNA_PREFIX
    #                         )
    #                     ), '{}.genes.txt'.format(
    #                         os.path.join(
    #                             out_dir, FILTERED_COUNTS_DIR, BUS_INTRON_PREFIX
    #                         )
    #                     )
    #                 ],
    #                 ec_paths=[None, None],
    #                 t2g_path=self.t2g_path,
    #                 txnames_path=txnames_path,
    #                 loom=True,
    #                 h5ad=False,
    #                 by_name=False,
    #                 tcc=False,
    #                 nucleus=False,
    #                 threads=threads,
    #             )
    #         ]
    #         self.assertEqual(args[0], convert_matrices.call_args_list[0])
    #         self.assertEqual(args[1], convert_matrices.call_args_list[1])
    # 
    # def test_count_velocity_strand(self):
    #     with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
    #         mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
    #         mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
    #         mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
    #         mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
    #         mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
    #         mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
    #         mock.patch('kb_python.count.bustools_count') as bustools_count,\
    #         mock.patch('kb_python.count.convert_matrices') as convert_matrices,\
    #         mock.patch('kb_python.count.filter_with_bustools') as filter_with_bustools,\
    #         mock.patch('kb_python.count.STATS') as STATS,\
    #         mock.patch('kb_python.count.render_report') as render_report,\
    #         mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata:
    #         out_dir = self.temp_dir
    #         temp_dir = self.temp_dir
    #         counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
    #         threads = 99999
    #         memory = 'TEST'
    #         bus_path = os.path.join(out_dir, BUS_FILENAME)
    #         ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
    #         txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
    #         info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    #         inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
    #         inspect_cdna_path = os.path.join(
    #             out_dir, f'inspect.{BUS_CDNA_PREFIX}.json'
    #         )
    #         inspect_intron_path = os.path.join(
    #             out_dir, f'inspect.{BUS_INTRON_PREFIX}.json'
    #         )
    #         bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
    #         bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
    #         bus_scs_path = os.path.join(out_dir, BUS_UNFILTERED_FILENAME)
    #         cdna_capture_path = os.path.join(
    #             temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
    #         )
    #         intron_capture_path = os.path.join(
    #             temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
    #         )
    #         cdna_s_path = os.path.join(
    #             out_dir, '{}{}'.format(BUS_CDNA_PREFIX, BUS_UNFILTERED_SUFFIX)
    #         )
    #         intron_s_path = os.path.join(
    #             out_dir,
    #             '{}{}'.format(BUS_INTRON_PREFIX, BUS_UNFILTERED_SUFFIX)
    #         )
    #         cdna_t2c_path = mock.MagicMock()
    #         intron_t2c_path = mock.MagicMock()
    #         stream_fastqs.return_value = self.fastqs
    #         kallisto_bus.return_value = {
    #             'bus': bus_path,
    #             'ecmap': ecmap_path,
    #             'txnames': txnames_path,
    #             'info': info_path
    #         }
    #         bustools_sort.side_effect = [{
    #             'bus': bus_s_path
    #         }, {
    #             'bus': bus_scs_path
    #         }, {
    #             'bus': cdna_s_path
    #         }, {
    #             'bus': intron_s_path
    #         }]
    #         bustools_inspect.side_effect = [{
    #             'inspect': inspect_path
    #         }, {
    #             'inspect': inspect_cdna_path
    #         }, {
    #             'inspect': inspect_intron_path
    #         }]
    #         bustools_capture.side_effect = [{
    #             'bus': cdna_capture_path
    #         }, {
    #             'bus': intron_capture_path
    #         }]
    #         bustools_correct.return_value = {'bus': bus_sc_path}
    #         bustools_count.side_effect = [{
    #             'mtx':
    #                 '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                 ),
    #         }, {
    #             'mtx':
    #                 '{}.mtx'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #             'genes':
    #                 '{}.genes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #             'barcodes':
    #                 '{}.barcodes.txt'.format(
    #                     os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                 ),
    #         }]
    #         STATS.save.return_value = 'stats'
    # 
    #         self.assertEqual({
    #             'stats': 'stats',
    #             'unfiltered': {
    #                 'bus': bus_path,
    #                 'ecmap': ecmap_path,
    #                 'txnames': txnames_path,
    #                 'info': info_path,
    #                 'inspect': inspect_path,
    #                 'bus_scs': bus_scs_path,
    #                 BUS_CDNA_PREFIX: {
    #                     'bus':
    #                         cdna_s_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(counts_dir, BUS_CDNA_PREFIX)
    #                         ),
    #                     'inspect':
    #                         inspect_cdna_path
    #                 },
    #                 BUS_INTRON_PREFIX: {
    #                     'bus':
    #                         intron_s_path,
    #                     'mtx':
    #                         '{}.mtx'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'genes':
    #                         '{}.genes.txt'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'barcodes':
    #                         '{}.barcodes.txt'.format(
    #                             os.path.join(counts_dir, BUS_INTRON_PREFIX)
    #                         ),
    #                     'inspect':
    #                         inspect_intron_path
    #                 }
    #             }
    #         },
    #                          count.count_nac(
    #                              self.index_path,
    #                              self.t2g_path,
    #                              cdna_t2c_path,
    #                              intron_t2c_path,
    #                              self.technology,
    #                              out_dir,
    #                              self.fastqs,
    #                              whitelist_path=self.whitelist_path,
    #                              temp_dir=temp_dir,
    #                              threads=threads,
    #                              memory=memory,
    #                              strand='unstranded'
    #                          ))
    #         stream_fastqs.assert_called_once_with(
    #             self.fastqs, temp_dir=temp_dir
    #         )
    #         kallisto_bus.assert_called_once_with(
    #             self.fastqs,
    #             self.index_path,
    #             self.technology,
    #             out_dir,
    #             threads=threads,
    #             paired=False,
    #             genomebam=False,
    #             strand='unstranded',
    #             gtf_path=None,
    #             chromosomes_path=None,
    #         )
    #         self.assertEqual(bustools_sort.call_count, 4)
    #         bustools_sort.assert_has_calls([
    #             call(
    #                 bus_path,
    #                 bus_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 bus_sc_path,
    #                 bus_scs_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 cdna_capture_path,
    #                 cdna_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             ),
    #             call(
    #                 intron_capture_path,
    #                 intron_s_path,
    #                 temp_dir=temp_dir,
    #                 threads=threads,
    #                 memory=memory,
    #                 store_num=False
    #             )
    #         ])
    #         self.assertEqual(3, bustools_inspect.call_count)
    #         bustools_inspect.assert_has_calls([
    #             call(
    #                 bus_s_path,
    #                 inspect_path,
    #                 whitelist_path=self.whitelist_path,
    #             ),
    #             call(
    #                 cdna_s_path,
    #                 inspect_cdna_path,
    #                 whitelist_path=self.whitelist_path,
    #             ),
    #             call(
    #                 intron_s_path,
    #                 inspect_intron_path,
    #                 whitelist_path=self.whitelist_path,
    #             )
    #         ])
    #         copy_or_create_whitelist.assert_not_called()
    #         bustools_correct.assert_called_once_with(
    #             bus_s_path, bus_sc_path, self.whitelist_path
    #         )
    #         self.assertEqual(2, bustools_count.call_count)
    #         bustools_count.assert_has_calls([
    #             call(
    #                 cdna_s_path,
    #                 os.path.join(counts_dir, BUS_CDNA_PREFIX),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 cm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             ),
    #             call(
    #                 intron_s_path,
    #                 os.path.join(counts_dir, BUS_INTRON_PREFIX),
    #                 self.t2g_path,
    #                 ecmap_path,
    #                 txnames_path,
    #                 tcc=False,
    #                 mm=False,
    #                 cm=False,
    #                 umi_gene=False,
    #                 em=False,
    #             )
    #         ])
    #         filter_with_bustools.assert_not_called()
    #         convert_matrices.assert_not_called()
    # 
    #         STATS.start.assert_called_once()
    #         STATS.end.assert_called_once()
    #         STATS.to_dict.assert_not_called()
    #         import_matrix_as_anndata.assert_not_called()
    #         render_report.assert_not_called()
