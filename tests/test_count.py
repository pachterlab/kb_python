import os
import tempfile
import uuid
from unittest import mock, TestCase
from unittest.mock import call

import kb_python.count as count
from kb_python.constants import (
    ADATA_PREFIX,
    BUS_CDNA_PREFIX,
    BUS_FILENAME,
    BUS_INTRON_PREFIX,
    BUS_S_FILENAME,
    BUS_SC_FILENAME,
    BUS_SCS_FILENAME,
    COUNTS_DIR,
    COUNTS_PREFIX,
    ECMAP_FILENAME,
    FILTERED_DIR,
    FILTER_WHITELIST_FILENAME,
    INSPECT_FILENAME,
    TXNAMES_FILENAME,
    UNFILTERED_DIR,
    WHITELIST_FILENAME,
)
from tests.mixins import TestMixin


class TestCount(TestMixin, TestCase):

    def setUp(self):
        makedirs_mock = mock.patch('kb_python.count.os.makedirs')
        makedirs_mock.start()
        self.addCleanup(makedirs_mock.stop)

    def test_kallisto_bus(self):
        out_dir = tempfile.mkdtemp()
        result = count.kallisto_bus(
            self.fastqs, self.index_path, self.technology, out_dir, threads=1
        )
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_sort(self):
        out_dir = tempfile.mkdtemp()
        out_path = os.path.join(out_dir, BUS_S_FILENAME)
        result = count.bustools_sort(
            self.bus_path, out_path, threads=1, memory='1G'
        )
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_inspect(self):
        out_dir = tempfile.mkdtemp()
        out_path = os.path.join(out_dir, INSPECT_FILENAME)
        result = count.bustools_inspect(
            self.bus_s_path, out_path, self.whitelist_path, self.ecmap_path
        )
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_correct(self):
        out_dir = tempfile.mkdtemp()
        out_path = os.path.join(out_dir, BUS_SC_FILENAME)
        result = count.bustools_correct(
            self.bus_s_path, out_path, self.whitelist_path
        )
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_count(self):
        out_dir = tempfile.mkdtemp()
        counts_path = os.path.join(out_dir, COUNTS_PREFIX)
        result = count.bustools_count(
            self.bus_scs_path, counts_path, self.t2g_path, self.ecmap_path,
            self.txnames_path
        )
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_capture(self):
        out_dir = tempfile.mkdtemp()
        out_path = os.path.join(out_dir, 'capture.bus')
        result = count.bustools_capture(
            self.lamanno_bus_scs_path, out_path, self.lamanno_cdna_t2c_path,
            self.lamanno_txnames_path, self.lamanno_txnames_path
        )
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_whitelist(self):
        out_dir = tempfile.mkdtemp()
        out_path = os.path.join(out_dir, 'whitelist.txt')
        result = count.bustools_whitelist(self.bus_s_path, out_path)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_convert_matrix_to_loom(self):
        out_path = os.path.join(
            tempfile.mkdtemp(), '{}.loom'.format(uuid.uuid4())
        )
        count.convert_matrix_to_loom(
            self.matrix_path, self.barcodes_path, self.genes_path, out_path
        )
        self.assertTrue(os.path.exists(out_path))

    def test_convert_matrix_to_h5ad(self):
        out_path = os.path.join(
            tempfile.mkdtemp(), '{}.h5ad'.format(uuid.uuid4())
        )
        count.convert_matrix_to_h5ad(
            self.matrix_path, self.barcodes_path, self.genes_path, out_path
        )
        self.assertTrue(os.path.exists(out_path))

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

    def test_count_with_whitelist(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix_to_loom') as convert_matrix_to_loom,\
            mock.patch('kb_python.count.convert_matrix_to_h5ad') as convert_matrix_to_h5ad:
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_prefix = os.path.join(out_dir, COUNTS_DIR, COUNTS_PREFIX)
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_SCS_FILENAME)
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
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

            self.assertEqual({
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'inspect': inspect_path,
                'bus_scs': bus_scs_path,
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
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
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            bustools_count.assert_called_once_with(
                bus_scs_path, counts_prefix, self.t2g_path, ecmap_path,
                txnames_path
            )
            convert_matrix_to_loom.assert_not_called()
            convert_matrix_to_h5ad.assert_not_called()

    def test_count_loom(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix_to_loom') as convert_matrix_to_loom,\
            mock.patch('kb_python.count.convert_matrix_to_h5ad') as convert_matrix_to_h5ad:
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_prefix = os.path.join(out_dir, COUNTS_DIR, COUNTS_PREFIX)
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_SCS_FILENAME)
            loom_path = mock.MagicMock()
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
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
            convert_matrix_to_loom.return_value = loom_path

            self.assertEqual({
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'inspect': inspect_path,
                'bus_scs': bus_scs_path,
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
                'loom': loom_path,
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
                bus_scs_path, counts_prefix, self.t2g_path, ecmap_path,
                txnames_path
            )
            convert_matrix_to_loom.assert_called_once_with(
                '{}.mtx'.format(counts_prefix),
                '{}.barcodes.txt'.format(counts_prefix),
                '{}.genes.txt'.format(counts_prefix),
                os.path.join(
                    out_dir, COUNTS_DIR, '{}.loom'.format(ADATA_PREFIX)
                )
            )
            convert_matrix_to_h5ad.assert_not_called()

    def test_count_h5ad(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix_to_loom') as convert_matrix_to_loom,\
            mock.patch('kb_python.count.convert_matrix_to_h5ad') as convert_matrix_to_h5ad:
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_prefix = os.path.join(out_dir, COUNTS_DIR, COUNTS_PREFIX)
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_SCS_FILENAME)
            h5ad_path = mock.MagicMock()
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
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
            convert_matrix_to_h5ad.return_value = h5ad_path

            self.assertEqual({
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'inspect': inspect_path,
                'bus_scs': bus_scs_path,
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
                'h5ad': h5ad_path,
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
                                 h5ad=True,
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
                bus_scs_path, counts_prefix, self.t2g_path, ecmap_path,
                txnames_path
            )
            convert_matrix_to_loom.assert_not_called()
            convert_matrix_to_h5ad.assert_called_once_with(
                '{}.mtx'.format(counts_prefix),
                '{}.barcodes.txt'.format(counts_prefix),
                '{}.genes.txt'.format(counts_prefix),
                os.path.join(
                    out_dir, COUNTS_DIR, '{}.h5ad'.format(ADATA_PREFIX)
                )
            )

    def test_count_filter(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix_to_loom') as convert_matrix_to_loom,\
            mock.patch('kb_python.count.convert_matrix_to_h5ad') as convert_matrix_to_h5ad,\
            mock.patch('kb_python.count.bustools_whitelist') as bustools_whitelist,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture:
            out_dir = tempfile.mkdtemp()
            unfiltered_dir = os.path.join(out_dir, UNFILTERED_DIR)
            filtered_dir = os.path.join(out_dir, FILTERED_DIR)
            temp_dir = tempfile.mkdtemp()
            counts_prefix = os.path.join(
                unfiltered_dir, COUNTS_DIR, COUNTS_PREFIX
            )
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(unfiltered_dir, BUS_FILENAME)
            ecmap_path = os.path.join(unfiltered_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(unfiltered_dir, TXNAMES_FILENAME)
            inspect_path = os.path.join(unfiltered_dir, INSPECT_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(unfiltered_dir, BUS_SCS_FILENAME)
            filtered_counts_prefix = os.path.join(
                filtered_dir, COUNTS_DIR, COUNTS_PREFIX
            )
            filter_whitelist_path = os.path.join(
                filtered_dir, FILTER_WHITELIST_FILENAME
            )
            filtered_temp_bus_path = os.path.join(temp_dir, BUS_SCS_FILENAME)
            filtered_bus_path = os.path.join(filtered_dir, BUS_SCS_FILENAME)
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
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

            self.assertEqual({
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'inspect': inspect_path,
                    'bus_scs': bus_scs_path,
                    'mtx': '{}.mtx'.format(counts_prefix),
                    'genes': '{}.genes.txt'.format(counts_prefix),
                    'barcodes': '{}.barcodes.txt'.format(counts_prefix),
                },
                'filtered': {
                    'bus_scs':
                        filtered_bus_path,
                    'whitelist':
                        filter_whitelist_path,
                    'mtx':
                        '{}.mtx'.format(filtered_counts_prefix),
                    'genes':
                        '{}.genes.txt'.format(filtered_counts_prefix),
                    'barcodes':
                        '{}.barcodes.txt'.format(filtered_counts_prefix),
                }
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
                unfiltered_dir,
                threads=threads
            )
            self.assertEqual(bustools_sort.call_count, 3)
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
                    filtered_temp_bus_path,
                    filtered_bus_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path, inspect_path, self.whitelist_path, ecmap_path
            )
            copy_or_create_whitelist.assert_not_called()
            bustools_whitelist.assert_called_once_with(
                bus_scs_path, filter_whitelist_path
            )
            bustools_capture.assert_called_once_with(
                bus_scs_path,
                filtered_temp_bus_path,
                filter_whitelist_path,
                ecmap_path,
                txnames_path,
                capture_type='barcode',
            )
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            self.assertEqual(2, bustools_count.call_count)
            bustools_count.assert_has_calls([
                call(
                    bus_scs_path, counts_prefix, self.t2g_path, ecmap_path,
                    txnames_path
                ),
                call(
                    filtered_bus_path, filtered_counts_prefix, self.t2g_path,
                    ecmap_path, txnames_path
                ),
            ])
            convert_matrix_to_loom.assert_not_called()
            convert_matrix_to_h5ad.assert_not_called()

    def test_count_filter_loom(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix_to_loom') as convert_matrix_to_loom,\
            mock.patch('kb_python.count.convert_matrix_to_h5ad') as convert_matrix_to_h5ad,\
            mock.patch('kb_python.count.bustools_whitelist') as bustools_whitelist,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture:
            out_dir = tempfile.mkdtemp()
            unfiltered_dir = os.path.join(out_dir, UNFILTERED_DIR)
            filtered_dir = os.path.join(out_dir, FILTERED_DIR)
            temp_dir = tempfile.mkdtemp()
            counts_prefix = os.path.join(
                unfiltered_dir, COUNTS_DIR, COUNTS_PREFIX
            )
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(unfiltered_dir, BUS_FILENAME)
            ecmap_path = os.path.join(unfiltered_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(unfiltered_dir, TXNAMES_FILENAME)
            inspect_path = os.path.join(unfiltered_dir, INSPECT_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(unfiltered_dir, BUS_SCS_FILENAME)
            loom_path = mock.MagicMock()
            filtered_counts_prefix = os.path.join(
                filtered_dir, COUNTS_DIR, COUNTS_PREFIX
            )
            filter_whitelist_path = os.path.join(
                filtered_dir, FILTER_WHITELIST_FILENAME
            )
            filtered_temp_bus_path = os.path.join(temp_dir, BUS_SCS_FILENAME)
            filtered_bus_path = os.path.join(filtered_dir, BUS_SCS_FILENAME)
            filtered_loom_path = mock.MagicMock()
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
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
            convert_matrix_to_loom.side_effect = [loom_path, filtered_loom_path]

            self.assertEqual({
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'inspect': inspect_path,
                    'bus_scs': bus_scs_path,
                    'mtx': '{}.mtx'.format(counts_prefix),
                    'genes': '{}.genes.txt'.format(counts_prefix),
                    'barcodes': '{}.barcodes.txt'.format(counts_prefix),
                    'loom': loom_path,
                },
                'filtered': {
                    'bus_scs':
                        filtered_bus_path,
                    'whitelist':
                        filter_whitelist_path,
                    'mtx':
                        '{}.mtx'.format(filtered_counts_prefix),
                    'genes':
                        '{}.genes.txt'.format(filtered_counts_prefix),
                    'barcodes':
                        '{}.barcodes.txt'.format(filtered_counts_prefix),
                    'loom':
                        filtered_loom_path,
                }
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
                unfiltered_dir,
                threads=threads
            )
            self.assertEqual(bustools_sort.call_count, 3)
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
                    filtered_temp_bus_path,
                    filtered_bus_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path, inspect_path, self.whitelist_path, ecmap_path
            )
            copy_or_create_whitelist.assert_not_called()
            bustools_whitelist.assert_called_once_with(
                bus_scs_path, filter_whitelist_path
            )
            bustools_capture.assert_called_once_with(
                bus_scs_path,
                filtered_temp_bus_path,
                filter_whitelist_path,
                ecmap_path,
                txnames_path,
                capture_type='barcode',
            )
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            self.assertEqual(2, bustools_count.call_count)
            bustools_count.assert_has_calls([
                call(
                    bus_scs_path, counts_prefix, self.t2g_path, ecmap_path,
                    txnames_path
                ),
                call(
                    filtered_bus_path, filtered_counts_prefix, self.t2g_path,
                    ecmap_path, txnames_path
                ),
            ])
            self.assertEqual(2, convert_matrix_to_loom.call_count)
            convert_matrix_to_loom.assert_has_calls([
                call(
                    '{}.mtx'.format(counts_prefix),
                    '{}.barcodes.txt'.format(counts_prefix),
                    '{}.genes.txt'.format(counts_prefix),
                    os.path.join(
                        unfiltered_dir, COUNTS_DIR,
                        '{}.loom'.format(ADATA_PREFIX)
                    )
                ),
                call(
                    '{}.mtx'.format(filtered_counts_prefix),
                    '{}.barcodes.txt'.format(filtered_counts_prefix),
                    '{}.genes.txt'.format(filtered_counts_prefix),
                    os.path.join(
                        filtered_dir, COUNTS_DIR,
                        '{}.loom'.format(ADATA_PREFIX)
                    )
                )
            ])
            convert_matrix_to_h5ad.assert_not_called()

    def test_count_filter_h5ad(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix_to_loom') as convert_matrix_to_loom,\
            mock.patch('kb_python.count.convert_matrix_to_h5ad') as convert_matrix_to_h5ad,\
            mock.patch('kb_python.count.bustools_whitelist') as bustools_whitelist,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture:
            out_dir = tempfile.mkdtemp()
            unfiltered_dir = os.path.join(out_dir, UNFILTERED_DIR)
            filtered_dir = os.path.join(out_dir, FILTERED_DIR)
            temp_dir = tempfile.mkdtemp()
            counts_prefix = os.path.join(
                unfiltered_dir, COUNTS_DIR, COUNTS_PREFIX
            )
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(unfiltered_dir, BUS_FILENAME)
            ecmap_path = os.path.join(unfiltered_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(unfiltered_dir, TXNAMES_FILENAME)
            inspect_path = os.path.join(unfiltered_dir, INSPECT_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(unfiltered_dir, BUS_SCS_FILENAME)
            h5ad_path = mock.MagicMock()
            filtered_counts_prefix = os.path.join(
                filtered_dir, COUNTS_DIR, COUNTS_PREFIX
            )
            filter_whitelist_path = os.path.join(
                filtered_dir, FILTER_WHITELIST_FILENAME
            )
            filtered_temp_bus_path = os.path.join(temp_dir, BUS_SCS_FILENAME)
            filtered_bus_path = os.path.join(filtered_dir, BUS_SCS_FILENAME)
            filtered_h5ad_path = mock.MagicMock()
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
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
            convert_matrix_to_h5ad.side_effect = [h5ad_path, filtered_h5ad_path]

            self.assertEqual({
                'unfiltered': {
                    'bus': bus_path,
                    'ecmap': ecmap_path,
                    'txnames': txnames_path,
                    'inspect': inspect_path,
                    'bus_scs': bus_scs_path,
                    'mtx': '{}.mtx'.format(counts_prefix),
                    'genes': '{}.genes.txt'.format(counts_prefix),
                    'barcodes': '{}.barcodes.txt'.format(counts_prefix),
                    'h5ad': h5ad_path,
                },
                'filtered': {
                    'bus_scs':
                        filtered_bus_path,
                    'whitelist':
                        filter_whitelist_path,
                    'mtx':
                        '{}.mtx'.format(filtered_counts_prefix),
                    'genes':
                        '{}.genes.txt'.format(filtered_counts_prefix),
                    'barcodes':
                        '{}.barcodes.txt'.format(filtered_counts_prefix),
                    'h5ad':
                        filtered_h5ad_path,
                }
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
                                 memory=memory,
                                 h5ad=True
                             ))
            stream_fastqs.assert_called_once_with(
                self.fastqs, temp_dir=temp_dir
            )
            kallisto_bus.assert_called_once_with(
                self.fastqs,
                self.index_path,
                self.technology,
                unfiltered_dir,
                threads=threads
            )
            self.assertEqual(bustools_sort.call_count, 3)
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
                    filtered_temp_bus_path,
                    filtered_bus_path,
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )
            ])
            bustools_inspect.assert_called_once_with(
                bus_s_path, inspect_path, self.whitelist_path, ecmap_path
            )
            copy_or_create_whitelist.assert_not_called()
            bustools_whitelist.assert_called_once_with(
                bus_scs_path, filter_whitelist_path
            )
            bustools_capture.assert_called_once_with(
                bus_scs_path,
                filtered_temp_bus_path,
                filter_whitelist_path,
                ecmap_path,
                txnames_path,
                capture_type='barcode',
            )
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            self.assertEqual(2, bustools_count.call_count)
            bustools_count.assert_has_calls([
                call(
                    bus_scs_path, counts_prefix, self.t2g_path, ecmap_path,
                    txnames_path
                ),
                call(
                    filtered_bus_path, filtered_counts_prefix, self.t2g_path,
                    ecmap_path, txnames_path
                ),
            ])
            self.assertEqual(2, convert_matrix_to_h5ad.call_count)
            convert_matrix_to_h5ad.assert_has_calls([
                call(
                    '{}.mtx'.format(counts_prefix),
                    '{}.barcodes.txt'.format(counts_prefix),
                    '{}.genes.txt'.format(counts_prefix),
                    os.path.join(
                        unfiltered_dir, COUNTS_DIR,
                        '{}.h5ad'.format(ADATA_PREFIX)
                    )
                ),
                call(
                    '{}.mtx'.format(filtered_counts_prefix),
                    '{}.barcodes.txt'.format(filtered_counts_prefix),
                    '{}.genes.txt'.format(filtered_counts_prefix),
                    os.path.join(
                        filtered_dir, COUNTS_DIR,
                        '{}.h5ad'.format(ADATA_PREFIX)
                    )
                )
            ])
            convert_matrix_to_loom.assert_not_called()

    def test_count_without_whitelist(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.convert_matrix_to_loom') as convert_matrix_to_loom,\
            mock.patch('kb_python.count.convert_matrix_to_h5ad') as convert_matrix_to_h5ad:
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_prefix = os.path.join(out_dir, COUNTS_DIR, COUNTS_PREFIX)
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_SCS_FILENAME)
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
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

            self.assertEqual({
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'whitelist': self.whitelist_path,
                'inspect': inspect_path,
                'bus_scs': bus_scs_path,
                'mtx': '{}.mtx'.format(counts_prefix),
                'genes': '{}.genes.txt'.format(counts_prefix),
                'barcodes': '{}.barcodes.txt'.format(counts_prefix),
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
                bus_scs_path, counts_prefix, self.t2g_path, ecmap_path,
                txnames_path
            )
            convert_matrix_to_loom.assert_not_called()
            convert_matrix_to_h5ad.assert_not_called()

    def test_count_lamanno_with_whitelist(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata,\
            mock.patch('kb_python.count.overlay_anndatas') as overlay_anndatas:
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_dir = os.path.join(out_dir, COUNTS_DIR)
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_SCS_FILENAME)
            cdna_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
            )
            intron_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
            )
            cdna_s_path = os.path.join(
                out_dir, '{}.s.bus'.format(BUS_CDNA_PREFIX)
            )
            intron_s_path = os.path.join(
                out_dir, '{}.s.bus'.format(BUS_INTRON_PREFIX)
            )
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
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
            bustools_inspect.return_value = {'inspect': inspect_path}
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

            self.assertEqual({
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'inspect': inspect_path,
                'bus_scs': bus_scs_path,
                BUS_CDNA_PREFIX: {
                    'bus_s':
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
                },
                BUS_INTRON_PREFIX: {
                    'bus_s':
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
                }
            },
                             count.count_lamanno(
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
            bustools_inspect.assert_called_once_with(
                bus_s_path, inspect_path, self.whitelist_path, ecmap_path
            )
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            self.assertEqual(2, bustools_count.call_count)
            bustools_count.assert_has_calls([
                call(
                    cdna_s_path, os.path.join(counts_dir, BUS_CDNA_PREFIX),
                    self.t2g_path, ecmap_path, txnames_path
                ),
                call(
                    intron_s_path, os.path.join(counts_dir, BUS_INTRON_PREFIX),
                    self.t2g_path, ecmap_path, txnames_path
                )
            ])
            import_matrix_as_anndata.assert_not_called()
            overlay_anndatas.assert_not_called()

    def test_count_lamanno_loom(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata,\
            mock.patch('kb_python.count.overlay_anndatas') as overlay_anndatas:
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_dir = os.path.join(out_dir, COUNTS_DIR)
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_SCS_FILENAME)
            cdna_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
            )
            intron_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
            )
            cdna_s_path = os.path.join(
                out_dir, '{}.s.bus'.format(BUS_CDNA_PREFIX)
            )
            intron_s_path = os.path.join(
                out_dir, '{}.s.bus'.format(BUS_INTRON_PREFIX)
            )
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            adata_spliced = mock.MagicMock()
            adata_unspliced = mock.MagicMock()
            adata = mock.MagicMock()
            loom_path = os.path.join(counts_dir, '{}.loom'.format(ADATA_PREFIX))
            adata.write_loom.return_value = loom_path
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
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
            bustools_inspect.return_value = {'inspect': inspect_path}
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
            import_matrix_as_anndata.side_effect = [
                adata_spliced, adata_unspliced
            ]
            overlay_anndatas.return_value = adata

            self.assertEqual({
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'inspect': inspect_path,
                'bus_scs': bus_scs_path,
                'loom': loom_path,
                BUS_CDNA_PREFIX: {
                    'bus_s':
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
                },
                BUS_INTRON_PREFIX: {
                    'bus_s':
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
                }
            },
                             count.count_lamanno(
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
            bustools_inspect.assert_called_once_with(
                bus_s_path, inspect_path, self.whitelist_path, ecmap_path
            )
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            self.assertEqual(2, bustools_count.call_count)
            bustools_count.assert_has_calls([
                call(
                    cdna_s_path, os.path.join(counts_dir, BUS_CDNA_PREFIX),
                    self.t2g_path, ecmap_path, txnames_path
                ),
                call(
                    intron_s_path, os.path.join(counts_dir, BUS_INTRON_PREFIX),
                    self.t2g_path, ecmap_path, txnames_path
                )
            ])
            self.assertEqual(2, import_matrix_as_anndata.call_count)
            import_matrix_as_anndata.assert_has_calls([
                call(
                    '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ), '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    )
                ),
                call(
                    '{}.mtx'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ), '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ), '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    )
                ),
            ])
            overlay_anndatas.assert_called_once_with(
                adata_spliced, adata_unspliced
            )
            adata.write_loom.assert_called_once_with(loom_path)
            adata.write.assert_not_called()

    def test_count_lamanno_h5ad(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata,\
            mock.patch('kb_python.count.overlay_anndatas') as overlay_anndatas:
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_dir = os.path.join(out_dir, COUNTS_DIR)
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_SCS_FILENAME)
            cdna_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
            )
            intron_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
            )
            cdna_s_path = os.path.join(
                out_dir, '{}.s.bus'.format(BUS_CDNA_PREFIX)
            )
            intron_s_path = os.path.join(
                out_dir, '{}.s.bus'.format(BUS_INTRON_PREFIX)
            )
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            adata_spliced = mock.MagicMock()
            adata_unspliced = mock.MagicMock()
            adata = mock.MagicMock()
            h5ad_path = os.path.join(counts_dir, '{}.h5ad'.format(ADATA_PREFIX))
            adata.write.return_value = h5ad_path
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
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
            bustools_inspect.return_value = {'inspect': inspect_path}
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
            import_matrix_as_anndata.side_effect = [
                adata_spliced, adata_unspliced
            ]
            overlay_anndatas.return_value = adata

            self.assertEqual({
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'inspect': inspect_path,
                'bus_scs': bus_scs_path,
                'h5ad': h5ad_path,
                BUS_CDNA_PREFIX: {
                    'bus_s':
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
                },
                BUS_INTRON_PREFIX: {
                    'bus_s':
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
                }
            },
                             count.count_lamanno(
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
                                 h5ad=True
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
            bustools_inspect.assert_called_once_with(
                bus_s_path, inspect_path, self.whitelist_path, ecmap_path
            )
            copy_or_create_whitelist.assert_not_called()
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            self.assertEqual(2, bustools_count.call_count)
            bustools_count.assert_has_calls([
                call(
                    cdna_s_path, os.path.join(counts_dir, BUS_CDNA_PREFIX),
                    self.t2g_path, ecmap_path, txnames_path
                ),
                call(
                    intron_s_path, os.path.join(counts_dir, BUS_INTRON_PREFIX),
                    self.t2g_path, ecmap_path, txnames_path
                )
            ])
            self.assertEqual(2, import_matrix_as_anndata.call_count)
            import_matrix_as_anndata.assert_has_calls([
                call(
                    '{}.mtx'.format(os.path.join(counts_dir, BUS_CDNA_PREFIX)),
                    '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    ), '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_CDNA_PREFIX)
                    )
                ),
                call(
                    '{}.mtx'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ), '{}.barcodes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    ), '{}.genes.txt'.format(
                        os.path.join(counts_dir, BUS_INTRON_PREFIX)
                    )
                ),
            ])
            overlay_anndatas.assert_called_once_with(
                adata_spliced, adata_unspliced
            )
            adata.write_loom.assert_not_called()
            adata.write.assert_called_once_with(h5ad_path)

    def test_count_lamanno_without_whitelist(self):
        with mock.patch('kb_python.count.stream_fastqs') as stream_fastqs,\
            mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
            mock.patch('kb_python.count.bustools_count') as bustools_count,\
            mock.patch('kb_python.count.import_matrix_as_anndata') as import_matrix_as_anndata,\
            mock.patch('kb_python.count.overlay_anndatas') as overlay_anndatas:
            out_dir = tempfile.mkdtemp()
            temp_dir = tempfile.mkdtemp()
            counts_dir = os.path.join(out_dir, COUNTS_DIR)
            threads = 99999
            memory = 'TEST'
            bus_path = os.path.join(out_dir, BUS_FILENAME)
            ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
            txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
            inspect_path = os.path.join(out_dir, INSPECT_FILENAME)
            bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
            bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
            bus_scs_path = os.path.join(out_dir, BUS_SCS_FILENAME)
            cdna_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_CDNA_PREFIX)
            )
            intron_capture_path = os.path.join(
                temp_dir, '{}.bus'.format(BUS_INTRON_PREFIX)
            )
            cdna_s_path = os.path.join(
                out_dir, '{}.s.bus'.format(BUS_CDNA_PREFIX)
            )
            intron_s_path = os.path.join(
                out_dir, '{}.s.bus'.format(BUS_INTRON_PREFIX)
            )
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            stream_fastqs.return_value = self.fastqs
            kallisto_bus.return_value = {
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
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
            bustools_inspect.return_value = {'inspect': inspect_path}
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

            self.assertEqual({
                'bus': bus_path,
                'ecmap': ecmap_path,
                'txnames': txnames_path,
                'inspect': inspect_path,
                'bus_scs': bus_scs_path,
                'whitelist': self.whitelist_path,
                BUS_CDNA_PREFIX: {
                    'bus_s':
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
                },
                BUS_INTRON_PREFIX: {
                    'bus_s':
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
                }
            },
                             count.count_lamanno(
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
            bustools_inspect.assert_called_once_with(
                bus_s_path, inspect_path, self.whitelist_path, ecmap_path
            )
            copy_or_create_whitelist.assert_called_once_with(
                self.technology, bus_s_path, out_dir
            )
            bustools_correct.assert_called_once_with(
                bus_s_path, bus_sc_path, self.whitelist_path
            )
            self.assertEqual(2, bustools_count.call_count)
            bustools_count.assert_has_calls([
                call(
                    cdna_s_path, os.path.join(counts_dir, BUS_CDNA_PREFIX),
                    self.t2g_path, ecmap_path, txnames_path
                ),
                call(
                    intron_s_path, os.path.join(counts_dir, BUS_INTRON_PREFIX),
                    self.t2g_path, ecmap_path, txnames_path
                )
            ])
            import_matrix_as_anndata.assert_not_called()
            overlay_anndatas.assert_not_called()
