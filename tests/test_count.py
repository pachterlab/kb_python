import os
import tempfile
from unittest import mock, TestCase
from unittest.mock import call

import kb_python.count as count
from kb_python.constants import (
    BUS_CDNA_PREFIX,
    BUS_FILENAME,
    BUS_INTRON_PREFIX,
    BUS_S_FILENAME,
    BUS_SC_FILENAME,
    BUS_SCS_FILENAME,
    COUNTS_DIR,
    COUNTS_PREFIX,
    ECMAP_FILENAME,
    INSPECT_FILENAME,
    TXNAMES_FILENAME,
    WHITELIST_FILENAME,
)
from tests.mixins import TestMixin


class TestCount(TestMixin, TestCase):

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
        counts_dir = os.path.join(out_dir, COUNTS_DIR)
        os.makedirs(counts_dir, exist_ok=True)
        counts_path = os.path.join(counts_dir, COUNTS_PREFIX)
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
            self.velocity_bus_scs_path, out_path, self.velocity_cdna_t2c_path,
            self.velocity_txnames_path, self.velocity_txnames_path
        )
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_whitelist(self):
        out_dir = tempfile.mkdtemp()
        out_path = os.path.join(out_dir, 'whitelist.txt')
        result = count.bustools_whitelist(self.bus_s_path, out_path)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

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
        with mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count:
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

    def test_count_without_whitelist(self):
        with mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_count') as bustools_count:
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

    def test_count_velocity_with_whitelist(self):
        with mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
            mock.patch('kb_python.count.bustools_count') as bustools_count:
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

    def test_count_velocity_without_whitelist(self):
        with mock.patch('kb_python.count.kallisto_bus') as kallisto_bus,\
            mock.patch('kb_python.count.bustools_sort') as bustools_sort,\
            mock.patch('kb_python.count.bustools_inspect') as bustools_inspect,\
            mock.patch('kb_python.count.copy_or_create_whitelist') as copy_or_create_whitelist,\
            mock.patch('kb_python.count.bustools_correct') as bustools_correct,\
            mock.patch('kb_python.count.bustools_capture') as bustools_capture,\
            mock.patch('kb_python.count.bustools_count') as bustools_count:
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
