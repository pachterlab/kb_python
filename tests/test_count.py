import os
import tempfile
from unittest import TestCase

import kb_python
from kb_python.constants import (
    BUS_S_FILENAME,
    BUS_SC_FILENAME,
    COUNTS_DIR,
    COUNTS_PREFIX,
)
from tests.mixins import TestMixin


class TestCount(TestMixin, TestCase):

    def test_kallisto_bus(self):
        out_dir = tempfile.mkdtemp()
        result = kb_python.kallisto_bus(
            self.fastqs, self.index_path, self.technology, out_dir, threads=1
        )
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_sort(self):
        out_dir = tempfile.mkdtemp()
        out_path = os.path.join(out_dir, BUS_S_FILENAME)
        result = kb_python.bustools_sort(
            self.bus_path, out_path, threads=1, memory='1G'
        )
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_correct(self):
        out_dir = tempfile.mkdtemp()
        out_path = os.path.join(out_dir, BUS_SC_FILENAME)
        result = kb_python.bustools_correct(
            self.bus_s_path, out_path, self.whitelist_path
        )
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_count(self):
        out_dir = tempfile.mkdtemp()
        counts_dir = os.path.join(out_dir, COUNTS_DIR)
        os.makedirs(counts_dir, exist_ok=True)
        counts_path = os.path.join(counts_dir, COUNTS_PREFIX)
        result = kb_python.bustools_count(
            self.bus_scs_path, counts_path, self.t2g_path, self.ecmap_path,
            self.txnames_path
        )
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_bustools_whitelist(self):
        out_dir = tempfile.mkdtemp()
        out_path = os.path.join(out_dir, 'whitelist.txt')
        result = kb_python.bustools_whitelist(self.bus_s_path, out_path)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))
