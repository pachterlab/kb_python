import os
import tempfile
import uuid
from unittest import TestCase

import kb_python.ref as ref
from tests.mixins import TestMixin


class TestRef(TestMixin, TestCase):

    def test_kallisto_index(self):
        index_path = os.path.join(
            tempfile.gettempdir(), '{}.idx'.format(uuid.uuid4())
        )
        self.assertFalse(os.path.exists(index_path))
        result = ref.kallisto_index(self.fasta_path, index_path)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_create_t2g(self):
        t2g_path = os.path.join(
            tempfile.gettempdir(), '{}.t2g'.format(uuid.uuid4())
        )
        self.assertFalse(os.path.exists(t2g_path))
        result = ref.create_t2g(self.gtf_path, t2g_path)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))
