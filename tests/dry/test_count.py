from unittest import TestCase

import kb_python.dry.count as count
from tests.mixins import TestMixin


class TestCount(TestMixin, TestCase):

    def test_stream_batch(self):
        with self.assertRaises(Exception):
            count.stream_batch(self.smartseq3_remote_batch_path, self.temp_dir)
        self.assertEqual(
            self.smartseq3_paired_batch_path,
            count.stream_batch(self.smartseq3_paired_batch_path, self.temp_dir)
        )
