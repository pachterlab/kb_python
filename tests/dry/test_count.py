from unittest import TestCase

import kb_python.dry.count as count


class TestCount(TestCase):

    def test_stream_batch(self):
        with self.assertRaises(Exception):
            count.stream_batch(self.smartseq3_remote_batch_path, self.temp_dir)
