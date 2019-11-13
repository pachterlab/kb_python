from unittest import TestCase

import kb_python.dry.count as count


class TestCount(TestCase):

    def test_count_bus_records(self):
        self.assertTrue(count.count_bus_records('path'))
