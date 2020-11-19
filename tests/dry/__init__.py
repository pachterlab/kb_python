from unittest import mock, TestCase

import kb_python.dry as dry
from tests.mixins import TestMixin


def dry_dummy_function(i):
    return i + 1


@dry.dryable(dry_dummy_function)
def dummy_function(i):
    return i


class TestDry(TestMixin, TestCase):

    def test_dummy_function(self):
        dry.dummy_function()

    def test_dryable_not_dry(self):
        with mock.patch('kb_python.dry.is_dry') as is_dry:
            is_dry.return_value = False
            self.assertEqual(1, dummy_function(1))

    def test_dryable_dry(self):
        with mock.patch('kb_python.dry.is_dry') as is_dry:
            is_dry.return_value = True
            self.assertEqual(2, dummy_function(1))
