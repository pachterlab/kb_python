import os
from unittest import mock, TestCase

import kb_python.config as config


class TestConfig(TestCase):

    def test_get_provided_kallisto_path(self):
        with mock.patch('kb_python.config.PLATFORM', 'linux'),\
            mock.patch('kb_python.config.os.path.isfile', return_value=True):
            self.assertTrue(
                config.get_provided_kallisto_path().endswith('kallisto')
            )

    def test_get_provided_kallisto_path_windows(self):
        with mock.patch('kb_python.config.PLATFORM', 'windows'),\
            mock.patch('kb_python.config.os.path.isfile', return_value=True):
            self.assertTrue(
                config.get_provided_kallisto_path().endswith('kallisto.exe')
            )

    def test_get_provided_kallisto_path_unsupported(self):
        with mock.patch('kb_python.config.os.path.isfile', return_value=False):
            self.assertEqual(None, config.get_provided_kallisto_path())

    def test_set_kallisto_binary_path(self):
        with mock.patch('kb_python.config.KALLISTO_PATH', None),\
            mock.patch('kb_python.config.shutil.which', return_value=None),\
            mock.patch('kb_python.config.os.path.isfile', return_value=True),\
            mock.patch('kb_python.config.os.access', return_value=True):
            config.set_kallisto_binary_path('kallisto')
            self.assertEqual(
                os.path.abspath('kallisto'), config.get_kallisto_binary_path()
            )

    def test_set_kallisto_binary_path_in_path(self):
        with mock.patch('kb_python.config.KALLISTO_PATH', None),\
            mock.patch('kb_python.config.shutil.which', return_value='path/to/kallisto'),\
            mock.patch('kb_python.config.os.access', return_value=True):
            config.set_kallisto_binary_path('kallisto')
            self.assertEqual(
                os.path.abspath('path/to/kallisto'),
                config.get_kallisto_binary_path()
            )

    def test_set_kallisto_binary_path_not_isfile(self):
        with mock.patch('kb_python.config.KALLISTO_PATH', None),\
            mock.patch('kb_python.config.shutil.which', return_value=None),\
            mock.patch('kb_python.config.os.path.isfile', return_value=False):
            with self.assertRaises(Exception):
                config.set_kallisto_binary_path('kallisto')

    def test_set_kallisto_binary_path_not_executable(self):
        with mock.patch('kb_python.config.KALLISTO_PATH', None),\
            mock.patch('kb_python.config.shutil.which', return_value='path/to/kallisto'),\
            mock.patch('kb_python.config.os.access', return_value=False):
            with self.assertRaises(config.ConfigError):
                config.set_kallisto_binary_path('kallisto')

    def test_get_provided_bustools_path(self):
        with mock.patch('kb_python.config.PLATFORM', 'linux'),\
            mock.patch('kb_python.config.os.path.isfile', return_value=True):
            self.assertTrue(
                config.get_provided_bustools_path().endswith('bustools')
            )

    def test_get_provided_bustools_path_windows(self):
        with mock.patch('kb_python.config.PLATFORM', 'windows'),\
            mock.patch('kb_python.config.os.path.isfile', return_value=True):
            self.assertTrue(
                config.get_provided_bustools_path().endswith('bustools.exe')
            )

    def test_get_provided_bustools_path_unsupported(self):
        with mock.patch('kb_python.config.os.path.isfile', return_value=False):
            self.assertEqual(None, config.get_provided_bustools_path())

    def test_set_bustools_binary_path(self):
        with mock.patch('kb_python.config.BUSTOOLS_PATH', None),\
            mock.patch('kb_python.config.shutil.which', return_value=None),\
            mock.patch('kb_python.config.os.path.isfile', return_value=True),\
            mock.patch('kb_python.config.os.access', return_value=True):
            config.set_bustools_binary_path('bustools')
            self.assertEqual(
                os.path.abspath('bustools'), config.get_bustools_binary_path()
            )

    def test_set_bustools_binary_path_in_path(self):
        with mock.patch('kb_python.config.BUSTOOLS_PATH', None),\
            mock.patch('kb_python.config.shutil.which', return_value='path/to/bustools'),\
            mock.patch('kb_python.config.os.access', return_value=True):
            config.set_bustools_binary_path('bustools')
            self.assertEqual(
                os.path.abspath('path/to/bustools'),
                config.get_bustools_binary_path()
            )

    def test_set_bustools_binary_path_not_isfile(self):
        with mock.patch('kb_python.config.BUSTOOLS_PATH', None),\
            mock.patch('kb_python.config.shutil.which', return_value=None),\
            mock.patch('kb_python.config.os.path.isfile', return_value=False):
            with self.assertRaises(Exception):
                config.set_bustools_binary_path('bustools')

    def test_set_bustools_binary_path_not_executable(self):
        with mock.patch('kb_python.config.BUSTOOLS_PATH', None),\
            mock.patch('kb_python.config.shutil.which', return_value='path/to/bustools'),\
            mock.patch('kb_python.config.os.access', return_value=False):
            with self.assertRaises(config.ConfigError):
                config.set_bustools_binary_path('bustools')
