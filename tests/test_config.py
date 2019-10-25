from unittest import mock, TestCase

import kb_python.config as config


class TestConfig(TestCase):

    def test_get_kallisto_binary_path(self):
        with mock.patch('kb_python.config.PLATFORM', 'linux'),\
            mock.patch('kb_python.config.os.path.exists') as exists:
            exists.return_value = True
            self.assertTrue(
                config.get_kallisto_binary_path().endswith('kallisto')
            )

    def test_get_kallisto_binary_path_windows(self):
        with mock.patch('kb_python.config.PLATFORM', 'windows'),\
            mock.patch('kb_python.config.os.path.exists') as exists:
            exists.return_value = True
            self.assertTrue(
                config.get_kallisto_binary_path().endswith('kallisto.exe')
            )

    def test_get_kallisto_binary_path_unsupported(self):
        with mock.patch('kb_python.config.os.path.exists') as exists:
            exists.return_value = False
            with self.assertRaises(config.UnsupportedOSException):
                config.get_kallisto_binary_path()

    def test_get_bustools_binary_path(self):
        with mock.patch('kb_python.config.PLATFORM', 'linux'),\
            mock.patch('kb_python.config.os.path.exists') as exists:
            exists.return_value = True
            self.assertTrue(
                config.get_bustools_binary_path().endswith('bustools')
            )

    def test_get_bustools_binary_path_windows(self):
        with mock.patch('kb_python.config.PLATFORM', 'windows'),\
            mock.patch('kb_python.config.os.path.exists') as exists:
            exists.return_value = True
            self.assertTrue(
                config.get_bustools_binary_path().endswith('bustools.exe')
            )

    def test_get_bustools_binary_path_unsupported(self):
        with mock.patch('kb_python.config.os.path.exists') as exists:
            exists.return_value = False
            with self.assertRaises(config.UnsupportedOSException):
                config.get_bustools_binary_path()
