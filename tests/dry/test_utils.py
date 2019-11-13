from unittest import mock, TestCase

import kb_python.dry.utils as utils
from kb_python.config import UnsupportedOSException


class TestUtils(TestCase):

    def test_run_executable(self):
        with mock.patch('kb_python.dry.utils.print') as p:
            self.assertIsNone(utils.run_executable(['1', '2', 3]))
            p.assert_called_once_with('1 2 3')

    def test_run_executable_quiet(self):
        with mock.patch('kb_python.dry.utils.print') as p:
            self.assertIsNone(utils.run_executable(['1', '2', 3], quiet=True))
            p.assert_not_called()

    def test_make_directory(self):
        with mock.patch('kb_python.dry.utils.print') as p:
            self.assertIsNone(utils.make_directory('path'))
            p.assert_called()

    def test_remove_directory(self):
        with mock.patch('kb_python.dry.utils.print') as p:
            self.assertIsNone(utils.remove_directory('path'))
            p.assert_called()

    def test_stream_file(self):
        with mock.patch('kb_python.dry.utils.print') as p,\
            mock.patch('kb_python.dry.utils.PLATFORM', 'darwin'):
            self.assertIsNone(utils.stream_file('url', 'path'))
            p.assert_called()

    def test_stream_file_windows(self):
        with mock.patch('kb_python.dry.utils.print') as p,\
            mock.patch('kb_python.dry.utils.PLATFORM', 'windows'):
            with self.assertRaises(UnsupportedOSException):
                utils.stream_file('url', 'path')
            p.assert_not_called()
