import os
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
        with mock.patch('kb_python.dry.utils.print') as p,\
            mock.patch('kb_python.dry.utils.PLATFORM', 'darwin'):
            self.assertIsNone(utils.make_directory('path'))
            p.assert_called()

    def test_make_directory_windows(self):
        with mock.patch('kb_python.dry.utils.print') as p,\
            mock.patch('kb_python.dry.utils.PLATFORM', 'windows'):
            self.assertIsNone(utils.make_directory('path'))
            p.assert_called()

    def test_remove_directory(self):
        with mock.patch('kb_python.dry.utils.print') as p,\
            mock.patch('kb_python.dry.utils.PLATFORM', 'darwin'):
            self.assertIsNone(utils.remove_directory('path'))
            p.assert_called()

    def test_remove_directory_windows(self):
        with mock.patch('kb_python.dry.utils.print') as p,\
            mock.patch('kb_python.dry.utils.PLATFORM', 'windows'):
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

    def test_move_file(self):
        with mock.patch('kb_python.dry.utils.print') as p,\
            mock.patch('kb_python.dry.utils.PLATFORM', 'darwin'):
            utils.move_file('source', 'destination')
            p.assert_called()

    def test_move_file_windows(self):
        with mock.patch('kb_python.dry.utils.print') as p,\
            mock.patch('kb_python.dry.utils.PLATFORM', 'windows'):
            utils.move_file('source', 'destination')
            p.assert_called()

    def test_copy_whitelist(self):
        with mock.patch('kb_python.dry.utils.print') as p:
            self.assertEquals(
                'path/10xv2_whitelist.txt',
                utils.copy_whitelist('10xv2', 'path')
            )
            p.assert_called()

    def test_copy_map(self):
        with mock.patch('kb_python.dry.utils.print') as p:
            self.assertEqual(
                'path/10xv3_feature_barcode_map.txt',
                utils.copy_map('10xv3', 'path')
            )
            p.assert_called()

    def test_get_temporary_filename(self):
        temp = utils.get_temporary_filename('temp')
        self.assertTrue(temp.startswith(os.path.join('temp', 'tmp')))
