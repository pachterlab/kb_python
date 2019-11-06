import gzip
import os
import subprocess as sp
import tempfile
import uuid
from unittest import mock, TestCase
from unittest.mock import call

import kb_python.utils as utils
from kb_python.config import UnsupportedOSException
from tests.mixins import TestMixin


class TestUtils(TestMixin, TestCase):

    def test_open_as_text_textfile(self):
        path = os.path.join(
            tempfile.gettempdir(), '{}.txt'.format(uuid.uuid4())
        )
        with utils.open_as_text(path, 'w') as f:
            f.write('TESTING')
        self.assertTrue(os.path.exists(path))
        with utils.open_as_text(path, 'r') as f:
            self.assertEqual(f.read(), 'TESTING')

    def test_open_as_text_gzip(self):
        path = os.path.join(tempfile.gettempdir(), '{}.gz'.format(uuid.uuid4()))
        with utils.open_as_text(path, 'w') as f:
            f.write('TESTING')
        self.assertTrue(os.path.exists(path))
        with utils.open_as_text(path, 'r') as f:
            self.assertEqual(f.read(), 'TESTING')

    def test_run_executable(self):
        p = utils.run_executable(['echo', 'TEST'], stream=False)
        self.assertEqual(p.stdout.read(), 'TEST\n')

    def test_run_exectuable_raises_exception(self):
        with self.assertRaises(sp.SubprocessError):
            utils.run_executable(['bash', 'nonexistent option'])

    def test_run_exectuable_with_returncode(self):
        utils.run_executable(['bash', 'nonexistent option'], returncode=127)

    def test_run_executable_no_wait(self):
        with mock.patch('kb_python.utils.sp') as sp_mock:
            sp_mock.Popen().returncode = 0
            utils.run_executable(['echo', 'TEST'], wait=False)
            sp_mock.Popen().poll.assert_not_called()

    def test_run_executable_with_stream(self):
        with mock.patch('kb_python.utils.logger.debug') as debug_mock:
            utils.run_executable(['echo', 'TEST'], stream=True)
            debug_mock.assert_has_calls([call('TEST')])

    def test_run_chain(self):
        ps = utils.run_chain(['echo', 'TEST'], ['grep', 'T'])
        self.assertEqual(ps[1].stdout.read(), 'TEST\n')

    def test_run_chain_fails_single_command(self):
        with self.assertRaises(AssertionError):
            utils.run_chain(['echo', 'TEST'])

    def test_run_chain_raises_exception_when_dead(self):
        with self.assertRaises(sp.SubprocessError):
            utils.run_chain(['sleep', '5'], ['grep', 'TEST'], ['ls'])

    def test_get_kallisto_version(self):
        with mock.patch('kb_python.utils.run_executable') as run_executable:
            run_executable().stdout.read.return_value = 'kallisto 1.2.3'
            self.assertEqual((1, 2, 3), utils.get_kallisto_version())

    def test_get_bustools_version(self):
        with mock.patch('kb_python.utils.run_executable') as run_executable:
            run_executable().stdout.read.return_value = 'bustools 1.2.3'
            self.assertEqual((1, 2, 3), utils.get_bustools_version())

    def test_parse_technologies(self):
        lines = [
            'short name       description',
            '----------       -----------',
            '10xv1            10x version 1 chemistry',
            '10xv2            10x version 2 chemistry',
        ]
        self.assertEqual(utils.parse_technologies(lines), {'10xv1', '10xv2'})

    def test_get_supported_technologies(self):
        with mock.patch('kb_python.utils.run_executable') as run_executable,\
           mock.patch('kb_python.utils.parse_technologies') as parse_technologies:
            run_executable().stdout = 'TEST'
            utils.get_supported_technologies()
            parse_technologies.assert_called_once_with('TEST')

    def test_whitelist_provided(self):
        self.assertTrue(utils.whitelist_provided('10xv2'))
        self.assertFalse(utils.whitelist_provided('UNSUPPORTED'))

    def test_stream_file(self):
        with mock.patch('kb_python.utils.PLATFORM', 'linux'),\
            mock.patch('kb_python.utils.os') as os,\
            mock.patch('kb_python.utils.threading') as threading,\
            mock.patch('kb_python.utils.urlretrieve') as urlretrieve:
            url = mock.MagicMock()
            path = mock.MagicMock()
            utils.stream_file(url, path)
            os.mkfifo.assert_called_once_with(path)
            threading.Thread.assert_called_once_with(
                target=urlretrieve, args=(url, path), daemon=True
            )
            threading.Thread().start.assert_called_once_with()

    def test_stream_file_windows(self):
        with mock.patch('kb_python.utils.PLATFORM', 'windows'),\
            mock.patch('kb_python.utils.os') as os,\
            mock.patch('kb_python.utils.threading') as threading,\
            mock.patch('kb_python.utils.urlretrieve') as urlretrieve:
            url = mock.MagicMock()
            path = mock.MagicMock()
            with self.assertRaises(UnsupportedOSException):
                utils.stream_file(url, path)
            os.mkfifo.assert_not_called()
            threading.thread.assert_not_called()
            urlretrieve.assert_not_called()

    def test_import_matrix_as_anndata(self):
        utils.import_matrix_as_anndata(
            self.matrix_path, self.barcodes_path, self.genes_path
        )

    def test_overlay_anndatas(self):
        adata_spliced = utils.import_matrix_as_anndata(
            self.spliced_matrix_path, self.spliced_barcodes_path,
            self.spliced_genes_path
        )
        adata_unspliced = utils.import_matrix_as_anndata(
            self.unspliced_matrix_path, self.unspliced_barcodes_path,
            self.unspliced_genes_path
        )
        adata = utils.overlay_anndatas(adata_spliced, adata_unspliced)
        self.assertEqual({'spliced', 'unspliced'}, set(adata.layers.keys()))

    def test_sum_anndatas(self):
        adata_spliced = utils.import_matrix_as_anndata(
            self.spliced_matrix_path, self.spliced_barcodes_path,
            self.spliced_genes_path
        )
        adata_unspliced = utils.import_matrix_as_anndata(
            self.unspliced_matrix_path, self.unspliced_barcodes_path,
            self.unspliced_genes_path
        )
        adata = utils.sum_anndatas(adata_spliced, adata_unspliced)
        self.assertEqual(2.0, adata.X[5, 15])

    def test_copy_whitelist(self):
        whitelist_path = utils.copy_whitelist('10xv1', tempfile.mkdtemp())
        self.assertTrue(os.path.exists(whitelist_path))

    def test_concatenate_files(self):
        temp_dir = tempfile.mkdtemp()
        file1_path = os.path.join(temp_dir, str(uuid.uuid4()))
        file2_path = os.path.join(temp_dir, '{}.gz'.format(uuid.uuid4()))

        with open(file1_path, 'w') as f:
            f.write('TEST1')
        with gzip.open(file2_path, 'wt') as f:
            f.write('TEST2')

        out_path = utils.concatenate_files(
            file1_path,
            file2_path,
            out_path=os.path.join(temp_dir, str(uuid.uuid4())),
            temp_dir=tempfile.mkdtemp()
        )

        with open(out_path, 'r') as f:
            self.assertEqual(f.read(), 'TEST1\nTEST2\n')
