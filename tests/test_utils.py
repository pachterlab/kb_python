import gzip
import os
import subprocess as sp
import tempfile
import uuid
from unittest import mock, TestCase

import kb_python.utils as utils
from tests.mixins import TestMixin


class TestUtils(TestMixin, TestCase):

    def test_run_executable(self):
        with mock.patch('kb_python.utils.logger.info') as info_mock:
            p = utils.run_executable(['echo', 'TEST'], stream=False)
            info_mock.assert_called_once_with('echo TEST')
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
            debug_mock.assert_called_once_with('TEST')

    def test_run_chain(self):
        ps = utils.run_chain(['echo', 'TEST'], ['grep', 'T'])
        self.assertEqual(ps[1].stdout.read(), 'TEST\n')

    def test_run_chain_fails_single_command(self):
        with self.assertRaises(AssertionError):
            utils.run_chain(['echo', 'TEST'])

    def test_run_chain_raises_exception_when_dead(self):
        with self.assertRaises(sp.SubprocessError):
            utils.run_chain(['sleep', '5'], ['grep', 'TEST'], ['ls'])

    def test_generate_cdna_fasta(self):
        out_path = os.path.join(
            tempfile.gettempdir(), '{}.fa'.format(uuid.uuid4())
        )
        utils.generate_cdna_fasta(
            self.sorted_fasta_path, self.sorted_gtf_path, out_path
        )
        with open(out_path, 'r') as f, open(self.split_cdna_fasta_path,
                                            'r') as split:
            self.assertEqual(f.read(), split.read())

    def test_generate_intron_fasta(self):
        out_path = os.path.join(
            tempfile.gettempdir(), '{}.fa'.format(uuid.uuid4())
        )
        utils.generate_intron_fasta(
            self.sorted_fasta_path, self.sorted_gtf_path, out_path
        )
        with open(out_path, 'r') as f, open(self.split_intron_fasta_path,
                                            'r') as split:
            self.assertEqual(f.read(), split.read())

    def test_get_version(self):
        with mock.patch('kb_python.utils.run_executable') as run_executable:
            run_executable().stdout.read.return_value = 'kallisto 1.2.3'
            self.assertEqual(utils.get_version('kallisto'), (1, 2, 3))

    def test_check_dependencies(self):
        with mock.patch('kb_python.utils.MINIMUM_REQUIREMENTS') as minimum_requirements,\
            mock.patch('kb_python.utils.get_version') as get_version:
            minimum_requirements.items.return_value = [('TEST', (0, 0, 0))]
            get_version.return_value = (0, 0, 1)
            utils.check_dependencies()

    def test_check_dependencies_fail(self):
        with mock.patch('kb_python.utils.MINIMUM_REQUIREMENTS') as minimum_requirements,\
            mock.patch('kb_python.utils.get_version') as get_version:
            minimum_requirements.items.return_value = [('TEST', (1, 0, 0))]
            get_version.return_value = (0, 0, 1)
            with self.assertRaises(utils.UnmetDependencyException):
                utils.check_dependencies()

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

    def test_import_matrix_as_anndata(self):
        utils.import_matrix_as_anndata(
            self.matrix_path, self.barcodes_path, self.genes_path
        )

    def test_convert_matrix_to_loom(self):
        out_path = os.path.join(
            tempfile.mkdtemp(), '{}.loom'.format(uuid.uuid4())
        )
        utils.convert_matrix_to_loom(
            self.matrix_path, self.barcodes_path, self.genes_path, out_path
        )
        self.assertTrue(os.path.exists(out_path))

    def test_convert_matrix_to_h5ad(self):
        out_path = os.path.join(
            tempfile.mkdtemp(), '{}.h5ad'.format(uuid.uuid4())
        )
        utils.convert_matrix_to_h5ad(
            self.matrix_path, self.barcodes_path, self.genes_path, out_path
        )
        self.assertTrue(os.path.exists(out_path))

    def test_copy_whitelist(self):
        whitelist_path = utils.copy_whitelist('10xv1', tempfile.mkdtemp())
        self.assertTrue(os.path.exists(whitelist_path))

    def test_get_transcripts_from_fasta(self):
        transcripts = utils.get_transcripts_from_fasta(self.cdna_small_path)
        self.assertEqual({
            'ENST00000456328.1',
            'ENST00000450305.2',
            'ENST00000488147.3',
        }, set(transcripts))

    def test_get_transcripts_from_fasta_gzip(self):
        transcripts = utils.get_transcripts_from_fasta(
            self.cdna_small_gzip_path
        )
        self.assertEqual({
            'ENST00000456328.1',
            'ENST00000450305.2',
            'ENST00000488147.3',
        }, set(transcripts))

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
