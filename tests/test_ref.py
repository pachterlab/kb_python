import os
import tempfile
import uuid
from unittest import mock, TestCase
from unittest.mock import call

import kb_python.ref as ref
from kb_python.constants import COMBINED_FILENAME
from tests.mixins import TestMixin


class TestRef(TestMixin, TestCase):

    def test_kallisto_index(self):
        index_path = os.path.join(
            tempfile.gettempdir(), '{}.idx'.format(uuid.uuid4())
        )
        self.assertFalse(os.path.exists(index_path))
        result = ref.kallisto_index(self.fasta_path, index_path)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_create_t2g(self):
        t2g_path = os.path.join(
            tempfile.gettempdir(), '{}.t2g'.format(uuid.uuid4())
        )
        self.assertFalse(os.path.exists(t2g_path))
        result = ref.create_t2g(self.gtf_path, t2g_path)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_create_t2c(self):
        t2c_path = os.path.join(tempfile.mkdtemp(), str(uuid.uuid4()))
        self.assertEqual({
            't2c':
                t2c_path,
            'transcripts':
                ['ENST00000456328.1', 'ENST00000450305.2', 'ENST00000488147.3']
        }, ref.create_t2c(self.cdna_small_path, t2c_path))
        self.assertTrue(os.path.exists(t2c_path))

    def test_ref(self):
        with mock.patch('kb_python.ref.create_t2g') as create_t2g,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = False
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            kallisto_index.return_value = {'index': index_path}
            create_t2g.return_value = {'t2g': t2g_path}
            self.assertEqual({
                't2g': t2g_path,
                'index': index_path,
            }, ref.ref(self.fasta_path, self.gtf_path, index_path, t2g_path))
            create_t2g.assert_called_once_with(self.gtf_path, t2g_path)
            kallisto_index.assert_called_once_with(self.fasta_path, index_path)

    def test_ref_exists(self):
        with mock.patch('kb_python.ref.create_t2g') as create_t2g,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = True
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            kallisto_index.return_value = {'index': index_path}
            create_t2g.return_value = {'t2g': t2g_path}
            self.assertEqual({
                't2g': t2g_path
            }, ref.ref(self.fasta_path, self.gtf_path, index_path, t2g_path))
            create_t2g.assert_called_once_with(self.gtf_path, t2g_path)
            kallisto_index.assert_not_called()

    def test_ref_overwrite(self):
        with mock.patch('kb_python.ref.create_t2g') as create_t2g,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = True
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            kallisto_index.return_value = {'index': index_path}
            create_t2g.return_value = {'t2g': t2g_path}
            self.assertEqual({
                't2g': t2g_path,
                'index': index_path,
            },
                             ref.ref(
                                 self.fasta_path,
                                 self.gtf_path,
                                 index_path,
                                 t2g_path,
                                 overwrite=True
                             ))
            create_t2g.assert_called_once_with(self.gtf_path, t2g_path)
            kallisto_index.assert_called_once_with(self.fasta_path, index_path)

    def test_ref_velocity(self):
        with mock.patch('kb_python.ref.create_t2g') as create_t2g,\
            mock.patch('kb_python.ref.create_t2c') as create_t2c,\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = False
            temp_dir = tempfile.mkdtemp()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            combined_path = mock.MagicMock()
            create_t2g.return_value = {'t2g': t2g_path}
            create_t2c.side_effect = [{
                't2c': cdna_t2c_path,
                'transcripts': ['A'],
            }, {
                't2c': intron_t2c_path,
                'transcripts': ['B'],
            }]
            concatenate_files.return_value = combined_path
            kallisto_index.return_value = {'index': index_path}
            self.assertEqual({
                't2g': t2g_path,
                'cdna_t2c': cdna_t2c_path,
                'intron_t2c': intron_t2c_path,
                'index': index_path,
            },
                             ref.ref_velocity(
                                 self.cdna_path,
                                 self.intron_path,
                                 self.gtf_path,
                                 index_path,
                                 t2g_path,
                                 cdna_t2c_path,
                                 intron_t2c_path,
                                 temp_dir=temp_dir
                             ))
            create_t2g.assert_called_once_with(
                self.gtf_path, t2g_path, transcripts={'A', 'B'}
            )
            self.assertEqual(create_t2c.call_count, 2)
            create_t2c.assert_has_calls([
                call(self.cdna_path, cdna_t2c_path),
                call(self.intron_path, intron_t2c_path)
            ])
            concatenate_files.assert_called_once_with(
                self.cdna_path,
                self.intron_path,
                out_path=os.path.join(temp_dir, COMBINED_FILENAME),
                temp_dir=temp_dir
            )
            kallisto_index.assert_called_once_with(combined_path, index_path)

    def test_ref_velocity_exists(self):
        with mock.patch('kb_python.ref.create_t2g') as create_t2g,\
            mock.patch('kb_python.ref.create_t2c') as create_t2c,\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = True
            temp_dir = tempfile.mkdtemp()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            combined_path = mock.MagicMock()
            create_t2g.return_value = {'t2g': t2g_path}
            create_t2c.side_effect = [{
                't2c': cdna_t2c_path,
                'transcripts': ['A'],
            }, {
                't2c': intron_t2c_path,
                'transcripts': ['B'],
            }]
            concatenate_files.return_value = combined_path
            kallisto_index.return_value = {'index': index_path}
            self.assertEqual({
                't2g': t2g_path,
                'cdna_t2c': cdna_t2c_path,
                'intron_t2c': intron_t2c_path,
            },
                             ref.ref_velocity(
                                 self.cdna_path,
                                 self.intron_path,
                                 self.gtf_path,
                                 index_path,
                                 t2g_path,
                                 cdna_t2c_path,
                                 intron_t2c_path,
                                 temp_dir=temp_dir
                             ))
            create_t2g.assert_called_once_with(
                self.gtf_path, t2g_path, transcripts={'A', 'B'}
            )
            self.assertEqual(create_t2c.call_count, 2)
            create_t2c.assert_has_calls([
                call(self.cdna_path, cdna_t2c_path),
                call(self.intron_path, intron_t2c_path)
            ])
            concatenate_files.assert_not_called()
            kallisto_index.assert_not_called()

    def test_ref_velocity_overwrite(self):
        with mock.patch('kb_python.ref.create_t2g') as create_t2g,\
            mock.patch('kb_python.ref.create_t2c') as create_t2c,\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = True
            temp_dir = tempfile.mkdtemp()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            combined_path = mock.MagicMock()
            create_t2g.return_value = {'t2g': t2g_path}
            create_t2c.side_effect = [{
                't2c': cdna_t2c_path,
                'transcripts': ['A'],
            }, {
                't2c': intron_t2c_path,
                'transcripts': ['B'],
            }]
            concatenate_files.return_value = combined_path
            kallisto_index.return_value = {'index': index_path}
            self.assertEqual({
                't2g': t2g_path,
                'cdna_t2c': cdna_t2c_path,
                'intron_t2c': intron_t2c_path,
                'index': index_path,
            },
                             ref.ref_velocity(
                                 self.cdna_path,
                                 self.intron_path,
                                 self.gtf_path,
                                 index_path,
                                 t2g_path,
                                 cdna_t2c_path,
                                 intron_t2c_path,
                                 temp_dir=temp_dir,
                                 overwrite=True
                             ))
            create_t2g.assert_called_once_with(
                self.gtf_path, t2g_path, transcripts={'A', 'B'}
            )
            self.assertEqual(create_t2c.call_count, 2)
            create_t2c.assert_has_calls([
                call(self.cdna_path, cdna_t2c_path),
                call(self.intron_path, intron_t2c_path)
            ])
            concatenate_files.assert_called_once_with(
                self.cdna_path,
                self.intron_path,
                out_path=os.path.join(temp_dir, COMBINED_FILENAME),
                temp_dir=temp_dir
            )
            kallisto_index.assert_called_once_with(combined_path, index_path)
