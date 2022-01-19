import os
import tempfile
from unittest import mock, TestCase
from unittest.mock import call

import kb_python.compile as compile
from tests.mixins import TestMixin


class TestCompile(TestMixin, TestCase):

    def test_get_latest_github_release_tag(self):
        with mock.patch('kb_python.compile.requests') as requests:
            response = mock.MagicMock()
            requests.get.return_value = response
            response.json.return_value = [{
                'tag_name': 'tag1'
            }, {
                'tag_name': 'tag2'
            }]
            self.assertEqual(
                'tag1', compile.get_latest_github_release_tag('url')
            )
            requests.get.assert_called_once_with('url')
            response.raise_for_status.assert_called_once_with()

    def test_get_filename_from_url_disposition(self):
        with mock.patch('kb_python.compile.requests') as requests:
            response = mock.MagicMock()
            requests.get.return_value = response
            response.headers = {
                'content-disposition':
                    'attachment; filename=pachterlab-kallisto-v0.46.2-0-gae81a86.tar.gz'
            }
            self.assertEqual(
                'pachterlab-kallisto-v0.46.2-0-gae81a86.tar.gz',
                compile.get_filename_from_url('url')
            )
            requests.get.assert_called_once_with('url')
            response.raise_for_status.assert_called_once_with()

    def test_get_filename_from_url_without_disposition(self):
        with mock.patch('kb_python.compile.requests') as requests:
            response = mock.MagicMock()
            requests.get.return_value = response
            response.headers = {}
            self.assertEqual(
                'file.tar.gz',
                compile.
                get_filename_from_url('https://test.com/file.tar.gz?test=test')
            )
            requests.get.assert_called_once_with(
                'https://test.com/file.tar.gz?test=test'
            )
            response.raise_for_status.assert_called_once_with()

    def test_get_kallisto_url(self):
        with mock.patch('kb_python.compile.get_latest_github_release_tag'
                        ) as get_latest_github_release_tag:
            get_latest_github_release_tag.return_value = 'test'
            self.assertEqual(
                'https://api.github.com/repos/pachterlab/kallisto/tarball/test',
                compile.get_kallisto_url()
            )
            get_latest_github_release_tag.assert_called_once_with(
                'https://api.github.com/repos/pachterlab/kallisto/releases'
            )

    def test_get_kallisto_url_ref(self):
        with mock.patch('kb_python.compile.get_latest_github_release_tag'
                        ) as get_latest_github_release_tag:
            self.assertEqual(
                'https://api.github.com/repos/pachterlab/kallisto/tarball/ref',
                compile.get_kallisto_url('ref')
            )
            get_latest_github_release_tag.assert_not_called()

    def test_get_bustools_url(self):
        with mock.patch('kb_python.compile.get_latest_github_release_tag'
                        ) as get_latest_github_release_tag:
            get_latest_github_release_tag.return_value = 'test'
            self.assertEqual(
                'https://api.github.com/repos/BUStools/bustools/tarball/test',
                compile.get_bustools_url()
            )
            get_latest_github_release_tag.assert_called_once_with(
                'https://api.github.com/repos/BUStools/bustools/releases'
            )

    def test_get_bustools_url_ref(self):
        with mock.patch('kb_python.compile.get_latest_github_release_tag'
                        ) as get_latest_github_release_tag:
            self.assertEqual(
                'https://api.github.com/repos/BUStools/bustools/tarball/ref',
                compile.get_bustools_url('ref')
            )
            get_latest_github_release_tag.assert_not_called()

    def test_find_git_root(self):
        d1 = tempfile.mkdtemp(dir=self.temp_dir)
        d2 = tempfile.mkdtemp(dir=d1)
        tempfile.mkdtemp(dir=self.temp_dir)
        with open(os.path.join(d2, '.gitignore'), 'w') as f:
            f.write('test')

        self.assertEqual(d2, compile.find_git_root(self.temp_dir))

    def test_compile_kallisto(self):
        with mock.patch('kb_python.compile.run_executable') as run_executable:
            source_dir = tempfile.mkdtemp(dir=self.temp_dir)
            os.makedirs(os.path.join(source_dir, 'ext', 'htslib'))
            os.makedirs(os.path.join(source_dir, 'build', 'src'))
            with open(os.path.join(source_dir, 'build', 'src', 'kallisto'),
                      'w') as f:
                f.write('test')
            with open(os.path.join(source_dir, 'license.txt'), 'w') as f:
                f.write('test')
            binary_path = os.path.join(self.temp_dir, 'kallisto')
            self.assertEqual(
                binary_path, compile.compile_kallisto(source_dir, binary_path)
            )
            run_executable.assert_has_calls([
                call(['autoheader']),
                call(['autoconf']),
                call(['cmake', '..']),
                call(['make'])
            ])

    def test_compile_bustools(self):
        with mock.patch('kb_python.compile.run_executable') as run_executable:
            source_dir = tempfile.mkdtemp(dir=self.temp_dir)
            os.makedirs(os.path.join(source_dir, 'build', 'src'))
            with open(os.path.join(source_dir, 'build', 'src', 'bustools'),
                      'w') as f:
                f.write('test')
            with open(os.path.join(source_dir, 'LICENSE'), 'w') as f:
                f.write('test')
            binary_path = os.path.join(self.temp_dir, 'bustools')
            self.assertEqual(
                binary_path, compile.compile_bustools(source_dir, binary_path)
            )
            run_executable.assert_has_calls([
                call(['cmake', '..']), call(['make'])
            ])

    def test_compile_full_kallisto(self):
        with mock.patch('kb_python.compile.get_kallisto_url', return_value='kallisto_url'),\
            mock.patch('kb_python.compile.get_bustools_url', return_value='bustools_url'),\
            mock.patch('kb_python.compile.download_file') as download_file,\
            mock.patch('kb_python.compile.get_filename_from_url', return_value='filename'),\
            mock.patch('kb_python.compile.shutil.unpack_archive') as unpack_archive,\
            mock.patch('kb_python.compile.find_git_root') as find_git_root,\
            mock.patch('kb_python.compile.compile_kallisto') as compile_kallisto,\
            mock.patch('kb_python.compile.compile_bustools') as compile_bustools:
            out_dir = self.temp_dir
            self.assertEqual({
                'kallisto': compile_kallisto.return_value
            }, compile.compile('kallisto', out_dir, temp_dir=self.temp_dir))
            download_file.assert_called_once_with(
                'kallisto_url', os.path.join(self.temp_dir, 'filename')
            )
            unpack_archive.assert_called_once_with(
                download_file.return_value, mock.ANY
            )
            find_git_root.assert_called_once_with(mock.ANY)
            compile_kallisto.assert_called_once_with(
                find_git_root.return_value,
                os.path.join(out_dir, 'kallisto'),
                cmake_arguments=None
            )
            compile_bustools.assert_not_called()

    def test_compile_full_bustools(self):
        with mock.patch('kb_python.compile.get_kallisto_url', return_value='kallisto_url'),\
            mock.patch('kb_python.compile.get_bustools_url', return_value='bustools_url'),\
            mock.patch('kb_python.compile.download_file') as download_file,\
            mock.patch('kb_python.compile.get_filename_from_url', return_value='filename'),\
            mock.patch('kb_python.compile.shutil.unpack_archive') as unpack_archive,\
            mock.patch('kb_python.compile.find_git_root') as find_git_root,\
            mock.patch('kb_python.compile.compile_kallisto') as compile_kallisto,\
            mock.patch('kb_python.compile.compile_bustools') as compile_bustools:
            out_dir = self.temp_dir
            self.assertEqual({
                'bustools': compile_bustools.return_value
            }, compile.compile('bustools', out_dir, temp_dir=self.temp_dir))
            download_file.assert_called_once_with(
                'bustools_url', os.path.join(self.temp_dir, 'filename')
            )
            unpack_archive.assert_called_once_with(
                download_file.return_value, mock.ANY
            )
            find_git_root.assert_called_once_with(mock.ANY)
            compile_bustools.assert_called_once_with(
                find_git_root.return_value,
                os.path.join(out_dir, 'bustools'),
                cmake_arguments=None
            )
            compile_kallisto.assert_not_called()

    def test_compile_full_all(self):
        with mock.patch('kb_python.compile.get_kallisto_url', return_value='kallisto_url'),\
            mock.patch('kb_python.compile.get_bustools_url', return_value='bustools_url'),\
            mock.patch('kb_python.compile.download_file', side_effect=['download1', 'download2']) as download_file,\
            mock.patch('kb_python.compile.get_filename_from_url', side_effect=['filename1', 'filename2']),\
            mock.patch('kb_python.compile.shutil.unpack_archive') as unpack_archive,\
            mock.patch('kb_python.compile.find_git_root', side_effect=['root1', 'root2']) as find_git_root,\
            mock.patch('kb_python.compile.compile_kallisto') as compile_kallisto,\
            mock.patch('kb_python.compile.compile_bustools') as compile_bustools:
            out_dir = self.temp_dir
            self.assertEqual({
                'kallisto': compile_kallisto.return_value,
                'bustools': compile_bustools.return_value,
            }, compile.compile('all', out_dir, temp_dir=self.temp_dir))
            download_file.assert_has_calls([
                call('kallisto_url', os.path.join(self.temp_dir, 'filename1')),
                call('bustools_url', os.path.join(self.temp_dir, 'filename2'))
            ])
            unpack_archive.assert_has_calls([
                call('download1', mock.ANY),
                call('download2', mock.ANY)
            ])
            find_git_root.assert_has_calls([call(mock.ANY), call(mock.ANY)])
            compile_kallisto.assert_called_once_with(
                'root1',
                os.path.join(out_dir, 'kallisto'),
                cmake_arguments=None
            )
            compile_bustools.assert_called_once_with(
                'root2',
                os.path.join(out_dir, 'bustools'),
                cmake_arguments=None
            )

    def test_compile_full_with_url(self):
        with mock.patch('kb_python.compile.get_kallisto_url', return_value='kallisto_url'),\
            mock.patch('kb_python.compile.get_bustools_url', return_value='bustools_url'),\
            mock.patch('kb_python.compile.download_file') as download_file,\
            mock.patch('kb_python.compile.get_filename_from_url', return_value='filename'),\
            mock.patch('kb_python.compile.shutil.unpack_archive') as unpack_archive,\
            mock.patch('kb_python.compile.find_git_root') as find_git_root,\
            mock.patch('kb_python.compile.compile_kallisto') as compile_kallisto,\
            mock.patch('kb_python.compile.compile_bustools') as compile_bustools:
            out_dir = self.temp_dir
            self.assertEqual({'kallisto': compile_kallisto.return_value},
                             compile.compile(
                                 'kallisto',
                                 out_dir,
                                 url='testurl',
                                 temp_dir=self.temp_dir
                             ))
            download_file.assert_called_once_with(
                'testurl', os.path.join(self.temp_dir, 'filename')
            )
            unpack_archive.assert_called_once_with(
                download_file.return_value, mock.ANY
            )
            find_git_root.assert_called_once_with(mock.ANY)
            compile_kallisto.assert_called_once_with(
                find_git_root.return_value,
                os.path.join(out_dir, 'kallisto'),
                cmake_arguments=None
            )
            compile_bustools.assert_not_called()

    def test_compile_full_with_ref(self):
        with mock.patch('kb_python.compile.get_kallisto_url', return_value='kallisto_url') as get_kallisto_url,\
            mock.patch('kb_python.compile.get_bustools_url', return_value='bustools_url'),\
            mock.patch('kb_python.compile.download_file') as download_file,\
            mock.patch('kb_python.compile.get_filename_from_url', return_value='filename'),\
            mock.patch('kb_python.compile.shutil.unpack_archive') as unpack_archive,\
            mock.patch('kb_python.compile.find_git_root') as find_git_root,\
            mock.patch('kb_python.compile.compile_kallisto') as compile_kallisto,\
            mock.patch('kb_python.compile.compile_bustools') as compile_bustools:
            out_dir = self.temp_dir
            self.assertEqual({'kallisto': compile_kallisto.return_value},
                             compile.compile(
                                 'kallisto',
                                 out_dir,
                                 ref='ref',
                                 temp_dir=self.temp_dir
                             ))
            get_kallisto_url.assert_called_once_with('ref')
            download_file.assert_called_once_with(
                'kallisto_url', os.path.join(self.temp_dir, 'filename')
            )
            unpack_archive.assert_called_once_with(
                download_file.return_value, mock.ANY
            )
            find_git_root.assert_called_once_with(mock.ANY)
            compile_kallisto.assert_called_once_with(
                find_git_root.return_value,
                os.path.join(out_dir, 'kallisto'),
                cmake_arguments=None
            )
            compile_bustools.assert_not_called()
