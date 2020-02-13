from unittest import mock, TestCase
from unittest.mock import ANY

import numpy as np

import kb_python.report as report
from tests.mixins import TestMixin


class TestReport(TestMixin, TestCase):

    def test_dict_to_table(self):
        with mock.patch('kb_python.report.go') as go:
            d = {'key': 'val', 'k': ['val1', 'val2']}
            self.assertEqual(go.Figure(), report.dict_to_table(d))
            go.Table.assert_called_once_with(
                columnwidth=[3, 7],
                header={'values': ['key', 'value']},
                cells={
                    'values': [['key', 'k', ''], ['val', 'val1', 'val2']],
                    'align': ['right', 'left']
                }
            )

    def test_knee_plot(self):
        with mock.patch('kb_python.report.go') as go:
            n_counts = [1]
            result = report.knee_plot(n_counts)
            go.Scattergl.assert_called_once_with(
                x=[1], y=np.arange(1), mode='lines'
            )
            go.Figure.assert_called_once_with(data=go.Scattergl())
            self.assertEqual(go.Figure(), result)

    def test_genes_detected_plot(self):
        with mock.patch('kb_python.report.go') as go:
            n_counts = [1]
            n_genes = [2]
            result = report.genes_detected_plot(n_counts, n_genes)
            go.Scattergl.assert_called_once_with(
                x=n_counts, y=n_genes, mode='markers'
            )
            go.Figure.assert_called_once_with(data=go.Scattergl())
            self.assertEqual(go.Figure(), result)

    def test_elbow_plot(self):
        with mock.patch('kb_python.report.go') as go:
            pca_variance_ratio = [[1, 2]]
            result = report.elbow_plot(pca_variance_ratio)
            go.Scattergl.assert_called_once_with(
                x=np.arange(1, 2), y=pca_variance_ratio, mode='markers'
            )
            go.Figure.assert_called_once_with(data=go.Scattergl())
            self.assertEqual(go.Figure(), result)

    def test_pca_plot(self):
        with mock.patch('kb_python.report.go') as go:
            pc = np.array([[1, 2]])
            result = report.pca_plot(pc)
            go.Scattergl.assert_called_once_with(x=[1], y=[2], mode='markers')
            go.Figure.assert_called_once_with(data=go.Scattergl())
            self.assertEqual(go.Figure(), result)

    def test_write_report(self):
        with mock.patch('kb_python.report.Template') as Template,\
            mock.patch('kb_python.report.open'):
            stats_path = mock.MagicMock()
            info_path = mock.MagicMock()
            inspect_path = mock.MagicMock()
            matrix_path = mock.MagicMock()
            barcodes_path = mock.MagicMock()
            genes_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            out_path = mock.MagicMock()

            self.assertEqual(
                out_path,
                report.write_report(
                    stats_path, info_path, inspect_path, out_path, matrix_path,
                    barcodes_path, genes_path, t2g_path
                )
            )
            Template.assert_called_once()
            Template().render.assert_called_once_with(
                packages=ANY,
                stats_path=stats_path,
                info_path=info_path,
                inspect_path=inspect_path,
                matrix_path=matrix_path,
                barcodes_path=barcodes_path,
                genes_path=genes_path,
                t2g_path=t2g_path
            )

    def test_execute_report(self):
        with mock.patch('kb_python.report.nbformat') as nbformat,\
            mock.patch('kb_python.report.ExecutePreprocessor') as ExecutePreprocessor,\
            mock.patch('kb_python.report.HTMLExporter') as HTMLExporter,\
            mock.patch('kb_python.report.open'):
            execute_path = mock.MagicMock()
            nb_path = mock.MagicMock()
            html_path = mock.MagicMock()

            HTMLExporter(
            ).from_notebook_node.return_value = ('html', 'resources')
            self.assertEqual(
                (nb_path, html_path),
                report.execute_report(execute_path, nb_path, html_path)
            )
            nbformat.read.assert_called_once()
            ExecutePreprocessor.assert_called_once()
            ExecutePreprocessor().preprocess.assert_called_once_with(
                nbformat.read()
            )
            nbformat.write.assert_called_once_with(nbformat.read(), ANY)
            HTMLExporter().from_notebook_node.assert_called_once_with(
                nbformat.read()
            )

    def test_render_report(self):
        with mock.patch('kb_python.report.write_report') as write_report,\
            mock.patch('kb_python.report.execute_report') as execute_report,\
            mock.patch('kb_python.report.get_temporary_filename') as get_temporary_filename:
            stats_path = mock.MagicMock()
            info_path = mock.MagicMock()
            inspect_path = mock.MagicMock()
            matrix_path = mock.MagicMock()
            barcodes_path = mock.MagicMock()
            genes_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            nb_path = mock.MagicMock()
            html_path = mock.MagicMock()
            temp_dir = mock.MagicMock()

            self.assertEqual({
                'report_notebook': nb_path,
                'report_html': html_path
            },
                             report.render_report(
                                 stats_path, info_path, inspect_path, nb_path,
                                 html_path, matrix_path, barcodes_path,
                                 genes_path, t2g_path, temp_dir
                             ))
            get_temporary_filename.assert_called_once_with(temp_dir)
            write_report.assert_called_once_with(
                stats_path, info_path, inspect_path, get_temporary_filename(),
                matrix_path, barcodes_path, genes_path, t2g_path
            )
            execute_report.assert_called_once_with(
                write_report(), nb_path, html_path
            )
