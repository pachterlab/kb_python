from unittest import mock, TestCase
from unittest.mock import ANY, call

import numpy as np

import kb_python.report as report
from tests.mixins import TestMixin


class TestReport(TestMixin, TestCase):

    def test_dict_to_table(self):
        with mock.patch('kb_python.report.go') as go,\
             mock.patch('kb_python.report.py') as py:
            d = {'key': 'val', 'k': ['val1', 'val2']}
            self.assertEqual(py.to_html(), report.dict_to_table(d))
            go.Table.assert_called_once_with(
                columnwidth=[3, 7],
                header={
                    'line_color': 'white',
                    'fill_color': 'white'
                },
                cells={
                    'values': [['key', 'k', ''], ['val', 'val1', 'val2']],
                    'align': ['right', 'left']
                }
            )

    def test_format_card(self):
        title = 'title'
        subtitle = 'subtitle'
        body = 'body'
        self.assertEqual({
            'title': 'title',
            'subtitle': 'subtitle',
            'body': 'body'
        }, report.format_card(title, subtitle, body))

    def test_kb_info(self):
        with mock.patch('kb_python.report.format_card') as format_card,\
            mock.patch('kb_python.report.dict_to_table') as dict_to_table:
            stats = {'some': 'stats'}
            result = report.kb_info(stats)
            dict_to_table.assert_called_once_with(stats)
            format_card.assert_called_once_with(ANY, ANY, dict_to_table())
            self.assertEqual(format_card(), result)

    def test_kallisto_info(self):
        with mock.patch('kb_python.report.format_card') as format_card,\
            mock.patch('kb_python.report.dict_to_table') as dict_to_table:
            info = {'some': 'stats'}
            result = report.kallisto_info(info)
            dict_to_table.assert_called_once_with(info)
            format_card.assert_called_once_with(ANY, ANY, dict_to_table())
            self.assertEqual(format_card(), result)

    def test_bus_info(self):
        with mock.patch('kb_python.report.format_card') as format_card,\
            mock.patch('kb_python.report.dict_to_table') as dict_to_table:
            inspect = {'some': 'stats'}
            result = report.bus_info(inspect)
            dict_to_table.assert_called_once_with(inspect)
            format_card.assert_called_once_with(ANY, ANY, dict_to_table())
            self.assertEqual(format_card(), result)

    def test_gene_info(self):
        with mock.patch('kb_python.report.format_card') as format_card,\
            mock.patch('kb_python.report.dict_to_table') as dict_to_table:
            n_counts = [1]
            n_genes = [2]
            result = report.gene_info(n_counts, n_genes)
            dict_to_table.assert_called_once()
            format_card.assert_called_once_with(ANY, ANY, dict_to_table())
            self.assertEqual(format_card(), result)

    def test_knee_plot(self):
        with mock.patch('kb_python.report.format_card') as format_card,\
            mock.patch('kb_python.report.go') as go,\
            mock.patch('kb_python.report.py') as py:
            n_counts = [1]
            result = report.knee_plot(n_counts)
            go.Scattergl.assert_called_once_with(
                x=[1], y=np.arange(1), mode='lines'
            )
            go.Figure.assert_called_once_with(data=go.Scattergl())
            py.to_html.assert_called_once_with(
                go.Figure(), include_plotlyjs=False, full_html=False
            )
            format_card.assert_called_once_with(ANY, ANY, py.to_html())
            self.assertEqual(format_card(), result)

    def test_genes_detected_plot(self):
        with mock.patch('kb_python.report.format_card') as format_card,\
            mock.patch('kb_python.report.go') as go,\
            mock.patch('kb_python.report.py') as py:
            n_counts = [1]
            n_genes = [2]
            result = report.genes_detected_plot(n_counts, n_genes)
            go.Scattergl.assert_called_once_with(
                x=n_counts, y=n_genes, mode='markers'
            )
            go.Figure.assert_called_once_with(data=go.Scattergl())
            py.to_html.assert_called_once_with(
                go.Figure(), include_plotlyjs=False, full_html=False
            )
            format_card.assert_called_once_with(ANY, ANY, py.to_html())
            self.assertEqual(format_card(), result)

    def test_elbow_plot(self):
        with mock.patch('kb_python.report.format_card') as format_card,\
            mock.patch('kb_python.report.go') as go,\
            mock.patch('kb_python.report.py') as py:
            pca_variance_ratio = [[1, 2]]
            result = report.elbow_plot(pca_variance_ratio)
            go.Scattergl.assert_called_once_with(
                x=np.arange(1, 2), y=pca_variance_ratio, mode='markers'
            )
            go.Figure.assert_called_once_with(data=go.Scattergl())
            py.to_html.assert_called_once_with(
                go.Figure(), include_plotlyjs=False, full_html=False
            )
            format_card.assert_called_once_with(ANY, ANY, py.to_html())
            self.assertEqual(format_card(), result)

    def test_pca_plot(self):
        with mock.patch('kb_python.report.format_card') as format_card,\
            mock.patch('kb_python.report.go') as go,\
            mock.patch('kb_python.report.py') as py:
            pc = np.array([[1, 2]])
            result = report.pca_plot(pc)
            go.Scattergl.assert_called_once_with(x=[1], y=[2], mode='markers')
            go.Figure.assert_called_once_with(data=go.Scattergl())
            py.to_html.assert_called_once_with(
                go.Figure(), include_plotlyjs=False, full_html=False
            )
            format_card.assert_called_once_with(ANY, ANY, py.to_html())
            self.assertEqual(format_card(), result)

    def test_render_report(self):
        with mock.patch('kb_python.report.kb_info') as kb_info,\
            mock.patch('kb_python.report.kallisto_info') as kallisto_info,\
            mock.patch('kb_python.report.bus_info') as bus_info,\
            mock.patch('kb_python.report.gene_info') as gene_info,\
            mock.patch('kb_python.report.knee_plot') as knee_plot,\
            mock.patch('kb_python.report.genes_detected_plot') as genes_detected_plot,\
            mock.patch('kb_python.report.elbow_plot') as elbow_plot,\
            mock.patch('kb_python.report.pca_plot') as pca_plot,\
            mock.patch('kb_python.report.sc') as sc,\
            mock.patch('kb_python.report.Template') as Template,\
            mock.patch('kb_python.report.open'),\
            mock.patch('kb_python.report.json') as json,\
            mock.patch('kb_python.report.PCA') as PCA:
            adata = mock.MagicMock()
            stats = mock.MagicMock()
            info_path = mock.MagicMock()
            inspect_path = mock.MagicMock()
            out_path = mock.MagicMock()

            self.assertEqual({'report': out_path},
                             report.render_report(
                                 stats,
                                 info_path,
                                 inspect_path,
                                 out_path,
                                 adata=adata
                             ))
            kb_info.assert_called_once_with(stats)
            kallisto_info.assert_called_once_with(json.load())
            bus_info.assert_called_once_with(json.load())

            self.assertEqual(2, sc.pp.filter_cells.call_count)
            sc.pp.filter_cells.assert_has_calls([
                call(adata, min_genes=0),
                call(adata, min_counts=0)
            ])
            sc.pp.normalize_total.assert_called_once_with(adata, target_sum=1e4)
            sc.pp.log1p.assert_called_once_with(adata)
            PCA.assert_called_once_with(n_components=10)
            PCA().fit_transform.assert_called_once_with(adata.X.todense())
            gene_info.assert_called_once_with(
                adata.obs['n_counts'], adata.obs['n_genes']
            )
            knee_plot.assert_called_once_with(adata.obs['n_counts'])
            genes_detected_plot.assert_called_once_with(
                adata.obs['n_counts'], adata.obs['n_genes']
            )
            elbow_plot.assert_called_once_with(PCA().explained_variance_ratio_)
            pca_plot.assert_called_once_with(PCA().fit_transform())
            Template.assert_called_once()
            Template().render.assert_called_once()

    def test_render_report_no_adata(self):
        with mock.patch('kb_python.report.kb_info') as kb_info,\
            mock.patch('kb_python.report.kallisto_info') as kallisto_info,\
            mock.patch('kb_python.report.bus_info') as bus_info,\
            mock.patch('kb_python.report.gene_info') as gene_info,\
            mock.patch('kb_python.report.knee_plot') as knee_plot,\
            mock.patch('kb_python.report.genes_detected_plot') as genes_detected_plot,\
            mock.patch('kb_python.report.elbow_plot') as elbow_plot,\
            mock.patch('kb_python.report.pca_plot') as pca_plot,\
            mock.patch('kb_python.report.sc') as sc,\
            mock.patch('kb_python.report.Template') as Template,\
            mock.patch('kb_python.report.open'),\
            mock.patch('kb_python.report.json') as json,\
            mock.patch('kb_python.report.PCA') as PCA:
            stats = mock.MagicMock()
            info_path = mock.MagicMock()
            inspect_path = mock.MagicMock()
            out_path = mock.MagicMock()

            self.assertEqual({'report': out_path},
                             report.render_report(
                                 stats,
                                 info_path,
                                 inspect_path,
                                 out_path,
                                 adata=None
                             ))
            kb_info.assert_called_once_with(stats)
            kallisto_info.assert_called_once_with(json.load())
            bus_info.assert_called_once_with(json.load())

            sc.pp.filter_cells.assert_not_called()
            sc.pp.normalize_total.assert_not_called()
            sc.pp.log1p.assert_not_called()
            PCA.assert_not_called()
            PCA().fit_transform.assert_not_called()
            gene_info.assert_not_called()
            knee_plot.assert_not_called()
            genes_detected_plot.assert_not_called()
            elbow_plot.assert_not_called()
            pca_plot.assert_not_called()
            Template.assert_called_once()
            Template().render.assert_called_once()
