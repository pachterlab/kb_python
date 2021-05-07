import os

import nbformat
import numpy as np
import plotly.graph_objects as go
from nbconvert import HTMLExporter
from nbconvert.preprocessors import ExecutePreprocessor
from jinja2 import Template

from . import __version__
from .config import PACKAGE_PATH
from .dry import dryable
from .utils import get_temporary_filename

REPORT_DIR = os.path.join(PACKAGE_PATH, 'report')
BASIC_TEMPLATE_PATH = os.path.join(REPORT_DIR, 'report_basic.ipynb')
MATRIX_TEMPLATE_PATH = os.path.join(REPORT_DIR, 'report_matrix.ipynb')

MARGIN = go.layout.Margin(t=5, r=5, b=0, l=5)  # noqa: E741


def dict_to_table(d, column_ratio=[3, 7], column_align=['right', 'left']):
    """Convert a dictionary to a Plot.ly table of key-value pairs.

    :param d: dictionary to convert
    :type d: dict
    :param column_ratio: relative column widths, represented as a ratio,
                         defaults to `[3, 7]`
    :type column_ratio: list, optional
    :param column_align: column text alignments, defaults to `['right', 'left']`
    :type column_align: list, optional

    :return: figure
    :rtype: plotly.graph_objs.Figure
    """
    keys = []
    values = []
    for key, value in d.items():
        if isinstance(value, list):
            keys.append(key)
            values.append(value[0])
            for val in value[1:]:
                keys.append('')
                values.append(val)
        else:
            keys.append(key)
            values.append(value)

    table = go.Table(
        columnwidth=column_ratio,
        header={'values': ['key', 'value']},
        cells={
            'values': [keys, values],
            'align': column_align
        }
    )
    figure = go.Figure(data=table)
    figure.update_layout(
        margin=MARGIN,
        xaxis_automargin=True,
        yaxis_automargin=True,
        autosize=True
    )
    return figure


def knee_plot(n_counts):
    """Generate knee plot card.

    :param n_counts: list of UMI counts
    :type n_counts: list

    :return: figure
    :rtype: plotly.graph_objs.Figure
    """
    knee = np.sort(n_counts)[::-1]
    scatter = go.Scattergl(x=knee, y=np.arange(len(knee)), mode='lines')
    figure = go.Figure(data=scatter)
    figure.update_layout(
        margin=MARGIN,
        xaxis_title='UMI counts',
        yaxis_title='Number of barcodes',
        xaxis_type='log',
        yaxis_type='log',
        xaxis_automargin=True,
        yaxis_automargin=True,
        autosize=True
    )
    return figure


def genes_detected_plot(n_counts, n_genes):
    """Generate genes detected plot card.

    :param n_counts: list of UMI counts
    :type n_counts: list
    :param n_genes: list of gene counts
    :type n_genes: list

    :return: figure
    :rtype: plotly.graph_objs.Figure
    """
    scatter = go.Scattergl(x=n_counts, y=n_genes, mode='markers')
    figure = go.Figure(data=scatter)
    figure.update_layout(
        margin=MARGIN,
        xaxis_title='UMI counts',
        yaxis_title='Genes detected',
        xaxis_type='log',
        yaxis_type='log',
        xaxis_automargin=True,
        yaxis_automargin=True,
        autosize=True
    )
    return figure


def elbow_plot(pca_variance_ratio):
    """Generate elbow plot card.

    :param pca_variance_ratio: list PCA variance ratios
    :type pca_variance_ratio: list

    :return: figure
    :rtype: plotly.graph_objs.Figure
    """
    scatter = go.Scattergl(
        x=np.arange(1,
                    len(pca_variance_ratio) + 1),
        y=pca_variance_ratio,
        mode='markers'
    )
    figure = go.Figure(data=scatter)
    figure.update_layout(
        margin=MARGIN,
        xaxis_title='PC',
        yaxis_title='Explained variance ratio',
        xaxis_automargin=True,
        yaxis_automargin=True,
        autosize=True
    )
    return figure


def pca_plot(pc):
    """Generate PCA plot card.

    :param pc: embeddings
    :type pc: list

    :return: figure
    :rtype: plotly.graph_objs.Figure
    """
    scatter = go.Scattergl(x=pc[:, 0], y=pc[:, 1], mode='markers')
    figure = go.Figure(data=scatter)
    figure.update_layout(
        margin=MARGIN,
        xaxis_title='PC 1',
        yaxis_title='PC 2',
        xaxis_automargin=True,
        yaxis_automargin=True,
        autosize=True
    )
    return figure


def write_report(
    stats_path,
    info_path,
    inspect_path,
    out_path,
    matrix_path=None,
    barcodes_path=None,
    genes_path=None,
    t2g_path=None
):
    """Render the Jupyter notebook report with Jinja2.

    :param stats_path: path to kb stats JSON
    :type stats_path: str
    :param info_path: path to run_info.json
    :type info_path: str
    :param inspect_path: path to inspect.json
    :type inspect_path: str
    :param out_path: path to Jupyter notebook to generate
    :type out_path: str
    :param matrix_path: path to matrix
    :type matrix_path: str
    :param barcodes_path: list of paths to barcodes.txt
    :type barcodes_path: str
    :param genes_path: path to genes.txt, defaults to `None`
    :type genes_path: str, optional
    :param t2g_path: path to transcript-to-gene mapping
    :type t2g_path: str

    :return: path to notebook generated
    :rtype: str
    """
    template_path = MATRIX_TEMPLATE_PATH if all(
        p is not None
        for p in [matrix_path, barcodes_path, genes_path, t2g_path]
    ) else BASIC_TEMPLATE_PATH
    with open(template_path, 'r') as f, open(out_path, 'w') as out:
        template = Template(f.read())
        out.write(
            template.render(
                packages=f'#!pip install kb-python>={__version__}',
                stats_path=stats_path,
                info_path=info_path,
                inspect_path=inspect_path,
                matrix_path=matrix_path,
                barcodes_path=barcodes_path,
                genes_path=genes_path,
                t2g_path=t2g_path
            )
        )

    return out_path


def execute_report(execute_path, nb_path, html_path):
    """Execute the report and write the results as a Jupyter notebook and HTML.

    :param execute_path: path to Jupyter notebook to execute
    :type execute_path: str
    :param nb_path: path to Jupyter notebook to generate
    :type nb_path: str
    :param html_path: path to HTML to generate
    :type html_path: str

    :return: tuple containing executed notebook and HTML
    :rtype: tuple
    """
    with open(execute_path, 'r') as f:
        nb = nbformat.read(f, as_version=4)

    ep = ExecutePreprocessor(timeout=600)
    ep.preprocess(nb)

    with open(nb_path, 'w') as f:
        nbformat.write(nb, f)

    with open(html_path, 'w') as f:
        html_exporter = HTMLExporter()
        html, resources = html_exporter.from_notebook_node(nb)
        f.write(html)

    return nb_path, html_path


@dryable(lambda *args, **kwargs: {})
def render_report(
    stats_path,
    info_path,
    inspect_path,
    nb_path,
    html_path,
    matrix_path=None,
    barcodes_path=None,
    genes_path=None,
    t2g_path=None,
    temp_dir='tmp'
):
    """Render and execute the report.

    :param stats_path: path to kb stats JSON
    :type stats_path: str
    :param info_path: path to run_info.json
    :type info_path: str
    :param inspect_path: path to inspect.json
    :type inspect_path: str
    :param nb_path: path to Jupyter notebook to generate
    :type nb_path: str
    :param html_path: path to HTML to generate
    :type html_path: str
    :param matrix_path: path to matrix
    :type matrix_path: str
    :param barcodes_path: list of paths to barcodes.txt
    :type barcodes_path: str
    :param genes_path: path to genes.txt, defaults to `None`
    :type genes_path: str, optional
    :param t2g_path: path to transcript-to-gene mapping
    :type t2g_path: str
    :param temp_dir: path to temporary directory, defaults to `tmp`
    :type temp_dir: str, optional

    :return: dictionary containing notebook and HTML paths
    :rtype: dict
    """
    temp_path = write_report(
        stats_path, info_path, inspect_path, get_temporary_filename(temp_dir),
        matrix_path, barcodes_path, genes_path, t2g_path
    )
    execute_report(temp_path, nb_path, html_path)

    return {'report_notebook': nb_path, 'report_html': html_path}
