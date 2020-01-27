import json
import logging
import os

import numpy as np
import plotly.graph_objects as go
import plotly.io as py
import scanpy as sc
from jinja2 import Template
from sklearn.decomposition import PCA

from .config import PACKAGE_PATH

logger = logging.getLogger(__name__)

REPORT_DIR = os.path.join(PACKAGE_PATH, 'report')
REPORT_TEMPLATE_PATH = os.path.join(REPORT_DIR, 'base.html')
REPORT_CSS_DIR = os.path.join(REPORT_DIR, 'css')
REPORT_JS_DIR = os.path.join(REPORT_DIR, 'js')
MARGIN = go.layout.Margin(t=5, r=5, b=5, l=5)  # noqa: E741


def dict_to_table(d, column_ratio=[3, 7], column_align=['right', 'left']):
    """Convert a dictionary to a Plot.ly table of key-value pairs.

    :param d: dictionary to convert
    :type d: dict
    :param column_ratio: relative column widths, represented as a ratio,
                         defaults to `[3, 7]`
    :type column_ratio: list, optional
    :param column_align: column text alignments, defaults to `['right', 'left']`
    :type column_align: list, optional

    :return: HTML string
    :rtype: str
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
        header={
            'line_color': 'white',
            'fill_color': 'white'
        },
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
    return py.to_html(figure, include_plotlyjs=False, full_html=False)


def format_card(title, subtitle, body):
    """Format card title, subtitle, and body texts to a format usable by the
    Jinja2 template.

    :param title: card title
    :type title: str
    :param subtitle: card subtitle
    :type subtitle: str
    :param body: card body
    :type body: str

    :return: dictionary representation of the card
    :rtype: dict
    """
    return {'title': title, 'subtitle': subtitle, 'body': body}


def kb_info(stats):
    """Generate the kb run info card.

    :param stats: dictionary of run statistics to display
    :type stats: dict

    :return: dictionary representation of the card
    :rtype: dict
    """
    return format_card(
        'kb run info', 'Overall run statistics', dict_to_table(stats)
    )


def kallisto_info(info):
    """Generate the kallisto info card.

    :param info: dictionary of kallisto run info, read from run_info.json
    :type info: dict

    :return: dictionary representation of the card
    :rtype: dict
    """
    return format_card(
        'kallisto run info', 'From kallisto log (run_info.json)',
        dict_to_table(info)
    )


def bus_info(inspect):
    """Generate the bus info card.

    :param info: dictionary of bus info, read from inspect.json
    :type info: dict

    :return: dictionary representation of the card
    :rtype: dict
    """
    return format_card(
        'Bus file info',
        'From <code>bustools inspect</code> command (inspect.json)',
        dict_to_table(inspect)
    )


def gene_info(n_counts, n_genes):
    """Generate the gene info card.

    :param n_counts: list of UMI counts
    :type n_counts: list
    :param n_genes: list of gene counts
    :type n_genes: list

    :return: dictionary representation of the card
    :rtype: dict
    """
    d = {
        'Median UMIs per gene': np.median(n_counts),
        'Mean UMIs per gene': np.mean(n_counts),
        'Median genes per cell': np.median(n_genes),
        'Mean genes per cell': np.mean(n_genes),
    }
    return format_card(
        'Gene matrix info', 'Statistics calculated from gene matrix',
        dict_to_table(d)
    )


def knee_plot(n_counts):
    """Generate knee plot card.

    :param n_counts: list of UMI counts
    :type n_counts: list

    :return: dictionary representation of the card
    :rtype: dict
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
    return format_card(
        'Knee plot', (
            'For a given UMI count (x-axis), the number of cells that '
            'contain at least that many UMI counts (y-axis).'
        ), py.to_html(figure, include_plotlyjs=False, full_html=False)
    )


def genes_detected_plot(n_counts, n_genes):
    """Generate genes detected plot card.

    :param n_counts: list of UMI counts
    :type n_counts: list
    :param n_genes: list of gene counts
    :type n_genes: list

    :return: dictionary representation of the card
    :rtype: dict
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
    return format_card(
        'Genes detected',
        'Number of genes detected as a function of distinct UMI counts per cell.',
        py.to_html(figure, include_plotlyjs=False, full_html=False)
    )


def elbow_plot(pca_variance_ratio):
    """Generate elbow plot card.

    :param pca_variance_ratio: list PCA variance ratios
    :type pca_variance_ratio: list

    :return: dictionary representation of the card
    :rtype: dict
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
    return format_card(
        'Elbow plot',
        'Ratio of variance in data explained by first ten principal components.',
        py.to_html(figure, include_plotlyjs=False, full_html=False)
    )


def pca_plot(pc):
    """Generate PCA plot card.

    :param pc: embeddings
    :type pc: list

    :return: dictionary representation of the card
    :rtype: dict
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
    return format_card(
        'Principal component analysis', 'First two principal components.',
        py.to_html(figure, include_plotlyjs=False, full_html=False)
    )


def render_report(adata, stats, info_path, inspect_path, out_path, tcc=False):
    """Render the HTML report with Jinja2.

    :param adata: anndata to generate report for
    :type adata: anndata
    :param stats: overall run stats
    :type stats: dict
    :param info_path: path to run_info.json
    :type info_path: str
    :param inspect_path: path to inspect.json
    :type inspect_path: str
    :param out_path: path to report to generate
    :type out_path: str
    :param tcc: whether the anndatas are TCC matrices, defaults to `False`
    :type tcc: bool, optional

    :return: dictionary containing path to generated report
    :rtype: dict
    """
    # Construct list of cards to render.
    # Each sub-list is a row.
    cards = []
    cards.append([kb_info(stats)])
    with open(info_path, 'r') as inf, open(inspect_path, 'r') as ins:
        info = json.load(inf)
        inspect = json.load(ins)
    cards.append([kallisto_info(info), bus_info(inspect)])

    if not tcc:
        sc.pp.filter_cells(adata, min_genes=0)
        sc.pp.filter_cells(adata, min_counts=0)
        n_counts = adata.obs['n_counts']
        n_genes = adata.obs['n_genes']
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        pca = PCA(n_components=10)
        pc = pca.fit_transform(adata.X.todense())

        cards.append([gene_info(n_counts, n_genes)])
        cards.append([
            knee_plot(n_counts),
            genes_detected_plot(n_counts, n_genes)
        ])
        cards.append([elbow_plot(pca.explained_variance_ratio_), pca_plot(pc)])
    else:
        logger.warning((
            'Plots for TCC matrices are not yet supported. '
            'The HTML report will not contain any plots.'
        ))

    # Load Jinja template and render.
    with open(REPORT_TEMPLATE_PATH, 'r') as f, open(out_path, 'w') as out:
        template = Template(f.read())
        out.write(template.render(cards=cards))

    return {'report': out_path}
