# kb-python
![github version](https://img.shields.io/badge/Version-0.27.2-informational)
[![pypi version](https://img.shields.io/pypi/v/kb-python)](https://pypi.org/project/kb-python/0.27.2/)
![python versions](https://img.shields.io/pypi/pyversions/kb_python)
![status](https://github.com/pachterlab/kb_python/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/pachterlab/kb_python/branch/master/graph/badge.svg)](https://codecov.io/gh/pachterlab/kb_python)
[![pypi downloads](https://img.shields.io/pypi/dm/kb-python)](https://pypi.org/project/kb-python/)
[![docs](https://readthedocs.org/projects/kb-python/badge/?version=latest)](https://kb-python.readthedocs.io/en/latest/?badge=latest)
[![license](https://img.shields.io/pypi/l/kb-python)](LICENSE)

kb-python is a python package that wraps the [kallisto | bustools](https://www.kallistobus.tools) single-cell RNA-seq workflow [2]. It was developed by Kyung Hoi (Joseph) Min and A. Sina Booeshaghi.

The wrapper simplifies downloading and running of the kallisto
[1] and bustools [2] programs. It was inspired by Sten Linnarsson’s `loompy
fromfq` command (http://linnarssonlab.org/loompy/kallisto/index.html)

The `kb` program consists of two parts:

The `kb ref` command builds or downloads a species-specific index for
pseudoalignment of reads. This command must be run prior to `kb count`, and it
runs the `kallisto index` [1].

The `kb count` command runs the kallisto [1] and bustools [2] programs. It can
be used for pre-processing of data from a variety of single-cell RNA-seq
technologies, and for a number of different workflows (e.g. production of gene
count matrices, RNA velocity analyses, etc.). The output can be saved in a
variety of formats including mtx and loom.

If you use `kb` we ask that you cite the following two papers:

[1] Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal
probabilistic RNA-seq quantification. Nature biotechnology, 34(5), 525.

[2] Melsted, P., Booeshaghi, A. S., Liu, L., Gao, F., Lu, L., Min, K. H., da Veiga Beltrame, E., Hjorleifsson, K. E., Gehring, J., & Pachter, L. (2021). Modular and efficient pre-processing of single-cell RNA-seq. Nature Biotechnology.

## Prerequisites
None. The kallisto and bustools binaries are included with the package.

## Documentation
- User documentation and tutorials are available [here](https://www.kallistobus.tools).
- Developer documentation is hosted on [Read the Docs](https://kb-python.readthedocs.io/en/latest/).
