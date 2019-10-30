# kb-python
![github version](https://img.shields.io/badge/Version-0.2.1-informational)
[![pypi version](https://img.shields.io/pypi/v/kb-python)](https://pypi.org/project/kb-python/0.2.1/)
![python versions](https://img.shields.io/pypi/pyversions/kb_python)
![status](https://github.com/pachterlab/kb_python/workflows/CI/badge.svg)
[![pypi downloads](https://img.shields.io/pypi/dm/kb-python)](https://pypi.org/project/kb-python/)
[![docs](https://readthedocs.org/projects/kb-python/badge/?version=latest)](https://kb-python.readthedocs.io/en/latest/?badge=latest)
[![license](https://img.shields.io/pypi/l/kb-python)](LICENSE)

A wrapper for the [kallisto | bustools](https://www.kallistobus.tools) single-cell RNA-seq workflow.

## Prerequisites
None. The kallisto and bustools binaries are included with the package.

## Getting Started
Visit the [Getting Started](https://www.kallistobus.tools/kb_getting_started) page.

## Tutorials
- WIP...

## Development
### Documentation
Developer documentation is hosted on [Read the Docs](https://kb-python.readthedocs.io/en/latest/).

### Code Quality
`kb-python` uses `flake8` and `yapf` to ensure code quality and `nose`
to run unittests. All necessary dependencies for development can be installed
by running `pip install -r dev-requirements.txt`.

The CI workflow ensures all code passes code quality checks and unit tests.
It is recommended to use `pre-commit` to make sure each commit satisfies
code quality specifications. To do so, first install `pre-commit` by running
`pip install pre-commit`, and then at the root run `pre-commit install`.
Every future commit will pass through `flake8` and `yapf`.

### Bumpversion
Bumping versions is done with `bumpversion`. This should be installed from the
`dev-requirements.txt`, but can be installed separately with `pip`. To bump
version and release the new version to Pypi,
1. Run `make bump_patch`, `make bump_minor` or `make bump_major` depending
on which version to bump. This will make a new commit and create a new tag
with the new version.
2. Push the commit and tag with `make push_release`.
3. Go to the `releases` tab on Github. Select the version that was just commited.
`Edit tag`, write a description, and `Publish release`.
4. A Github Actions workflow will be triggered to build and upload the updated
package to Pypi.
