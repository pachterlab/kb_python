# kb_python
Version: v0.0.2
![](https://github.com/pachterlab/kb_python/workflows/CI/badge.svg)

A friendly wrapper around the kallisto | bustools pipeline for scRNA-seq analysis.

## Prerequisites
- [kallisto](https://pachterlab.github.io/kallisto/) v0.46.0 and up.
Must be accessible from the command-line as `kallisto`.
- [bustools](https://bustools.github.io/) v0.39.3 and up.
Must be accessible from the command-line as `bustools`.

There are plans to include installers for both prerequisites.

## Development
`kb_python` uses `flake8` and `yapf` to ensure code quality and `nose`
to run unittests. All necessary dependencies for development can be installed
by running `pip install -r dev-requirements.txt`.

The CI workflow ensures all code passes code quality checks and unittests.
It is recommended to use `pre-commit` to make sure each commit satisfies
code quality specifications. To do so, first install `pre-commit` by running
`pip install pre-commit`, and then at the root run `pre-commit install`.
Every future commit will pass through `flake8` and `yapf`.
