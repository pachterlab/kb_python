# kb_python
![](https://img.shields.io/badge/Version-0.2.0-informational)
![](https://github.com/pachterlab/kb_python/workflows/CI/badge.svg)

A friendly wrapper around the kallisto | bustools pipeline for scRNA-seq analysis.

## Prerequisites
- ~~[kallisto](https://pachterlab.github.io/kallisto/) v0.46.0 and up.
Must be accessible from the command-line as `kallisto`.~~
- ~~[bustools](https://bustools.github.io/) v0.39.3 and up.
Must be accessible from the command-line as `bustools`.~~

~~There are plans to include installers for both prerequisites.~~

Binaries are included with the package starting version `0.0.8`. There are no prerequisites starting from this version.

## Development
### Code Quality
`kb_python` uses `flake8` and `yapf` to ensure code quality and `nose`
to run unittests. All necessary dependencies for development can be installed
by running `pip install -r dev-requirements.txt`.

The CI workflow ensures all code passes code quality checks and unittests.
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
2. Push the commit and tag with `git push --tag`.
3. Go to the `releases` tab on Github. Select the version that was just commited.
`Edit tag`, write a description, and `Publish release`.
4. A Github Actions workflow will be triggered to build and upload the updated
package to Pypi.
