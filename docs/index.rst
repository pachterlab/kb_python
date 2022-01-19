.. kb-python documentation master file, created by
   sphinx-quickstart on Wed Oct 30 10:42:01 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to kb-python's documentation!
=====================================

This page contains **DEVELOPER** documentation for ``kb-python`` version ``0.27.0``.
For user documentation and tutorials, please go to `kallisto | bustools <https://www.kallistobus.tools/>`_.

Development Prerequisites
^^^^^^^^^^^^^^^^^^^^^^^^^
There are a couple of things you must set up on your machine so that all of your
commits satisfy code quality and unit-testing requirements. First, install all
necessary packages by running::

  pip install -r requirements.txt
  pip install -r dev-requirements.txt

Code qualty and unit tests are strictly enforced for every pull request via
Github actions.

Code Quality
""""""""""""
``kb-python`` uses ``flake8`` and ``yapf`` to ensure code quality. The easiest
way to set these up so that they run automatically for every commit is to install
``pre-commit`` hooks by running::

  pre-commit install

at the root of the repository.

Unit-testing
""""""""""""
``kb-python`` uses ``nose`` to run unit tests. There is a convenient Makefile
rule in place to run all tests.::

  make test

Releasing New Versions
^^^^^^^^^^^^^^^^^^^^^^
This section walks you through, step-by-step, how to release a new version.

1. Make sure you are on the up-to-date ``master`` branch.
2. Run ``make bump_patch``, ``make bump_minor``, or ``make bump_major`` depending
   on what version you will be bumping.
3. Run ``make push_release``. This will push the new commit and tag.
4. Go to the `releases` tab on Github.
5. Select the new release, edit the release description, and `Publish release`.
6. A Github action will automatically trigger to upload the new release to PyPi.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
