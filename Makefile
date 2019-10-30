.PHONY : install test check build docs clean push_release

test:
	nosetests --verbose --with-coverage --cover-package kb_python

check:
	flake8 kb_python && echo OK
	yapf -r --diff kb_python && echo OK

build:
	python setup.py sdist bdist_wheel

docs:
	sphinx-build -a docs docs/_build

clean:
	rm -rf build
	rm -rf dist
	rm -rf kb_python.egg-info
	rm -rf docs/_build
	rm -rf docs/api

bump_patch:
	bumpversion patch

bump_minor:
	bumpversion minor

bump_major:
	bumpversion major

push_release:
	git push && git push --tags
