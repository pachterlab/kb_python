.PHONY : test check build clean upload-test

test:
	nosetests -v -s

check:
	flake8 kb_python && echo OK
	yapf -r --diff kb_python && echo OK

build:
	python3 setup.py sdist bdist_wheel

clean:
	rm -r build
	rm -r dist

bump_patch:
	bumpversion patch

bump_minor:
	bumpversion minor

bump_major:
	bumpversion major

upload-test:
	twine upload --repository-url https://test.pypi.org/legacy/ dist/*
