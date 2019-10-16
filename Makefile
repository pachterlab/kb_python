.PHONY : build test

test:
	nosetests -v

build:
	python3 setup.py sdist bdist_wheel

upload-test:
	twine upload --repository-url https://test.pypi.org/legacy/ dist/*
