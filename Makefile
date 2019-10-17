.PHONY : install test check build clean upload-test

KALLISTO_VERSION=v0.46.0
BUSTOOLS_VERSION=v0.39.3

install:
	wget https://github.com/pachterlab/kallisto/releases/download/v0.46.0/kallisto_linux-$(KALLISTO_VERSION).tar.gz -O kallisto.tar.gz
	tar -xvzf kallisto.tar.gz
	sudo cp kallisto/kallisto /usr/local/bin

	wget https://github.com/BUStools/bustools/releases/download/v0.39.3/bustools_linux-$(BUSTOOLS_VERSION).tar.gz -O bustools.tar.gz
	tar -xvzf bustools.tar.gz
	sudo cp bustools/bustools /usr/local/bin

test:
	nosetests --verbose --with-coverage --cover-package kb_python

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
