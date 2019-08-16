.PHONY: init install clean-pyc clean-build build test publish docs-init docs

init:
	conda install --file requirements.txt

install:
	pip install -e .

test:
	pytest

clean-pyc:
	find . -name '*.pyc' -exec rm --force {} +
	find . -name '*.pyo' -exec rm --force {} +
	find . -name '*~' -exec rm --force  {} +

clean-build:
	rm -rf build/
	rm -rf dist/

clean: clean-pyc clean-build

build: clean-build
	python setup.py sdist
	python setup.py bdist_wheel

publish: build
	twine upload dist/*

publish-test:
	twine upload --repository-url https://test.pypi.org/legacy/ dist/*

docs-init:
	conda install --file docs/requirements.txt

docs:
	cd docs && python make_cli_rst.py && make html
