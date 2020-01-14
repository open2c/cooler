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
	rm -f .coverage
	rm -f coverage.xml
	rm -rf htmlcov/

clean: clean-pyc clean-build

build: clean-build
	python setup.py sdist
	python setup.py bdist_wheel

publish: build
	twine upload dist/*

publish-test: build
	twine upload --repository-url https://test.pypi.org/legacy/ dist/*
	# pip install --extra-index-url https://test.pypi.org/simple/ cooler

docs-init:
	conda install --file docs/requirements.txt

docs:
	cd docs && python make_cli_rst.py && make html
