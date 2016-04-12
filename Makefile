.PHONY: docs

init:
	conda install --file requirements.txt

install:
	pip install -e .

clean:
	rm -r build/

test:
	nosetests

docs-init:
	conda install --file docs/requirements.txt

docs:
	cd docs && make html

publish:
	python setup.py register
	python setup.py sdist upload
	python setup.py bdist_wheel upload
