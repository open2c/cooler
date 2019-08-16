# Contributing

## General guidelines

If you haven't contributed to open-source before, we recommend you read [this excellent guide by GitHub on how to contribute to open source](https://opensource.guide/how-to-contribute). The guide is long, so you can gloss over things you're familiar with.

If you're not already familiar with it, we follow the [fork and pull model](https://help.github.com/articles/about-collaborative-development-models) on GitHub. Also, check out this recommended [git workflow](https://www.asmeurer.com/git-workflow/).


## Contributing Code

This project has a number of requirements for all code contributed.

* We follow the [PEP-8 style](https://www.python.org/dev/peps/pep-0008/) convention.
* We use [Numpy-style docstrings](https://numpydoc.readthedocs.io/en/latest/format.html).
* It's ideal if user-facing API changes or new features have documentation added.


## Setting up Your Development Environment

After forking and cloning the repository, install in "editable" (i.e. development) mode using the `-e` option:

```sh
git clone https://github.com/mirnylab/cooler.git
cd cooler
pip install -e .[all]
```

Editable mode installs the package by creating a "link" to the working (repo) directory.


## Running/Adding Unit Tests

It is best if all new functionality and/or bug fixes have unit tests added with each use-case.

We use [pytest](https://docs.pytest.org/en/latest) as our unit testing framework. Once you've configured your environment, you can just `cd` to the root of your repository and run

```sh
pytest
```

Unit tests are automatically run on Travis CI for pull requests.


## Adding/Building the Documentation

If a feature is stable and relatively finalized, it is time to add it to the documentation. If you are adding any private/public functions, it is best to add docstrings, to aid in reviewing code and also for the API reference.

We use [Numpy style docstrings](https://numpydoc.readthedocs.io/en/latest/format.html>) and [Sphinx](http://www.sphinx-doc.org/en/stable) to document this library. Sphinx, in turn, uses [reStructuredText](http://www.sphinx-doc.org/en/stable/rest.html) as its markup language for adding code.

We use the [Sphinx Autosummary extension](http://www.sphinx-doc.org/en/stable/ext/autosummary.html) to generate API references. You may want to look at `docs/api.rst` to see how these files look and where to add new functions, classes or modules.

To build the documentation:

```sh
make docs
```

After this, you can find an HTML version of the documentation in `docs/_build/html/index.html`.

Documentation from `master` and tagged releases is automatically built and hosted thanks to [readthedocs](https://readthedocs.org/).


## Acknowledgments

This document is based off of the [guidelines from the sparse project](https://github.com/pydata/sparse/blob/master/docs/contributing.rst).



<!-- with the `pytest-cov` extension to check code coverage and `pytest-flake8` to check code style. You don't need to configure these extensions yourself.
This automatically checks code style and functionality, and prints code coverage, even though it doesn't fail on low coverage. -->


<!-- ## Coverage

The `pytest` script automatically reports coverage, both on the terminal for missing line numbers, and in annotated HTML form in `htmlcov/index.html`.
Coverage is automatically checked on CodeCov for pull requests. -->

<!-- 
## Adding and Running Benchmarks

We use [Airspeed Velocity](https://asv.readthedocs.io/en/latest/)  to run benchmarks. We have it set up to use `conda`, but you can edit the configuration locally if you so wish.
 -->

<!-- 
[Pull requests](https://akrabat.com/the-beginners-guide-to-contributing-to-a-github-project/) are welcome. The current requirements for testing are `pytest` and `mock`.

For development, clone and install in "editable" (i.e. development) mode with the `-e` option. This way you can also pull changes on the fly.
```sh
$ git clone https://github.com/mirnylab/cooler.git
$ cd cooler
$ pip install -e . -->