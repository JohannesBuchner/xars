.PHONY: clean clean-test clean-pyc clean-build docs help
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys

try:
	from urllib import pathname2url
except:
	from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

PYTHON := python3

BROWSER := $(PYTHON) -c "$$BROWSER_PYSCRIPT"

help:
	@$(PYTHON) -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-build clean-pyc clean-test clean-doc ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +
	find . -name '*.so' -exec rm -f {} +
	find . -name '*.c' -exec rm -f {} +

clean-test: ## remove test and coverage artifacts
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/
	rm -fr .pytest_cache

clean-doc:
	rm -rf docs/build

lint: ## check style with flake8
	flake8 .

lint-fix: ## fix style with autopep8 and isort; ignores to not autofix tabs to spaces, but still warn when mixed
	autopep8 . --in-place --aggressive --aggressive --aggressive --recursive --ignore=W191,E101,E111,E122
	isort .

test: ## run tests quickly with the default Python
	PYTHONPATH=. pytest xars/binning/__init__.py xars/xsects/__init__.py xars/geometries/layeredconetorus.py xars/geometries/conetorus.py xars/geometries/wedgetorus.py xars/coordtrans.py
	echo "backend: Agg" > matplotlibrc
	$(PYTHON) scripts/vizfek2.py
	$(PYTHON) -m xars.xsects
	PYTHONPATH=. $(PYTHON) examples/torus2.py --log10nh=24.2 --opening-angle=0 --nevents=100 --output=examples/myoutput
	PYTHONPATH=. $(PYTHON) examples/disk.py --nevents=3 --output=examples/output-disk --plot-interactions --plot-paths --plot-every=40 --verbose
	cd examples/example-blobs && echo "backend: Agg" > matplotlibrc && $(PYTHON) generate_blobs.py
	PYTHONPATH=. $(PYTHON) examples/torusC.py --geometry=examples/example-blobs/torusblob23.0.hdf5 --nevents=1000
	OMP_NUM_THREADS=3 PYTHONPATH=. $(PYTHON) examples/torusC.py --geometry=examples/example-blobs/torusblob23.0.hdf5 --nevents=1000
	cd examples/example-grid && echo "backend: Agg" > matplotlibrc && $(PYTHON) generate_warpeddisk.py
	PYTHONPATH=. $(PYTHON) examples/torusG.py --geometry=examples/example-grid/warpeddisk_1.hdf5 --nevents=100
	OMP_NUM_THREADS=3 $(PYTHON) examples/torusG.py --geometry=examples/example-grid/warpeddisk_1.hdf5 --nevents=100

test-all: ## run tests on every Python version with tox
	tox

coverage: ## check code coverage quickly with the default Python
	PYTHONPATH=. coverage run --source xars -m pytest xars/binning/__init__.py xars/xsects/__init__.py xars/geometries/layeredconetorus.py xars/geometries/conetorus.py xars/geometries/wedgetorus.py xars/coordtrans.py
	echo "backend: Agg" > matplotlibrc
	PYTHONPATH=. coverage run --append --source xars scripts/vizfek2.py
	PYTHONPATH=. coverage run --append --source xars -m xars.xsects
	PYTHONPATH=. coverage run --append --source xars examples/torus2.py --log10nh=24.2 --opening-angle=0 --nevents=100 --output=examples/myoutput
	PYTHONPATH=. coverage run --append --source xars examples/disk.py --nevents=3 --output=examples/output-disk --plot-interactions --plot-paths --plot-every=40 --verbose
	cd examples/example-blobs && echo "backend: Agg" > matplotlibrc && PYTHONPATH=../.. coverage run --append --source xars generate_blobs.py
	PYTHONPATH=. coverage run --append --source xars examples/torusC.py --geometry=examples/example-blobs/torusblob23.0.hdf5 --nevents=1000
	OMP_NUM_THREADS=3 PYTHONPATH=. coverage run --append --source xars examples/torusC.py --geometry=examples/example-blobs/torusblob23.0.hdf5 --nevents=1000
	cd examples/example-grid && echo "backend: Agg" > matplotlibrc && PYTHONPATH=../.. coverage run --append --source xars generate_warpeddisk.py
	PYTHONPATH=. coverage run --append --source xars examples/torusG.py --geometry=examples/example-grid/warpeddisk_1.hdf5 --nevents=100
	OMP_NUM_THREADS=3 PYTHONPATH=. coverage run --append --source xars examples/torusG.py --geometry=examples/example-grid/warpeddisk_1.hdf5 --nevents=100
	coverage report -m
	coverage html
	$(BROWSER) htmlcov/index.html

docs: ## generate Sphinx HTML documentation, including API docs
	rm -f docs/xars.rst
	rm -f docs/modules.rst
	#nbstripout docs/*.ipynb
	sphinx-apidoc -H API -o docs/ xars
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	sed --in-place '/href="lightrayrider\/raytrace.html"/d' docs/build/html/_modules/index.html
	sed --in-place '/href="lightrayrider\/parallel.html"/d' docs/build/html/_modules/index.html
	$(BROWSER) docs/build/html/index.html

servedocs: docs ## compile the docs watching for changes
	watchmedo shell-command -p '*.rst' -c '$(MAKE) -C docs html' -R -D .

release: dist ## package and upload a release
	twine upload --verbose dist/*.tar.gz

dist: clean ## builds source and wheel package
	$(PYTHON) -m build
	ls -l dist

install: clean ## install the package to the active Python's site-packages
	$(PYTHON) -m pip install .
