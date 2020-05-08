POETRY_BIN := $(shell which poetry 2> /dev/null)

.PHONY: clean
clean: # get rid of all build artifacts and caches
	rm -rf arpy.egg-info build dist **/__pycache__ .pytest_cache

.PHONY: ensure-poetry
ensure-poetry:
ifndef POETRY_BIN
	@echo "fetching poetry for dependency management..."
	curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python3
endif

.PHONY: poetry-install
poetry-install: ensure-poetry # install local poetry dependencies in a venv
	poetry install

.PHONY: format
format: ensure-poetry # format the code base using isort and black
	poetry run isort -y
	poetry run black $(PWD)

.PHONY: test
test: ensure-poetry # run the test suite
	poetry run pytest --disable-warnings -v --black --cov=arpy --cov-config .coveragerc --cov-report term-missing:skip-covered

.PHONY: wheel
wheel: ensure-poetry # build a wheel using poetry
	rm -rf dist
	poetry build --format wheel

.PHONY: install
install: wheel # install arpy locally
	python3 -m pip install dist/arpy-*-py3-none-any.whl
