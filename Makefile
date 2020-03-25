.PHONY: clean
clean: # get rid of all build artifacts
	@rm -rf arpy.egg-info build dist

.PHONY: get-poetry
get-poetry: # locally install poetry
	@curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python3

.PHONY: poetry-install
poetry-install: # install local poetry dependencies in a venv
	@poetry install

.PHONY: format
format: # format the code base using isort and black
	@poetry run isort -y
	@poetry run black $(PWD)

.PHONY: test
test: # run the test suite
	@poetry run pytest --disable-warnings -v --black --cov=arpy --cov-config .coveragerc --cov-report term-missing:skip-covered

.PHONY: wheel
wheel: # build a wheel using poetry
	@poetry build --format wheel
