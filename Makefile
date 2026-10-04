# Myelovar build/run/test/demo targets.
PYTHON := .venv/bin/python
PYTEST := .venv/bin/pytest
PIP := .venv/bin/pip

.PHONY: env data cache run run-b test api web build-web demo docker clean-data clean-outputs

env:                                ## create venv + install package (idempotent)
	@test -x $(PYTHON) || (python3.12 -m venv .venv || python3 -m venv .venv)
	$(PIP) install -q -e ".[api,test]"

data: env                           ## download + verify all public data, write manifest
	$(PYTHON) scripts/download_data.py

cache: data                         ## precompute reference caches (AF extract, indexes)
	$(PYTHON) scripts/build_cache.py

run: cache                          ## Mode A end-to-end -> outputs/
	$(PYTHON) -m myelovar run --mode A

run-b: cache                        ## Mode B end-to-end -> outputs/mode_b/
	$(PYTHON) -m myelovar run --mode B

test: env
	$(PYTHON) -m pytest

build-web:                          ## generate OpenAPI types + build React app
	cd web && npm install --no-audit --no-fund && npm run build

web: build-web
	@echo "static build in web/dist (served by api at /)"

api: env
	$(PYTHON) -m uvicorn api.main:app --host 127.0.0.1 --port 8000

demo: cache build-web               ## run precomputed demo and open browser
	$(PYTHON) scripts/make_demo_bundle.py
	$(PYTHON) -m myelovar serve --open

docker:
	docker compose up --build

clean-data:                         ## remove downloaded bulk data (keeps manifest)
	rm -rf data/raw data/tracks/encode data/tracks/ccre data/reference/genome

clean-outputs:
	rm -rf outputs/* demo_bundle/*
