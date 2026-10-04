# Myelovar demo image: API + built web UI + precomputed demo bundle.
# Full data (data/, 4.9G) is NOT baked in; mount it via docker-compose for live runs.
FROM python:3.12-slim

ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    MYELOVAR_WEB_DIST=/app/web/dist \
    MYELOVAR_DEMO_DIR=/app/demo_bundle \
    MYELOVAR_VALIDATION_DIR=/app/validation \
    MYELOVAR_OUTPUTS_DIR=/app/outputs

WORKDIR /app

# Wheels for pysam/pyBigWig/scipy ship manylinux aarch64+x86_64; gcc only as fallback.
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc libc6-dev zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Dependency layer first (better cache): install the package itself.
COPY pyproject.toml README.md ./
COPY src ./src
RUN pip install -e ".[api]"

# App + precomputed artifacts.
COPY api ./api
COPY config.yaml ./
COPY web/dist ./web/dist
COPY validation ./validation
COPY demo_bundle ./demo_bundle
COPY scripts ./scripts

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
    CMD python -c "import urllib.request,sys; r=urllib.request.urlopen('http://127.0.0.1:8000/api/health', timeout=4); sys.exit(0 if r.status==200 else 1)"

CMD ["python", "-m", "uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000"]
