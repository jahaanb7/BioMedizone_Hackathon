# Myelovar

**Myeloma regulatory variant prioritization.** Takes a whole-genome VCF (Mode A) or a
set of GWAS loci (Mode B) and ranks *non-coding* variants by a complete evidence chain:

```
variant → active regulatory element in MM.1S myeloma cells → TF motif disrupted
        → target gene → myeloma relevance → plain-English explanation
```

Research use only. **Not a diagnostic.**

## Quick start

```bash
make env          # create .venv, install package (+api,test extras)
make data         # download + verify all public datasets -> data/, manifest
make cache        # precompute indexes (AF, ABC, super-enhancers, track indexes)
make run          # Mode A end-to-end -> outputs/
make test         # 33 unit/integration tests
make build-web    # OpenAPI types + React build -> web/dist
make api          # REST API + web UI on http://127.0.0.1:8000
make demo         # precomputed offline demo (no compute needed)
```

Fastest path (no downloads, no compute): `make api` and click **Demo A** in the UI —
the demo bundle ships precomputed results.

## The two modes

| | Mode A (genome funnel) | Mode B (GWAS loci) |
|---|---|---|
| Input | VCF (GIAB HG001 demo) | GWAS Catalog myeloma loci (37 loci / 64 unique leads) |
| Question | "scan a whole genome for rare regulatory variants" | "what causal variants hide inside known risk loci?" |
| Scope | whole genome | ±50 kb around each lead SNP + LD proxies (r² ≥ 0.6, EUR) |
| Output | `outputs/result.json`, `outputs/ranked_variants.csv` | `outputs/mode_b/…` |

Both modes run the same six-step funnel (Mode B skips the common-variant prune).

## Pipeline (what runs, in order)

1. **s0 load input** — read VCF/GWAS records, normalize alleles, autosomes only.
2. **s1 annotate + filter** — drop common variants (dbSNP `snp151Common` AF ≥ 0.01;
   missing AF = keep), apply DP/GQ when the VCF carries them. Mode A: 3.94M → 1.14M.
3. **s2 regulatory overlap** — keep variants inside any of **18 ENCODE MM.1S tracks**
   (ATAC/DNase, H3K27ac/H3K4me1/H3K4me3, super-enhancers, TF peaks, SCREEN cCREs).
   → 48,113 survivors. A variant outside all regulatory data has no story to tell.
4. **s3 motif effect** — extract ±30 bp reference window from the hg38 FASTA, build the
   alt sequence, score both strands against JASPAR PWMs for 10 myeloma TFs
   (IRF4, PRDM1, MAF, MAFB, MYC, MAX, XBP1, RELA, NFKB1, CTCF) and classify the
   change: `broken` / `created` / `changed` / `none` (needs ≥ 80% of PWM max on the
   stronger allele and ≥ 20-point loss).
5. **s4 gene link** — assign each element its target gene: **ABC model** first
   (activity × contact; B-cell proxy biosample — labelled as such), then distance
   within 100 kb, else nearest gene.
6. **s5 myeloma relevance** — DepMap CRISPR dependency + selectivity across 19 myeloma
   cell lines, known-myeloma-gene membership, expression, GWAS-locus proximity.
7. **s6 score + rank** — deterministic weighted sum (weights in `config.yaml`), sorted,
   written to `outputs/`.

Every variant carries its full evidence chain and a generated explanation
(`src/myelovar/explain.py`).

## Repository layout

```
src/myelovar/        importable, UI-free core
  steps/             s0..s6 pipeline stages
  io/                VCF, GTF, bigWig, BED, ABC, DepMap, GWAS/Ensembl readers
  reference.py       FASTA access + one-hot encoding + PWM scanning
  schemas.py         Pydantic result/variant/funnel schemas (round-trip tested)
  pipeline.py        orchestration + manifest/provenance
api/                 FastAPI thin layer (runs, variants, locus, validation, SPA)
web/                 Vite + React + TS UI (funnel, table, variant card, locus plot)
scripts/             download_data, build_cache, run_validation, make_demo_bundle
tests/               33 tests (motif brute-force cross-check, funnel invariants, API)
validation/          enrichment results + figures
demo_bundle/         precomputed offline demo results
config.yaml          all thresholds and weights — nothing is hard-coded in steps
```

## API (served by `make api`, SPA at `/`)

```
GET  /api/health                       readiness + data/demo flags
GET  /api/manifest                     data provenance (SHA hashes, statuses)
GET  /api/modes                        available demo modes
GET  /api/demo/{mode}                  load precomputed demo A or B
POST /api/runs                         start a pipeline run
GET  /api/runs/{id}                    job status
GET  /api/runs/{id}/funnel             stage-by-stage counts
GET  /api/runs/{id}/variants           paged + filtered variant list
GET  /api/runs/{id}/variants/{vid}     one variant with full evidence chain
GET  /api/runs/{id}/locus              tracks + genes + motifs for a region
GET  /api/runs/{id}/export.csv         CSV export
GET  /api/genes/{symbol}               gene card
GET  /api/validation, /api/validation/figures/{name}
```

## Data & honesty notes

`make data` writes `data/data_manifest.json` (content hash `62ead693212ebba9`):
**32 resources, 31 ok, 1 failed, 0 critical.** Deviations from the original plan,
recorded deliberately:

- **gnomAD → dbSNP `snp151Common`** for the common-variant prune (public gnomAD
  AF files were not fetched; dbSNP common is the shipped substitute).
- **DepMap from figshare** (Model.csv, CRISPRGeneEffect, expression) rather than the
  original portal URL scheme.
- **GWAS Catalog v2** schema (lead/proxy rows keyed by `variation1/variation2`).
- **ABC model uses a B-cell ENCODE proxy biosample** — no plasma-cell contact
  dataset exists; every payload labels it `proxy`.
- **GEO GSE160335 not downloaded**: only a 2.8 GB RAW.tar was reachable and individual
  files 404'd; refused rather than silently skipping. Marked `failed`, non-critical.

## Validation (independent, `make` → `scripts/run_validation.py`)

- cCRE enrichment of top-ranked variants: **OR = 3.85, p = 1.6e-78** (independent of scoring)
- known-myeloma-gene enrichment: **OR = 5.46, p = 2.4e-12** — *partially circular*
  (those genes inform scoring), so it is weighted 0.25 in the summary
- rank robustness under weight perturbation: **Spearman 0.976**

## Tests

```bash
make test          # 33 passed
```

Includes a brute-force cross-check of the vectorized motif scanner, both-strand
orientation tests, indel window grouping, schema round-trips against real output,
and API tests through FastAPI's TestClient.
