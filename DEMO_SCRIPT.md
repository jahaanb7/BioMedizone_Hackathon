# Myelovar demo script (~8 minutes)

Two paths. **Path 1 (offline demo)** uses precomputed results — no downloads, no
compute, works anywhere. **Path 2 (live run)** executes the real funnel in ~4 minutes
so you can show the numbers moving.

---

## Path 1 — offline demo (recommended for talks)

```bash
make env && make build-web && make api
# open http://127.0.0.1:8000
```

### 1. Dashboard — the funnel (30 s)

Click **Demo A**. Point at the funnel chart:

> "3.9 million raw genome variants in, 48,113 rare variants inside regulatory DNA
> that is actually active in myeloma cells. Every stage removes a class of noise:
> common → non-regulatory → non-disruptive."

### 2. Variants table — filters as story beats (90 s)

- Sort by score; show the top rows.
- Filter `motif_effect = broken` → **169 variants**. "These are the ones that break a
  transcription-factor binding site — IRF4, PRDM1, MYC, NF-κB — the factors that run
  myeloma."
- Filter `evidence.gwas = true` → risk-locus variants.
- Mention: paginated server-side, 48,113 rows, TanStack table, CSV export.

### 3. Variant detail card — the evidence chain (2 min)

Open a top variant. Walk the chain top to bottom:

1. **Element** — which of the 18 MM.1S tracks overlap + regulatory context score.
2. **Motif effect** — which TF, broken vs created, score delta, strand. "Both DNA
   strands are scored: forward sequence vs the PWM, and forward sequence vs the
   reversed-complemented PWM."
3. **Target gene** — symbol, distance, and *how* we linked it: ABC contact model first
   (labelled B-cell proxy), then distance, then nearest.
4. **Myeloma relevance** — DepMap dependency across 19 myeloma lines, known-gene flag,
   expression, GWAS support.
5. **Explanation** — the generated sentence. "This is the deliverable: a reviewer can
   challenge any link in the chain."

### 4. Locus plot (60 s)

Open the locus view: overlapping tracks as intervals, genes underneath, variants as
points sized by score. "Evidence lives in *regions*, not single positions — this is
what a reviewer sees when they ask 'what else is near this variant?'"

### 5. Validation page (60 s)

- cCRE enrichment **OR 3.85, p = 1.6e-78** — independent of scoring, so it's a fair test.
- Known-myeloma-gene **OR 5.46, p = 2.4e-12** — explicitly flagged *partially circular*
  and weighted down. "We label the circular check instead of hiding it."
- Weight-perturbation sweep: **Spearman 0.976** — the ranking isn't an artifact of one
  weight.

### 6. Manifest page (30 s)

32 resources, 31 ok, 1 failed (a 2.8 GB GEO tarball we refused to guess at), SHA hashes,
content hash `62ead693212ebba9`. "Every result embeds this hash — reproducible or
nothing."

### 7. Mode B (30 s)

Click **Demo B**: same funnel entered from 37 known GWAS loci + LD proxies instead of
a whole genome. "Mode A scans a genome; Mode B dissects known risk loci."

---

## Path 2 — live run

```bash
make data        # ~downloads + manifest
make cache       # indexes: AF, ABC, super-enhancers, 18 track indexes
make run         # Mode A, ~4 min on a laptop
make api
```

Or one mode at a time in the UI: hit **Run** on the dashboard, watch job status
transition, land on the run page when it finishes. Console equivalent:

```bash
.venv/bin/python -m myelovar run --mode A    # outputs/result.json + ranked_variants.csv
.venv/bin/python -m myelovar run --mode B    # outputs/mode_b/
.venv/bin/python -m pytest                   # 33 tests
```

---

## Numbers worth having in your head

| | |
|---|---|
| Raw variants | 3,935,734 |
| After common-variant prune | 1,137,368 |
| Inside MM.1S regulatory elements | 48,113 |
| Motifs broken / created | 169 / 140 |
| Target links: ABC / distance / nearest | 357 / 46,030 / 2,083 |
| MM.1S tracks | 18 |
| TFs scanned | 10 JASPAR PWMs, both strands |
| Myeloma cell lines (DepMap) | 19 |
| Validation | cCRE OR 3.85 · gene OR 5.46 (partially circular) · ρ 0.976 |
| Tests | 33 passed |

## Anticipated questions

- **"Is this a diagnostic?"** No. Research prioritization only; every page says so.
- **"Where does gene linking come from?"** ABC with a *B-cell proxy* — labelled in
  every payload, because no plasma-cell contact map exists.
- **"Why not gnomAD?"** The shipped prune uses dbSNP `snp151Common`; the substitution
  is recorded in the manifest and README.
- **"Isn't the known-gene enrichment circular?"** Partially — which is why it's
  flagged and down-weighted, while the independent cCRE test is the headline number.
