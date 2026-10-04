# How Myelovar works — behind the scenes

One paragraph: Myelovar downloads public genomics data, indexes it, filters a genome's
variants down to rare ones sitting inside regulatory DNA that is *active in myeloma
cells*, computes how each one disrupts a transcription-factor motif, traces which gene
the disrupted element controls, asks whether myeloma cells depend on that gene, then
ranks everything with a deterministic score and writes a human-readable evidence chain
for each result. Below is that story with the actual numbers and mechanisms.

---

## Stage 0 — Ingest (`s0_load_input.py`)

Read the input: a VCF for Mode A (demo: GIAB HG001, **3,893,341 records**; the funnel
sees **3,935,734** after normalization/dedup handling), or the Mode B GWAS table.
Records are normalized to `chrom/pos/ref/alt`, restricted to autosomes, and split so
indels get handled separately from SNVs later.

## Stage 1 — Annotate + filter (`s1_annotate_filter.py`)

Two questions per variant:

1. **Is it too common?** If population allele frequency ≥ 0.01 it is germline noise,
   not a candidate driver. AF comes from a prewarmed dbSNP `snp151Common` lookup
   (`cache/af`); **2,788,685 / 3,935,734** variants matched. Missing AF = keep
   (we don't drop what we can't measure).
2. **Is it good quality?** DP ≥ 10 and GQ ≥ 20 when the VCF carries them (GIAB
   does; the code warns and skips rather than fails if a VCF lacks FORMAT fields).

Mode A: **3,935,734 → 1,137,368** (2.79M removed). Mode B skips this stage — loci
were *chosen* for association, so common alleles inside them can still be causal.

## Stage 2 — Regulatory overlap (`s2_regulatory_overlap.py`)

A variant can only produce an evidence chain if some assay says its position does
something. We intersect every surviving variant with **18 ENCODE MM.1S (myeloma cell
line) tracks** — ATAC/DNase accessibility, H3K27ac/H3K4me1/H3K4me3 histone marks,
ROSE super-enhancers (591), TF ChIP peaks, and SCREEN cCREs. Tracks are loaded as
interval indexes (built by `scripts/build_cache.py`), so overlap is a sorted-merge
query per track, not a file scan.

Mode A: **1,137,368 → 48,113** inside ≥ 1 track. The surviving set also gets a
0–1 `regulatory_context_score` — a weighted blend of which evidence types are present
(accessibility 0.30, H3K27ac 0.25, super-enhancer 0.15, H3K4me1/H3K4me3 0.10 each,
TF peak 0.10).

## Stage 3 — Motif effect (`s3_motif_effect.py`) — the biology core

For each variant, with `W = 30 bp` of reference context on each side:

1. **Read the real DNA** from the hg38 FASTA (bgzip + faidx) — both strands via the
   reference, not by guessing.
2. **Build the alternate sequence** by splicing in the alt allele. SNVs keep the
   window length; indels change it, so windows are **grouped by length** and
   one-hot-encoded in batches (mixed-length windows cannot be tensored together).
3. **Score both strands.** For each of **10 JASPAR PWMs** (IRF4, PRDM1, MAF, MAFB,
   MYC, MAX, XBP1, RELA, NFKB1, CTCF): score the *forward* window against the PWM and
   the *forward* window against the reversed-complemented PWM — never reverse the
   sequence itself. Vectorized as a single `einsum("npb,bp->n")`.
4. **Classify the change** across all sliding windows on both strands, taking the
   best-scoring window per strand:
   - `broken` — reference scored ≥ 80% of PWM max and the alt loses ≥ 20 points
   - `created` — the mirror image (alt crosses that 80% bar, ref didn't)
   - `changed` — a real shift that misses the bar, `none` — noise

Mode A result: **169 broken, 140 created** motifs. This step shipped two bugs that
the test suite caught and now guards: swapped one-hot axes, and a reverse-complement
scan that scored the wrong strand. A brute-force reference scorer cross-checks the
vectorized scanner in `tests/test_motif.py`.

## Stage 4 — Gene link (`s4_gene_link.py`)

A disrupted motif only matters if you can name the gene it acts on. Priority order:

1. **ABC (activity-by-contact)** — element activity (stage-2 score) × genomic contact
   decay, matched to genes within the window: **357** variants. Caveat, stated
   everywhere: the contact maps are a **B-cell ENCODE proxy** (GM12878/Roadmap) —
   no plasma-cell Hi-C exists.
2. **Distance fallback** — target within 100 kb by TSS distance with exponential
   decay (10 kb scale): **46,030**.
3. **Nearest gene** regardless of distance: **2,083**.

Each link records *which method* assigned it, so a reviewer can tell strong links
from fallbacks.

## Stage 5 — Myeloma relevance (`s5_myeloma_relevance.py`)

Now ask whether the target gene matters *to myeloma*:

- **DepMap dependency** — mean Chronos gene effect ≤ −0.5 across **19 myeloma cell
  lines**, with a ≥ 0.30 selectivity margin vs. other lines (weight 0.40)
- **Known myeloma gene** — curated set (IRF4, PRDM1, MYC, …) (weight 0.30)
- **Expression** — TPM > 1 in ≥ 50% of myeloma lines (weight 0.15)
- **GWAS support** — target gene lies in/near a myeloma risk locus (weight 0.15)

Output: per-variant `depmap / known / expressed / gwas` flags plus a 0–1
`myeloma_relevance_score`.

## Stage 6 — Score + rank (`s6_score_rank.py`)

A deterministic weighted sum — no model, no randomness, same input → same output:

```
score = 0.25 · regulatory_context + 0.20 · motif_disruption + 0.20 · link
      + 0.25 · myeloma_relevance + 0.10 · gwas
```

All weights live in `config.yaml`. The 48,113 Mode A survivors are sorted and written
to `outputs/result.json` (5,000 shown in the API payload, all rows in
`outputs/ranked_variants.csv`, 61 MB). Perturbing every weight and re-ranking keeps
**Spearman ρ = 0.976** — the ordering is driven by evidence, not by any single knob.

## The evidence chain (the actual product)

Each ranked variant carries, and the UI renders:

```
variant chrX:pos REF>ALT
  ├─ regulatory element   which of the 18 tracks overlap + context score
  ├─ motif effect         which TF, direction (broken/created/changed), Δscore, strand
  ├─ target gene          gene symbol + link method (ABC/distance/nearest) + distance
  ├─ myeloma relevance    DepMap dependency/selectivity, known-gene, expression, GWAS
  └─ explanation          one generated sentence summarizing the above
```

## Mode B — same funnel, different entry

`src/myelovar/io/gwas.py` builds a table from **37 GWAS Catalog myeloma loci / 64
unique lead variants**: leads are labelled `lead`, LD proxies from Ensembl
(r² ≥ 0.6, EUR, ≤ 20 per lead, polite rate limiting) are labelled `proxy`, and every
variant in a ±50 kb window gets `window` (1,861,752 rows total; 1,230 proxies).
hg19→hg38 via the UCSC chain file. The funnel then runs s0–s6 on those rows with
`remove_common` disabled.

## Validation — does the ranking mean anything?

`scripts/run_validation.py`:

- **Independent check:** top-ranked variants are enriched in independently defined
  cCREs — **OR = 3.85, p = 1.6e-78** (cCREs never touch scoring, so this is fair).
- **Partially circular check:** enrichment for known myeloma genes — **OR = 5.46,
  p = 2.4e-12**, but those genes *do* inform stage 5, so this number is labelled
  partially circular and weighted 0.25 in the summary. We show it because hiding it
  would be worse.
- **Robustness:** weight perturbation sweep (above).
- **Mode B:** lead-SNP recovery — do ranked lists contain the original lead variants?
- Figures: funnel, score histogram, enrichment — in `validation/figures/`.

## Provenance

`data/data_manifest.json` records every resource: URL, SHA hash, byte size, status,
and a content hash for the whole manifest (`62ead693212ebba9`, 32 resources,
31 ok / 1 failed / 0 critical). Every result payload embeds the manifest hash, config
used, and step counts, so any ranked list can be traced to exact inputs and settings.
The one failure (GEO GSE160335, 2.8 GB RAW.tar, individual files 404) is recorded
rather than hidden; it is optional data.
