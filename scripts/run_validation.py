"""Statistical validation of the ranked outputs, computed on REAL data only.

Writes validation/results.json + validation/figures/*.png, served by the API
at /api/validation. Every test states whether it is independent of the scoring
function (nothing circular is presented as independent evidence).

Tests
-----
1. cCRE enrichment: are top-1000 ranked variants more often inside ENCODE
   SCREEN cCREs than the rest of the ranked set? (cCRE registry is annotated
   but NOT part of any score component -> independent.)
2. Known-myeloma-gene enrichment: are genes linked to the top-500 variants
   enriched for IntOGen PCM drivers / DepMap-selective genes vs genes linked
   to the rest? (PARTIALLY CIRCULAR: the myeloma-relevance component uses
   known-gene/DepMap flags with weight 0.25 - labelled as such.)
3. Mode B lead-SNP recovery: is the published GWAS lead SNP ranked in the
   top-3 of its locus? (NOT circular: every candidate in a locus carries the
   same GWAS flag, so the GWAS component cannot favour the lead.)
4. Weight-sensitivity summary from the pipeline's own perturbation analysis.
"""
from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest, fisher_exact

log = logging.getLogger("validate")

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation"
FIG = OUT / "figures"


def _split_genes(val) -> list[str]:
    if not isinstance(val, str) or not val:
        return []
    return [g.strip() for g in val.split(",") if g.strip()]


def ccre_enrichment(df: pd.DataFrame, top_n: int = 1000) -> dict:
    """Fisher: top-N ranked in SCREEN cCRE vs the rest of the ranked set."""
    df = df.sort_values("score", ascending=False)
    top = df.head(top_n)
    rest = df.iloc[top_n:]
    a = int(top["in_encode_ccre"].fillna(False).astype(bool).sum())
    b = len(top) - a
    c = int(rest["in_encode_ccre"].fillna(False).astype(bool).sum())
    d = len(rest) - c
    odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return {
        "test": f"top{top_n}_vs_rest__inside_ENCODE_SCREEN_cCRE",
        "contingency": {"top_in": a, "top_out": b, "rest_in": c, "rest_out": d},
        "odds_ratio": float(odds),
        "p_value": float(p),
        "alternative": "greater",
        "independent_of_scoring": True,
        "note": "cCRE membership is annotated for display but is not a "
                "component of any score; peaks (not cCREs) drive selection.",
    }


def gene_enrichment(df: pd.DataFrame, top_n: int = 500) -> dict:
    """Fisher: known myeloma genes among genes linked to top-N vs the rest."""
    mm_genes = set()
    mm_path = ROOT / "data/reference/mm_genes.tsv"
    if mm_path.exists():
        tsv = pd.read_csv(mm_path, sep="\t")
        mm_genes |= set(tsv["gene"].astype(str).str.upper())
    dep_path = ROOT / "data/cache/depmap/depmap_summary.parquet"
    if dep_path.exists():
        dep = pd.read_parquet(dep_path)
        sel = dep[dep["depmap_mm_selective"].fillna(False)]
        mm_genes |= set(sel["gene"].astype(str).str.upper())
    if not mm_genes:
        raise FileNotFoundError("no known-myeloma gene source found (mm_genes.tsv / DepMap)")

    df = df.sort_values("score", ascending=False)
    top = df.head(top_n)
    rest = df.iloc[top_n:]

    def genes_of(frame: pd.DataFrame) -> set[str]:
        out: set[str] = set()
        for _, row in frame.iterrows():
            out.update(g.upper() for g in _split_genes(row.get("target_genes")))
            ng = row.get("nearest_gene")
            if isinstance(ng, str) and ng:
                out.add(ng.upper())
        return out

    gt, gr = genes_of(top), genes_of(rest)
    a = len(gt & mm_genes)
    b = len(gt - mm_genes)
    c = len(gr & mm_genes)
    d = len(gr - mm_genes)
    odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return {
        "test": f"top{top_n}_linked_genes_vs_rest__known_myeloma_genes",
        "contingency": {"top_known": a, "top_other": b,
                        "rest_known": c, "rest_other": d},
        "known_gene_pool_size": len(mm_genes),
        "odds_ratio": float(odds),
        "p_value": float(p),
        "alternative": "greater",
        "independent_of_scoring": False,
        "note": "PARTIALLY CIRCULAR: myeloma-relevance (weight 0.25) uses "
                "known-gene/DepMap flags; gene LINKAGE itself (s4) does not.",
    }


def lead_recovery(df_b: pd.DataFrame, k: int = 3) -> dict:
    """Per locus: is any published lead SNP among the top-k scored candidates?"""
    if "gwas_role" not in df_b.columns:
        raise KeyError("Mode B table lacks gwas_role column")
    per_locus = []
    recovered = 0
    expected_ps = []
    for lid, sub in df_b.groupby("locus_id"):
        sub = sub.sort_values("score", ascending=False).reset_index(drop=True)
        n = len(sub)
        if n < 2:
            continue
        lead_rows = sub[sub["gwas_role"] == "lead"]
        if lead_rows.empty:
            continue
        lead_rank = int(lead_rows.index.min()) + 1     # 1-based best rank
        kk = min(k, n)
        hit = lead_rank <= kk
        recovered += int(hit)
        expected_ps.append(kk / n)
        per_locus.append({"locus_id": lid, "n_candidates": n,
                          "lead_rank": lead_rank, "recovered_at_k": bool(hit)})
    n_loci = len(per_locus)
    if n_loci == 0:
        raise RuntimeError("no Mode B loci with lead SNPs scored")
    pbar = float(np.mean(expected_ps))
    p = float(binomtest(recovered, n_loci, pbar, alternative="greater").pvalue)
    return {
        "test": f"mode_b_published_lead_in_top{k}",
        "n_loci": n_loci,
        "k": k,
        "recovered": recovered,
        "recovery_rate": recovered / n_loci,
        "expected_rate_random_ranking": pbar,
        "p_value": p,
        "independent_of_scoring": True,
        "note": "all candidates in a locus share the same GWAS flag, so the "
                "GWAS component cannot favour the lead SNP",
        "per_locus": per_locus,
    }


# ------------------------------------------------------------- figures -----
def fig_funnel(funnel: list[dict], path: Path) -> None:
    labels = [f"{s['step']} {s['label']}" for s in funnel]
    vals = [s["rows_out"] for s in funnel]
    fig, ax = plt.subplots(figsize=(9, 4.2))
    y = np.arange(len(vals))
    ax.barh(y, vals, color="#4da3ff")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xscale("symlog", linthresh=1)
    for i, v in enumerate(vals):
        ax.text(v * 1.15, i, f"{v:,}", va="center", fontsize=8)
    ax.set_xlabel("variants remaining (symlog)")
    ax.set_title("Mode A funnel — GIAB HG001 whole genome (real counts)")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def fig_scores(df: pd.DataFrame, path: Path, top_n: int = 1000) -> None:
    df = df.sort_values("score", ascending=False)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(df["score"], bins=80, color="#2f6fb4")
    cut = float(df["score"].iloc[min(top_n, len(df) - 1)])
    ax.axvline(cut, color="#ffb454", linestyle="--",
               label=f"top-{top_n} cut = {cut:.3f}")
    ax.set_xlabel("combined score")
    ax.set_ylabel("variants")
    ax.legend(fontsize=8)
    ax.set_title("Score distribution of ranked regulatory variants")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def fig_enrichment(tests: list[dict], path: Path) -> None:
    usable = [t for t in tests if t.get("p_value") is not None]
    names = [t["test"][:44] for t in usable]
    mlog = [-np.log10(max(float(t["p_value"]), 1e-300)) for t in usable]
    fig, ax = plt.subplots(figsize=(8, 0.7 + 0.8 * len(usable)))
    colors = ["#7ee0a3" if t.get("independent_of_scoring") else "#ffb454"
              for t in usable]
    ax.barh(np.arange(len(usable)), mlog, color=colors)
    ax.set_yticks(np.arange(len(usable)))
    ax.set_yticklabels(names, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("-log10(p), Fisher exact (one-sided greater)")
    for i, t in enumerate(usable):
        if t.get("odds_ratio") is not None:
            stat = f"OR={t['odds_ratio']:.2f}"
        else:
            stat = f"rate={t.get('recovery_rate', 0):.2f}"
        ax.text(mlog[i] + 0.1, i,
                f"{stat} p={t['p_value']:.2e}",
                va="center", fontsize=8)
    ax.set_title("Enrichment tests (green = independent of scoring, "
                 "amber = partially circular)")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def fig_lead_ranks(recovery: dict, path: Path) -> None:
    rows = recovery["per_locus"]
    fig, ax = plt.subplots(figsize=(9, 4))
    lids = [r["locus_id"] for r in rows]
    ranks = [r["lead_rank"] for r in rows]
    colors = ["#7ee0a3" if r["recovered_at_k"] else "#ff6b6b" for r in rows]
    ax.bar(np.arange(len(rows)), ranks, color=colors)
    ax.set_xticks(np.arange(len(rows)))
    ax.set_xticklabels(lids, rotation=90, fontsize=7)
    ax.axhline(recovery["k"], color="#ffb454", linestyle="--",
               label=f"top-{recovery['k']} cutoff")
    ax.set_ylabel("rank of published lead SNP (1 = best)")
    ax.legend(fontsize=8)
    ax.set_title("Mode B: lead-SNP rank within each locus")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# --------------------------------------------------------------- main ------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--top-n", type=int, default=1000)
    args = ap.parse_args()

    result_a_path = ROOT / "outputs/result.json"
    csv_a_path = ROOT / "outputs/ranked_variants.csv"
    result_b_path = ROOT / "outputs/mode_b/result.json"
    csv_b_path = ROOT / "outputs/mode_b/ranked_variants.csv"

    if not result_a_path.exists() or not csv_a_path.exists():
        raise SystemExit("Mode A outputs missing - run `make run` first")

    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    result_a = json.loads(result_a_path.read_text())
    df = pd.read_csv(csv_a_path, low_memory=False)

    tests: list[dict] = [ccre_enrichment(df, top_n=args.top_n),
                         gene_enrichment(df, top_n=500)]

    mode_b: dict | None = None
    if result_b_path.exists() and csv_b_path.exists():
        df_b = pd.read_csv(csv_b_path, low_memory=False)
        mode_b = lead_recovery(df_b, k=3)
        tests.append({k: v for k, v in mode_b.items() if k != "per_locus"})
        tests[-1]["contingency"] = None
        fig_lead_ranks(mode_b, FIG / "fig_mode_b_lead_ranks.png")

    fig_funnel(result_a["funnel"], FIG / "fig_funnel.png")
    fig_scores(df, FIG / "fig_score_hist.png", top_n=args.top_n)
    fig_enrichment(tests, FIG / "fig_enrichment.png")

    s3 = next((f for f in result_a["funnel"] if f["step"] == "s3"), {})
    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode_a": {
            "run_id": result_a["run"]["run_id"],
            "total_ranked": result_a["total_ranked"],
            "elapsed_seconds": result_a["run"]["elapsed_seconds"],
            "manifest_hash": result_a["run"]["manifest_hash"],
            "funnel": result_a["funnel"],
            "sensitivity": result_a["sensitivity"],
            "motif_note": s3.get("note", ""),
            "warnings": result_a.get("warnings", []),
        },
        "enrichment_tests": tests,
        "mode_b": mode_b,
        "figures": sorted(p.name for p in FIG.glob("*.png")),
    }
    out_path = OUT / "results.json"
    out_path.write_text(json.dumps(payload, indent=1))
    log.info("wrote %s (%d tests, %d figures)",
             out_path, len(tests), len(payload["figures"]))
    for t in tests:
        log.info("  %-60s OR=%s p=%.3g independent=%s",
                 t["test"], round(t.get("odds_ratio") or 0, 2),
                 float(t["p_value"]), t.get("independent_of_scoring"))
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    raise SystemExit(main())
