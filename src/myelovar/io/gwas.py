"""Mode B (GWAS loci) and Mode C (DepMap somatic) input builders.

Mode B flow (all GRCh38, all API-cached on disk):
  GWAS Catalog v2 leads -> locus windows -> Ensembl region variations
  -> EUR LD proxies (r2 >= threshold, 1000G phase 3) -> candidate table.
Common-variant filtering is disabled for Mode B (GWAS signals are common).
"""
from __future__ import annotations

import json
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

from myelovar.config import AppConfig
from myelovar.io import ensembl

log = logging.getLogger(__name__)

GWAS_CACHE = Path("data/cache/gwas")


def load_gwas_payload(cache_dir: Path | str = GWAS_CACHE) -> dict:
    path = Path(cache_dir) / "myeloma_gwas.json"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} missing - run `make data` (GWAS Catalog v2 fetch failed?)")
    return json.loads(path.read_text())


def build_loci(config: AppConfig, cache_dir: Path | str = GWAS_CACHE,
               force: bool = False) -> list[dict]:
    """Merge lead-SNP +/-window intervals into loci; cached to loci.json."""
    cache_dir = Path(cache_dir)
    out = cache_dir / "loci.json"
    if out.exists() and not force:
        return json.loads(out.read_text())
    payload = load_gwas_payload(cache_dir)
    window = config.gwas.window_kb * 1000
    leads: list[dict] = []
    for a in payload["gws_leads"]:
        loc = a.get("locations") or []
        if not loc:
            continue
        chrom_s, pos_s = str(loc[0]).rsplit(":", 1)
        try:
            pos = int(pos_s)
        except ValueError:
            continue
        snp = (a.get("snp_allele") or [{}])[0]
        leads.append({
            "rsid": snp.get("rs_id", ""),
            "chrom": chrom_s if chrom_s.startswith("chr") else f"chr{chrom_s}",
            "pos": pos,
            "p_value": a.get("p_value"),
            "effect_allele": snp.get("effect_allele"),
            "reported_genes": a.get("mapped_genes") or [],
            "trait": a.get("_disease_trait"),
            "study": a.get("_study"),
            "pmid": a.get("pubmed_id"),
        })
    if not leads:
        raise RuntimeError("no usable lead SNPs in GWAS payload")
    # merge overlapping +/-window spans per chromosome
    spans: dict[str, list[tuple[int, int]]] = {}
    for L in leads:
        spans.setdefault(L["chrom"], []).append((L["pos"] - window, L["pos"] + window))
    merged: dict[str, list[tuple[int, int, int]]] = {}
    for chrom, sp in spans.items():
        sp.sort()
        out_sp = [list(sp[0])]
        for s, e in sp[1:]:
            if s <= out_sp[-1][1]:
                out_sp[-1][1] = max(out_sp[-1][1], e)
            else:
                out_sp.append([s, e])
        merged[chrom] = [(int(s), int(e), i) for i, (s, e) in enumerate(out_sp)]
    loci: list[dict] = []
    counter = 0
    for chrom, sp in merged.items():
        for s, e, _i in sp:
            counter += 1
            member = [L for L in leads if L["chrom"] == chrom and s <= L["pos"] <= e]
            loci.append({"locus_id": f"L{counter}", "chrom": chrom,
                         "start": max(0, s), "end": e, "leads": member})
    loci.sort(key=lambda x: x["locus_id"])
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(loci, indent=1))
    log.info("built %d loci from %d leads (window +/-%d kb)",
             len(loci), len(leads), config.gwas.window_kb)
    return loci


def _variation_rows(detail: dict, chrom_hint: str) -> list[dict]:
    """Rows from an Ensembl variation record (GRCh38 mapping, split by ALT)."""
    rows = []
    rsid = detail.get("name", "")
    for m in detail.get("mappings", []):
        if m.get("assembly_name") != "GRCh38":
            continue
        chrom = m.get("seq_region_name", "")
        chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
        pos = int(m.get("start", 0))
        allele_string = m.get("allele_string", "")
        # Ensembl uses "A/C/G/T" (variation endpoint) or "A|C" style strings
        alleles = re.split(r"[/|]", allele_string) if allele_string else []
        if not alleles or not pos:
            continue
        ref, alts = alleles[0], alleles[1:]
        if len(ref) > 50 or any(len(a) > 50 for a in alts):
            continue
        for alt in alts:
            if alt in ("<del>", "*", "") or not alt:
                continue
            rows.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                         "rsid": rsid, "_detail": detail})
        break  # one GRCh38 mapping per variation is enough
    return rows


def build_mode_b_variants(config: AppConfig, force: bool = False) -> pd.DataFrame:
    """Build the Mode B candidate table (locus windows + LD proxies)."""
    cache = GWAS_CACHE / "mode_b_variants.parquet"
    if cache.exists() and not force:
        df = pd.read_parquet(cache)
        log.info("Mode B variant table cache hit: %d rows", len(df))
        return df

    loci = build_loci(config, force=force)
    window = config.gwas.window_kb * 1000
    r2_thr = config.gwas.ld_r2
    pop = config.gwas.population
    # r2-sorted, so the best proxies come first; 20/lead keeps the Ensembl
    # detail traffic polite (~1.6k requests) while covering each locus
    max_proxies_per_lead = 20

    rows: list[dict] = []
    by_key: dict[tuple, int] = {}

    def add(rec: dict, role: str, r2=None, locus_id="", lead_rsid="") -> None:
        key = (rec["chrom"], rec["pos"], rec["ref"], rec["alt"])
        if key in by_key:
            # The variant already exists as a window row (leads fall inside
            # their own window by construction) - upgrade the role instead of
            # silently dropping the lead/proxy annotation.
            idx = by_key[key]
            if role in ("lead", "proxy") and rows[idx]["gwas_role"] == "window":
                rows[idx].update(gwas_flag=True, gwas_role=role, gwas_r2=r2,
                                 lead_rsid=lead_rsid or rows[idx]["lead_rsid"],
                                 locus_id=locus_id or rows[idx]["locus_id"])
                detail = rec.get("_detail")
                if detail and rows[idx].get("af") is None:
                    rows[idx]["af"] = ensembl.eur_af(detail, pop,
                                                     alt=rec["alt"])
            return
        detail = rec.pop("_detail", {})
        by_key[key] = len(rows)
        rows.append({**rec, "gwas_flag": role in ("lead", "proxy"),
                     "gwas_role": role, "gwas_r2": r2, "locus_id": locus_id,
                     "lead_rsid": lead_rsid,
                     "af": ensembl.eur_af(detail, pop, alt=rec["alt"])
                     if detail else None})

    for locus in loci:
        lid, chrom = locus["locus_id"], locus["chrom"]
        # --- window variants via Ensembl region overlap (GRCh38)
        feats = ensembl.region_variations(chrom, locus["start"], locus["end"])
        n_before = len(rows)
        for f in feats:
            alleles = f.get("alleles") or []
            if not alleles or not f.get("start"):
                continue
            ref = alleles[0]
            if len(ref) > 50:
                continue
            for alt in alleles[1:]:
                if not alt or alt in ("<del>", "*"):
                    continue
                if len(alt) > 50:
                    continue
                rec = {"chrom": chrom, "pos": int(f["start"]), "ref": ref,
                       "alt": alt, "rsid": f.get("id", ""),
                       "consequence_hint": f.get("consequence_type", "")}
                key = (rec["chrom"], rec["pos"], rec["ref"], rec["alt"])
                if key not in by_key:
                    by_key[key] = len(rows)
                    rows.append({**rec, "gwas_flag": False, "gwas_role": "window",
                                 "gwas_r2": None, "locus_id": lid,
                                 "lead_rsid": "", "af": None})
        log.info("locus %s: %d window variants", lid, len(rows) - n_before)

        # --- leads + LD proxies
        for lead in locus["leads"]:
            rsid = lead["rsid"]
            if not rsid:
                continue
            try:
                detail = ensembl.variation_detail(rsid, with_pops=True)
            except RuntimeError as exc:
                log.warning("variation detail failed for %s: %s", rsid, exc)
                detail = {}
            for rec in _variation_rows(detail, chrom):
                add(rec, "lead", r2=1.0, locus_id=lid, lead_rsid=rsid)
            proxies = ensembl.ld_proxies(rsid, pop, r2_thr)
            pairs = []
            for p in proxies:
                v1 = p.get("variation1") or p.get("rs_id1")
                v2 = p.get("variation2") or p.get("rs_id2")
                other = v2 if v1 == rsid else v1
                if other and other != rsid:
                    try:
                        pairs.append((float(p.get("r2", 0)), other))
                    except (TypeError, ValueError):
                        continue
            pairs.sort(reverse=True)
            for r2, other in pairs[:max_proxies_per_lead]:
                try:
                    pd_ = ensembl.variation_detail(other, with_pops=True)
                except RuntimeError:
                    continue
                for rec in _variation_rows(pd_, chrom):
                    add(rec, "proxy", r2=r2, locus_id=lid, lead_rsid=rsid)

    if not rows:
        raise RuntimeError("Mode B builder produced no variants")
    df = pd.DataFrame(rows)
    df["variant_id"] = (df["chrom"] + ":" + df["pos"].astype(str) + ":"
                        + df["ref"] + ":" + df["alt"])
    df = df.drop_duplicates(subset="variant_id").reset_index(drop=True)
    df["filter"] = "PASS"
    df.to_parquet(cache, index=False)
    log.info("Mode B table: %d variants across %d loci "
             "(leads=%d, proxies=%d, window=%d)",
             len(df), df["locus_id"].nunique(),
             int((df["gwas_role"] == "lead").sum()),
             int((df["gwas_role"] == "proxy").sum()),
             int((df["gwas_role"] == "window").sum()))
    return df


def build_mode_c_table(config: AppConfig) -> pd.DataFrame:
    """Mode C: DepMap somatic mutations in myeloma lines as a variant table."""
    from myelovar.io.depmap import find_depmap_files, myeloma_model_ids, somatic_mutation_lines
    files = find_depmap_files(config.paths.cache_dir)
    mut = somatic_mutation_lines(files)
    if mut is None:
        raise FileNotFoundError(
            "Mode C requires OmicsSomaticMutations.csv - rerun `make data -- --with-mode-c`")
    cols = {c.lower(): c for c in mut.columns}
    def pick(*names):
        for n in names:
            if n in cols:
                return cols[n]
        raise RuntimeError(f"somatic file missing column {names}; have {list(cols)[:20]}")
    chrom_c = pick("chromosome", "chrom", "#chromosome")
    pos_c = pick("start_position", "start_position_1based", "pos", "start")
    ref_c = pick("reference_allele", "ref", "reference_allele_1based")
    alt_c = pick("tumor_seq_allele2", "alt", "tumor_allele")
    mut = mut.copy()
    mut["chrom"] = mut[chrom_c].astype(str).apply(
        lambda c: c if c.startswith("chr") else f"chr{c}")
    df = pd.DataFrame({
        "chrom": mut["chrom"],
        "pos": pd.to_numeric(mut[pos_c], errors="coerce"),
        "ref": mut[ref_c].astype(str).str.upper(),
        "alt": mut[alt_c].astype(str).str.upper(),
    }).dropna()
    df["pos"] = df["pos"].astype(np.int64)
    df = df[df["chrom"].str.match(r"^chr([1-9]|1[0-9]|2[0-2]|X)$")]
    df["variant_id"] = (df["chrom"] + ":" + df["pos"].astype(str) + ":"
                        + df["ref"] + ":" + df["alt"])
    df = df.drop_duplicates(subset="variant_id").reset_index(drop=True)
    df["filter"] = "PASS"
    log.info("Mode C table: %d somatic variants in myeloma lines", len(df))
    return df
