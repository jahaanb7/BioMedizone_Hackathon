"""Template-based one-sentence explanations built ONLY from real computed
values; any clause whose data is missing is omitted (never invented)."""
from __future__ import annotations

from typing import Any, Mapping


def _fmt_num(x: Any, digits: int = 2) -> str:
    try:
        return f"{float(x):.{digits}f}"
    except (TypeError, ValueError):
        return ""


def _active_element_clause(row: Mapping[str, Any]) -> str:
    tracks = row.get("active_tracks") or ""
    names = [t for t in str(tracks).split(",") if t]
    if not names:
        return ""
    pretty = ", ".join(sorted({n.split("_")[0] for n in names})[:4])
    cell = "MM.1S"
    clause = f"lies in {cell} active element(s) ({pretty})"
    if row.get("in_super_enhancer"):
        clause += " inside a predicted super-enhancer"
    return clause


def _motif_clause(row: Mapping[str, Any]) -> str:
    tf = row.get("top_motif_tf")
    eff = row.get("motif_effect") or "none"
    if not tf or eff == "none":
        return ""
    r = _fmt_num(row.get("motif_ref_rel"))
    a = _fmt_num(row.get("motif_alt_rel"))
    if eff == "broken":
        verb = "weakens a predicted"
    elif eff == "created":
        verb = "creates a predicted"
    else:
        return ""
    return f"It {verb} {tf} site (relative score {r} -> {a})"


def _gene_clause(row: Mapping[str, Any]) -> str:
    genes = [g for g in str(row.get("target_genes") or "").split(",") if g]
    if not genes:
        return ""
    gene = genes[0]
    method = row.get("link_method") or "none"
    conf = row.get("link_confidence") or "low"
    parts = [f"The element is linked to {gene}"]
    if method and method != "none":
        parts.append(f"({method}")
        if method.startswith("ABC"):
            score = row.get("link_score")
            if score:
                parts.append(f"score {_fmt_num(score)}")
        parts.append(f"{conf} confidence)")
    return " ".join(parts).replace(") ", ") ")


def _relevance_clauses(row: Mapping[str, Any]) -> list[str]:
    out: list[str] = []
    genes = [g for g in str(row.get("target_genes") or "").split(",") if g]
    gene = genes[0] if genes else (row.get("nearest_gene") or "")
    if not gene:
        return out
    if row.get("depmap_mm_selective"):
        out.append(f"{gene} is a selective dependency in myeloma cell lines (DepMap)")
    elif row.get("depmap_mm_mean_effect") is not None:
        try:
            eff = float(row.get("depmap_mm_mean_effect"))
            if eff <= -0.3:
                out.append(f"{gene} shows a dependency signal in myeloma lines "
                           f"(mean gene effect {_fmt_num(eff)})")
        except (TypeError, ValueError):
            pass
    if row.get("is_known_mm_gene"):
        out.append(f"{gene} is a known myeloma driver gene")
    if row.get("expressed_in_mm"):
        out.append(f"{gene} is expressed in myeloma cell lines")
    if row.get("gwas_flag"):
        out.append("this variant is a published myeloma GWAS lead/proxy variant")
    elif row.get("in_mm_gwas_locus"):
        out.append("this variant falls in a myeloma GWAS locus")
    return out


def build_explanation(row: Mapping[str, Any], rank: int, total: int) -> str:
    """Compose the evidence-chain sentence(s) for one ranked variant."""
    vid = str(row.get("variant_id", "?"))
    ref, alt = str(row.get("ref", "?")), str(row.get("alt", "?"))
    clauses: list[str] = []

    head = f"{vid} {ref}>{alt}"
    elem = _active_element_clause(row)
    clauses.append(f"{head} {elem}" if elem else f"{head}")

    motif = _motif_clause(row)
    if motif:
        clauses.append(motif)

    gene = _gene_clause(row)
    if gene:
        clauses.append(gene)

    clauses.extend(_relevance_clauses(row))

    text = ". ".join(clauses)
    if not text.endswith("."):
        text += "."
    text += f" Rank {rank} of {total}."
    return text


def row_warnings(row: Mapping[str, Any]) -> list[str]:
    """Per-variant caveat list (missing/low-confidence evidence only)."""
    out: list[str] = []
    if (row.get("link_confidence") == "low"):
        method = row.get("link_method") or "nearest-gene"
        out.append(f"gene link is {method}-based, low confidence")
    if not row.get("top_motif_tf") or (row.get("motif_effect") == "none"):
        out.append("no strong motif disruption predicted")
    if row.get("af") is None and not row.get("in_common_track"):
        out.append("no population frequency found in the common-variant table (treated as rare)")
    if row.get("depmap_mm_mean_effect") is None:
        out.append("no DepMap gene-effect for the linked gene")
    return out
