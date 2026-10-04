"""Build the offline demo bundle from real pipeline outputs.

Produces:
  demo_bundle/result_{A,B}.json        full PipelineResult (capped variant list)
  demo_bundle/mode_{a,b}/ranked_*.csv  full ranked table (served by /variants)
  demo_bundle/loci/<variant>.json      precomputed locus payloads for the top
                                       N variants of each mode (fully offline)

The API serves these for run ids demo-A / demo-B without touching reference
data, so `make demo` works on a machine with no downloads.
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

log = logging.getLogger("demo_bundle")


def _copy_mode(mode: str, result_src: Path, csv_src: Path, demo_dir: Path) -> bool:
    if not (result_src.exists() and csv_src.exists()):
        log.warning("mode %s: outputs missing (%s, %s) - skipping",
                    mode, result_src, csv_src)
        return False
    dest = demo_dir / f"mode_{mode.lower()}"
    dest.mkdir(parents=True, exist_ok=True)
    shutil.copy2(result_src, demo_dir / f"result_{mode}.json")
    shutil.copy2(result_src, dest / "result.json")
    shutil.copy2(csv_src, dest / "ranked_variants.csv")
    log.info("mode %s: result + CSV copied (%.1f MB CSV)",
             mode, csv_src.stat().st_size / 1e6)
    return True


def _precompute_loci(mode: str, csv_src: Path, demo_dir: Path, top_n: int) -> int:
    from api.routes import _build_locus_payload, _find_variant

    df = pd.read_csv(csv_src, low_memory=False)
    df = df.sort_values("score", ascending=False).head(top_n)
    loci_dir = demo_dir / "loci"
    loci_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    for vid in df["variant_id"].tolist():
        rec = _find_variant(df, vid)
        payload = _build_locus_payload(rec, vid, 50_000)
        safe = re.sub(r"[^A-Za-z0-9]+", "_", vid)
        (loci_dir / f"{safe}.json").write_text(
            json.dumps(json.loads(payload.model_dump_json()), indent=1))
        n += 1
    log.info("mode %s: %d locus payloads precomputed", mode, n)
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--top-loci", type=int, default=30,
                    help="locus payloads per mode for the top-N variants")
    args = ap.parse_args()

    demo_dir = ROOT / "demo_bundle"
    demo_dir.mkdir(parents=True, exist_ok=True)

    built = []
    plan = [
        ("A", ROOT / "outputs/result.json", ROOT / "outputs/ranked_variants.csv"),
        ("B", ROOT / "outputs/mode_b/result.json",
         ROOT / "outputs/mode_b/ranked_variants.csv"),
    ]
    for mode, result_src, csv_src in plan:
        if _copy_mode(mode, result_src, csv_src, demo_dir):
            _precompute_loci(mode, csv_src, demo_dir, args.top_loci)
            built.append(mode)

    if not built:
        log.error("no outputs found - run `make run` (and `make run-b`) first")
        return 1
    log.info("demo bundle ready for modes: %s", ", ".join(built))
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    raise SystemExit(main())
