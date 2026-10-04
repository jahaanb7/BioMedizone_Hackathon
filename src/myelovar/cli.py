"""Typer CLI: `myelovar run | data | cache | serve | export | demo`."""
from __future__ import annotations

import json
import logging
import subprocess
import sys
import webbrowser
from pathlib import Path
from typing import Optional

import typer

from myelovar import __version__
from myelovar.config import load_config
from myelovar.pipeline import load_result, run_pipeline

app = typer.Typer(help="Myelovar: myeloma regulatory variant prioritization "
                       "(research use only).", no_args_is_help=True)
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(name)s: %(message)s")


def _progress(p: dict) -> None:
    typer.echo(f"  [{int(p['fraction'] * 100):3d}%] {p['step']}: {p['message']}")


@app.command()
def run(mode: str = typer.Option("A", help="A: genome funnel, B: GWAS loci, C: DepMap somatic"),
        input: Optional[Path] = typer.Option(None, help="input VCF (.vcf/.vcf.gz), Mode A"),
        config: Optional[Path] = typer.Option(None, help="config.yaml path"),
        max_records: Optional[int] = typer.Option(None, help="cap input records (testing)"),
        out: Optional[Path] = typer.Option(None, help="outputs directory override")) -> None:
    """Run the full pipeline and write outputs/."""
    cfg = load_config(config)
    if out is not None:
        cfg.paths.outputs_dir = str(out)
    result = run_pipeline(str(input) if input else None, mode=mode, config=cfg,
                          progress_callback=_progress, max_records=max_records)
    typer.echo(f"\nrun {result.run.run_id}: {result.run.input_description}")
    typer.echo(f"  funnel: " + " -> ".join(
        f"{s.rows_out:,}" for s in result.funnel))
    typer.echo(f"  ranked CSV: {result.outputs['ranked_csv']}")
    typer.echo(f"  result:     {result.outputs['result_json']}")
    if result.warnings:
        typer.echo(f"  warnings ({len(result.warnings)}):")
        for wmsg in result.warnings[:10]:
            typer.echo(f"    - {wmsg}")


@app.command()
def data(fast: bool = typer.Option(False, help="skip slow optional downloads"),
         with_mode_c: bool = typer.Option(False, help="also fetch DepMap somatic")) -> None:
    """Download + verify all public datasets (writes data/data_manifest.json)."""
    args = [sys.executable, "scripts/download_data.py"]
    if fast:
        args.append("--fast")
    if with_mode_c:
        args.append("--with-mode-c")
    raise typer.Exit(subprocess.call(args))


@app.command()
def cache() -> None:
    """Build precomputed caches (AF extract, ROSE super-enhancers, DepMap/ABC)."""
    raise typer.Exit(subprocess.call([sys.executable, "scripts/build_cache.py"]))


@app.command()
def serve(host: str = typer.Option("127.0.0.1"),
          port: int = typer.Option(8000),
          open_browser: bool = typer.Option(False, "--open")) -> None:
    """Start the REST API (+ web UI) with uvicorn."""
    import uvicorn
    url = f"http://{host}:{port}/"
    if open_browser:
        webbrowser.open(url)
    uvicorn.run("api.main:app", host=host, port=port, log_level="info")


@app.command()
def export(result: Path = typer.Option("outputs/result.json"),
           csv: Path = typer.Option("outputs/exported_variants.csv")) -> None:
    """Export a stored result.json to CSV."""
    import pandas as pd
    r = load_result(result)
    rows = [v.model_dump() for v in r.variants]
    pd.DataFrame(rows).to_csv(csv, index=False)
    typer.echo(f"wrote {len(rows)} rows -> {csv}")


@app.command()
def version() -> None:
    """Print the package version."""
    typer.echo(__version__)


if __name__ == "__main__":
    app()
