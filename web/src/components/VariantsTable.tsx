/**
 * The ranked list — the product's primary screen.
 *
 * Server-side paging/sorting/filtering against GET /api/runs/{id}/variants
 * (limit ≤ 500). Filter + sort state lives in the URL query string, so any
 * filtered view is bookmarkable and shareable. Listens for the
 * "myelovar:filter" window event so the guided tour can drive it.
 */
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import {
  useReactTable,
  getCoreRowModel,
  getSortedRowModel,
  flexRender,
  type SortingState,
} from "@tanstack/react-table";
import { useSearchParams } from "react-router-dom";
import toast from "react-hot-toast";
import { Copy, RotateCcw, Star } from "lucide-react";
import { api, type VariantsPage } from "../api/client";
import {
  CountUp,
  EmptyState,
  EvidenceDots,
  Gauge,
  Help,
  MotifChip,
  ScoreBar,
  SkeletonRows,
} from "./ui";

const PAGE = 50;

type Row = Record<string, unknown> & {
  variant_id: string;
  score: number;
  rank: number;
  chrom: string;
  pos: number;
};

type Counts = { total: number; broken: number; high: number; se: number };

const FILTER_KEYS = ["gene", "tf", "chrom", "confidence", "motif_effect", "in_super_enhancer", "min_score"];

const FILTER_LABELS: Record<string, string> = {
  gene: "gene",
  tf: "TF",
  chrom: "chr",
  confidence: "link",
  motif_effect: "motif",
  in_super_enhancer: "super-enhancer",
  min_score: "min score",
};

export default function VariantsTable({
  runId,
  mode,
  onSelect,
}: {
  runId: string;
  mode: string;
  onSelect: (variantId: string) => void;
}) {
  const [sp, setSp] = useSearchParams();
  const [data, setData] = useState<VariantsPage | null>(null);
  const [counts, setCounts] = useState<Counts | null>(null);
  const [offset, setOffset] = useState(0);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [nonce, setNonce] = useState(0);
  const geneRef = useRef<HTMLInputElement>(null);

  /* ---- URL-derived state ------------------------------------------- */
  const gene = sp.get("gene") ?? "";
  const tf = sp.get("tf") ?? "";
  const chrom = sp.get("chrom") ?? "";
  const confidence = sp.get("confidence") ?? "";
  const motif = sp.get("motif_effect") ?? "";
  const se = sp.get("in_super_enhancer") ?? "";
  const minScore = sp.get("min_score") ?? "";
  const sortId = sp.get("sort") ?? "score";
  const order = (sp.get("order") ?? "desc") as "asc" | "desc";

  const activeFilters = FILTER_KEYS.filter((k) => sp.get(k));

  const setParam = useCallback(
    (key: string, value: string) => {
      const next = new URLSearchParams(sp);
      if (value) next.set(key, value);
      else next.delete(key);
      setSp(next, { replace: true });
      setOffset(0);
    },
    [sp, setSp],
  );

  const reset = useCallback(() => {
    const next = new URLSearchParams(sp);
    for (const k of FILTER_KEYS) next.delete(k);
    setSp(next, { replace: true });
    setOffset(0);
  }, [sp, setSp]);

  const sorting: SortingState = useMemo(
    () => [{ id: sortId, desc: order === "desc" }],
    [sortId, order],
  );

  const onSortingChange: (updater: SortingState | ((old: SortingState) => SortingState)) => void =
    useCallback(
      (updater) => {
        const next = typeof updater === "function" ? updater(sorting) : updater;
        const s = next[0];
        setParam("sort", s && s.id ? s.id : "score");
        setParam("order", s && s.id && !s.desc ? "asc" : "desc");
      },
      [sorting, setParam],
    );

  /* ---- query --------------------------------------------------------- */
  const query = useMemo(
    () => ({
      limit: PAGE,
      offset,
      sort: sortId,
      order,
      gene: gene || undefined,
      tf: tf || undefined,
      chrom: chrom || undefined,
      confidence: confidence || undefined,
      motif_effect: motif || undefined,
      in_super_enhancer: se || undefined,
      min_score: minScore || undefined,
      locus_id: mode === "B" ? (sp.get("locus_id") ?? undefined) : undefined,
    }),
    [offset, sortId, order, gene, tf, chrom, confidence, motif, se, minScore, mode, sp],
  );

  useEffect(() => {
    let cancelled = false;
    setLoading(true);
    setError(null);
    api
      .variants(runId, query)
      .then((p) => !cancelled && setData(p))
      .catch((e) => !cancelled && setError(String(e.message ?? e)))
      .finally(() => !cancelled && setLoading(false));
    return () => {
      cancelled = true;
    };
  }, [runId, query, nonce]);

  /* ---- unfiltered stats for the strip -------------------------------- */
  useEffect(() => {
    let dead = false;
    setCounts(null);
    const empty = { limit: 1 };
    Promise.all([
      api.variants(runId, empty),
      api.variants(runId, { ...empty, motif_effect: "broken" }),
      api.variants(runId, { ...empty, confidence: "high" }),
      api.variants(runId, { ...empty, in_super_enhancer: "true" }),
    ])
      .then(([t, b, h, s]) => {
        if (!dead)
          setCounts({ total: t.total, broken: b.total, high: h.total, se: s.total });
      })
      .catch(() => undefined);
    return () => {
      dead = true;
    };
  }, [runId]);

  /* ---- tour / external filter events --------------------------------- */
  useEffect(() => {
    const handler = (e: Event) => {
      const detail = (e as CustomEvent).detail as Record<string, string>;
      const next = new URLSearchParams(sp);
      for (const [k, v] of Object.entries(detail)) {
        if (v) next.set(k, v);
        else next.delete(k);
      }
      setSp(next, { replace: true });
      setOffset(0);
      toast("Filtered by the tour", { icon: "✨", duration: 1600 });
    };
    window.addEventListener("myelovar:filter", handler);
    return () => window.removeEventListener("myelovar:filter", handler);
  }, [sp, setSp]);

  /* ---- keyboard shortcuts -------------------------------------------- */
  useEffect(() => {
    const handler = (e: KeyboardEvent) => {
      const tag = (e.target as HTMLElement | null)?.tagName ?? "";
      const typing = /INPUT|SELECT|TEXTAREA/.test(tag);
      if (e.key === "/" && !typing) {
        e.preventDefault();
        geneRef.current?.focus();
        geneRef.current?.select();
      }
      if (e.key === "Escape" && FILTER_KEYS.some((k) => sp.get(k))) reset();
    };
    window.addEventListener("keydown", handler);
    return () => window.removeEventListener("keydown", handler);
  }, [sp, reset]);

  /* ---- columns -------------------------------------------------------- */
  const columns = useMemo(
    () => [
      {
        header: "Rank",
        accessorKey: "rank",
        enableSorting: true,
        cell: (c: { getValue: () => number }) => (
          <span className="cell-rank">{c.getValue()}</span>
        ),
      },
      {
        header: "Variant",
        accessorKey: "variant_id",
        enableSorting: true,
        cell: (c: { getValue: () => string }) => {
          const vid = c.getValue() as string;
          return (
            <span
              className="mono"
              onClick={(e) => {
                e.stopPropagation();
                void navigator.clipboard?.writeText(vid);
                toast.success("Variant ID copied", { duration: 1400 });
              }}
              title="click to copy"
              style={{ cursor: "copy" }}
            >
              {vid} <Copy size={10} style={{ opacity: 0.45, verticalAlign: -1 }} />
            </span>
          );
        },
      },
      {
        header: () => (
          <span>
            Score
            <Help text="Weighted sum: regulatory context 0.25 · motif disruption 0.20 · gene link 0.20 · myeloma relevance 0.25 · GWAS 0.10. Deterministic and reproducible." />
          </span>
        ),
        accessorKey: "score",
        enableSorting: true,
        cell: (c: { getValue: () => number }) => <ScoreBar value={c.getValue() as number} />,
      },
      {
        header: "Motif effect",
        accessorKey: "motif_effect",
        enableSorting: false,
        cell: (c: { getValue: () => string }) => <MotifChip effect={c.getValue() as string} />,
      },
      {
        header: "Top TF",
        accessorKey: "top_motif_tf",
        enableSorting: false,
        cell: (c: { getValue: () => unknown }) => (
          <span className="mono">{(c.getValue() as string) ?? "—"}</span>
        ),
      },
      {
        header: "Target genes",
        accessorKey: "target_genes",
        enableSorting: false,
        cell: (c: { getValue: () => unknown; row: { original: Row } }) => {
          const raw = c.getValue();
          const list =
            typeof raw === "string" && raw
              ? raw.split(",")
              : Array.isArray(raw)
                ? raw.map(String)
                : [];
          const first = list[0] ?? (c.row.original.nearest_gene as string) ?? "—";
          const conf = (c.row.original.link_confidence as string) ?? "low";
          return (
            <span className="gene-cell">
              {first}
              {list.length > 1 && <span className="plus">+{list.length - 1}</span>}{" "}
              <span className={`pill ${conf}`} title={`${conf} link confidence`}>
                {conf}
              </span>
            </span>
          );
        },
      },
      {
        header: "Context",
        accessorKey: "regulatory_context_score",
        enableSorting: false,
        cell: (c: { getValue: () => number }) => <Gauge value={(c.getValue() as number) ?? 0} />,
      },
      {
        header: () => (
          <span>
            Evidence
            <Help text="Dots: DepMap myeloma-selective dependency · known myeloma gene · inside a myeloma GWAS locus." />
          </span>
        ),
        id: "evidence",
        enableSorting: false,
        cell: (c: { row: { original: Row } }) => {
          const o = c.row.original;
          return (
            <EvidenceDots
              depmap={Boolean(o.depmap_mm_selective)}
              known={Boolean(o.is_known_mm_gene)}
              gwas={Boolean(o.in_mm_gwas_locus ?? o.gwas_flag)}
            />
          );
        },
      },
      {
        header: "Super-enh.",
        accessorKey: "in_super_enhancer",
        enableSorting: false,
        cell: (c: { getValue: () => unknown }) =>
          c.getValue() ? (
            <Star size={13} className="star" fill="currentColor" />
          ) : (
            <span className="muted">—</span>
          ),
      },
    ],
    [],
  );

  const table = useReactTable({
    columns,
    data: (data?.rows ?? []) as unknown as Row[],
    state: { sorting },
    onSortingChange,
    getCoreRowModel: getCoreRowModel(),
    getSortedRowModel: getSortedRowModel(),
    manualSorting: true,
  });

  const rows = data?.rows ?? [];
  const hasAny = rows.length > 0;

  /* ---- render --------------------------------------------------------- */
  const statChip = (
    key: keyof Counts | "total",
    label: string,
    tip: string,
    on: boolean,
    apply: () => void,
  ) => (
    <button
      className={`stat-chip ${on ? "on" : ""}`}
      onClick={apply}
      title={tip}
      type="button"
    >
      <span className="n">{counts ? <CountUp value={counts[key]} /> : "…"}</span>
      <span className="l">
        {label}
        <Help text={tip} />
      </span>
    </button>
  );

  return (
    <div className="card">
      <div className="stats-strip" data-tour="stats-strip">
        {statChip("total", "Total ranked", "All ranked variants in this run, unfiltered.", activeFilters.length === 0, reset)}
        {statChip(
          "broken",
          "Break TF motifs",
          "Variants that destroy a transcription-factor binding site (motif effect = broken).",
          motif === "broken",
          () => setParam("motif_effect", motif === "broken" ? "" : "broken"),
        )}
        {statChip(
          "high",
          "High-confidence links",
          "Variants whose target gene comes from the ABC contact model (high link confidence).",
          confidence === "high",
          () => setParam("confidence", confidence === "high" ? "" : "high"),
        )}
        {statChip(
          "se",
          "In super-enhancers",
          "Variants inside a ROSE super-enhancer of MM.1S myeloma cells.",
          se === "true",
          () => setParam("in_super_enhancer", se === "true" ? "" : "true"),
        )}
      </div>

      <div className="filters">
        <div>
          <label>Gene</label>
          <input
            ref={geneRef}
            value={gene}
            onChange={(e) => setParam("gene", e.target.value)}
            placeholder="IRF4  (press /)"
          />
        </div>
        <div>
          <label>TF</label>
          <input value={tf} onChange={(e) => setParam("tf", e.target.value)} placeholder="MYC" />
        </div>
        <div>
          <label>Chrom</label>
          <input
            value={chrom}
            onChange={(e) => setParam("chrom", e.target.value)}
            placeholder="chr1"
            style={{ width: 74 }}
          />
        </div>
        <div data-tour="filter-motif">
          <label>
            Motif effect
            <Help text="broken = strong site destroyed · created = new strong site formed · changed = score shifted below the thresholds." />
          </label>
          <select value={motif} onChange={(e) => setParam("motif_effect", e.target.value)}>
            <option value="">any</option>
            <option value="broken">broken</option>
            <option value="created">created</option>
            <option value="changed">changed</option>
            <option value="none">none</option>
          </select>
        </div>
        <div>
          <label>Link confidence</label>
          <select value={confidence} onChange={(e) => setParam("confidence", e.target.value)}>
            <option value="">any</option>
            <option value="high">high (ABC)</option>
            <option value="medium">medium (distance)</option>
            <option value="low">low (nearest)</option>
          </select>
        </div>
        <div>
          <label>
            Min score
            <Help text="Minimum composite score (0–1) for display." />
          </label>
          <input
            value={minScore}
            onChange={(e) => setParam("min_score", e.target.value)}
            placeholder="0.3"
            style={{ width: 80 }}
          />
        </div>
      </div>

      <div className="filter-row">
        {activeFilters.map((k) => (
          <span className="chip-filter" key={k}>
            {FILTER_LABELS[k]}: {sp.get(k)}
            <button aria-label={`clear ${k}`} onClick={() => setParam(k, "")}>
              ×
            </button>
          </span>
        ))}
        <span className="match-count">
          {loading && !data ? (
            "matching…"
          ) : (
            <>
              <b>{(data?.total ?? 0).toLocaleString()}</b> of{" "}
              {(counts?.total ?? data?.total ?? 0).toLocaleString()} variants match
            </>
          )}
        </span>
        {activeFilters.length > 0 && (
          <button className="ghost small" onClick={reset}>
            <RotateCcw size={12} /> Reset
          </button>
        )}
      </div>

      {error && (
        <div className="empty-state">
          <h4 className="error">Could not load variants</h4>
          <p>{error}</p>
          <button className="secondary" onClick={() => setNonce((n) => n + 1)}>
            Retry
          </button>
        </div>
      )}

      {!error && (
        <div className="table-wrap">
          <table className="var-table">
            <thead>
              {table.getHeaderGroups().map((hg) => (
                <tr key={hg.id}>
                  {hg.headers.map((h) => {
                    const sorted = h.column.getIsSorted();
                    return (
                      <th
                        key={h.id}
                        className={h.column.getCanSort() ? "sortable" : undefined}
                        onClick={h.column.getToggleSortingHandler()}
                        aria-sort={
                          sorted === "asc" ? "ascending" : sorted === "desc" ? "descending" : undefined
                        }
                      >
                        {flexRender(h.column.columnDef.header, h.getContext())}
                        {sorted === "asc" ? " ↑" : sorted === "desc" ? " ↓" : ""}
                      </th>
                    );
                  })}
                </tr>
              ))}
            </thead>
            <tbody style={{ opacity: loading && data ? 0.55 : 1 }}>
              {!data && loading && <SkeletonRows rows={10} cols={columns.length} />}
              {hasAny &&
                table.getRowModel().rows.map((row, i) => {
                  const r = row.original as Row;
                  return (
                    <tr
                      key={r.variant_id}
                      tabIndex={0}
                      className={(r.rank as number) <= 3 ? "top" : undefined}
                      data-tour={i === 0 ? "first-row" : undefined}
                      data-vid={r.variant_id}
                      onClick={() => onSelect(r.variant_id)}
                      onKeyDown={(e) => {
                        if (e.key === "Enter") onSelect(r.variant_id);
                      }}
                    >
                      {row.getVisibleCells().map((cell) => (
                        <td key={cell.id}>
                          {flexRender(cell.column.columnDef.cell, cell.getContext())}
                        </td>
                      ))}
                    </tr>
                  );
                })}
              {data && !hasAny && !loading && (
                <tr>
                  <td colSpan={columns.length} style={{ padding: 0, border: "none" }}>
                    <EmptyState
                      title="No variants match these filters"
                      body="Loosen or clear the filters to see the ranked list again."
                      action={reset}
                    />
                  </td>
                </tr>
              )}
            </tbody>
          </table>
        </div>
      )}

      {!error && data && (
        <div className="pager">
          <button
            className="secondary small"
            disabled={offset === 0}
            onClick={() => setOffset(Math.max(0, offset - PAGE))}
          >
            ← Prev
          </button>
          <span className="muted">
            {data.total ? (data.offset + 1).toLocaleString() : 0}–
            {Math.min(data.offset + PAGE, data.total).toLocaleString()} of{" "}
            {data.total.toLocaleString()}
          </span>
          <button
            className="secondary small"
            disabled={data.offset + PAGE >= data.total}
            onClick={() => setOffset(offset + PAGE)}
          >
            Next →
          </button>
          <span className="muted" style={{ marginLeft: "auto" }}>
            click a row for its evidence chain
          </span>
        </div>
      )}
    </div>
  );
}
