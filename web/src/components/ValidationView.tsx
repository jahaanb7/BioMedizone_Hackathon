/**
 * Validation: enrichment statistics with explicit independence flags,
 * weight-perturbation robustness, the Mode-A funnel and publication figures.
 * Every number comes from /api/validation (scripts/run_validation.py).
 */
import { useEffect, useState } from "react";
import { motion } from "motion/react";
import { api } from "../api/client";
import FunnelView from "./FunnelView";
import { Help, Skeleton } from "./ui";

type Enrich = {
  test?: string;
  odds_ratio?: number | null;
  p_value?: number | null;
  independent_of_scoring?: boolean;
  note?: string;
  contingency?: Record<string, unknown> | null;
  recovered?: number;
  n_loci?: number;
  recovery_rate?: number;
  expected_rate_random_ranking?: number;
};

type Payload = {
  generated_at?: string;
  mode_a?: {
    run_id?: string;
    total_ranked?: number;
    elapsed_seconds?: number;
    manifest_hash?: string;
    funnel?: unknown[];
    sensitivity?: {
      n_perturbations?: number;
      top_n?: number;
      jaccard_mean?: number;
      jaccard_min?: number;
      spearman_mean?: number;
      weight_scheme?: string;
    };
    warnings?: string[];
  };
  enrichment_tests?: Enrich[];
  mode_b?: Enrich & { per_locus?: Record<string, unknown>[] };
  figures?: string[];
};

const FRIENDLY_FIG: Record<string, string> = {
  "fig_funnel.png": "Filtering funnel — where the 3.9M input variants went",
  "fig_score_hist.png": "Score distribution across ranked variants",
  "fig_enrichment.png": "Top-ranked variants vs. the rest (enrichment)",
  "fig_mode_b_lead_ranks.png": "Mode B — ranks of published lead SNPs",
};

const TEST_NAME: Record<string, string> = {
  top1000_vs_rest__inside_ENCODE_SCREEN_cCRE:
    "Top 1,000 variants sit in independent cCREs",
  top500_linked_genes_vs_rest__known_myeloma_genes:
    "Top-500 linked genes are known myeloma genes",
  mode_b_published_lead_in_top3: "Mode B: published lead SNPs recovered in top 3",
};

function fmtP(p?: number | null): string {
  if (p === null || p === undefined) return "n/a";
  if (p === 0) return "< 1e-300";
  return p.toExponential(2);
}

function statName(t?: string): string {
  if (!t) return "test";
  return TEST_NAME[t] ?? t.replace(/__/g, " · ").replace(/_/g, " ");
}

function TestCard({ t }: { t: Enrich }) {
  const isLead = typeof t.recovered === "number";
  const big = isLead
    ? `${t.recovered}/${t.n_loci ?? 0}`
    : t.odds_ratio != null
      ? t.odds_ratio.toFixed(2)
      : "—";
  const unit = isLead ? "loci in top 3" : "odds ratio";
  return (
    <motion.div
      className="stat-card"
      initial={{ opacity: 0, y: 10 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.35, ease: [0.22, 1, 0.36, 1] }}
    >
      <div className="test-name">{statName(t.test)}</div>
      <div className="or">{big}</div>
      <div className="p">
        {unit} · p = {fmtP(t.p_value)}
        {isLead && t.recovery_rate != null && (
          <div style={{ marginTop: 4 }}>
            recovery {(t.recovery_rate * 100).toFixed(1)}% vs{" "}
            {((t.expected_rate_random_ranking ?? 0) * 100).toFixed(3)}% expected at random
          </div>
        )}
      </div>
      <div className="flags">
        <span className={`badge ${t.independent_of_scoring ? "ok" : "warn"}`}>
          {t.independent_of_scoring ? "independent of scoring" : "partially circular"}
        </span>
        {t.contingency && (
          <span className="badge mono">
            {Object.entries(t.contingency)
              .map(([k, v]) => `${k}=${typeof v === "number" ? v.toLocaleString() : v}`)
              .join(" · ")}
          </span>
        )}
      </div>
      {t.note && <div className="note">{t.note}</div>}
    </motion.div>
  );
}

function renderValue(v: unknown): React.ReactNode {
  if (v === null || v === undefined) return <span className="muted">n/a</span>;
  if (typeof v === "number") return Number.isInteger(v) ? v.toLocaleString() : v.toFixed(4);
  if (typeof v === "boolean") return v ? "yes" : "no";
  if (typeof v === "string") return v;
  if (Array.isArray(v))
    return (
      <ul style={{ margin: 0, paddingLeft: 18 }}>
        {v.map((x, i) => (
          <li key={i}>{renderValue(x)}</li>
        ))}
      </ul>
    );
  if (typeof v === "object") {
    return (
      <table>
        <tbody>
          {Object.entries(v as Record<string, unknown>).map(([k, val]) => (
            <tr key={k}>
              <td className="muted">{k}</td>
              <td>{renderValue(val)}</td>
            </tr>
          ))}
        </tbody>
      </table>
    );
  }
  return String(v);
}

export default function ValidationView() {
  const [data, setData] = useState<Payload | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    api
      .validation()
      .then((d) => setData(d as Payload))
      .catch((e) => setError(String(e.message ?? e)));
  }, []);

  if (error)
    return (
      <div className="card">
        <h2>Validation</h2>
        <p className="error">{error}</p>
        <p className="muted">Run scripts/run_validation.py to compute it.</p>
      </div>
    );

  if (!data)
    return (
      <div className="card">
        <div className="skel" style={{ height: 20, width: "35%", marginBottom: 16 }} />
        <div className="stat-cards">
          {[0, 1, 2].map((i) => (
            <Skeleton key={i} w="100%" />
          ))}
        </div>
        <p className="muted" style={{ marginTop: 14 }}>
          Loading validation results…
        </p>
      </div>
    );

  const a = data.mode_a ?? {};
  const sens = a.sensitivity ?? {};
  const tests = data.enrichment_tests ?? [];
  const figures = data.figures ?? [];
  const mb = data.mode_b;

  return (
    <>
      <div className="ctx-bar">
        <div className="ctx-title">
          Validation{" "}
          <span className="badge accent">statistical, on real data</span>
          {a.total_ranked && (
            <span className="badge mono">{a.total_ranked.toLocaleString()} ranked</span>
          )}
          {a.manifest_hash && (
            <span className="badge mono">{a.manifest_hash.slice(0, 16)}</span>
          )}
        </div>
        <div className="spacer" />
        {data.generated_at && (
          <span className="muted">generated {String(data.generated_at).slice(0, 19)}</span>
        )}
      </div>

      <div className="card" data-tour="validation-stats">
        <span className="section-kicker">Enrichment tests</span>
        <h2 style={{ marginTop: 6 }}>
          Do the top-ranked variants look like real regulatory signal?
          <Help text="Each test compares the top of the ranking against the rest of the list. Independence flags tell you whether the test re-uses any input of the score." />
        </h2>
        <div className="stat-cards" style={{ marginTop: 12 }}>
          {tests.map((t, i) => (
            <TestCard key={i} t={t} />
          ))}
        </div>
      </div>

      <div className="grid2">
        <div className="card">
          <span className="section-kicker">Robustness</span>
          <h2 style={{ marginTop: 6 }}>Weight perturbation sweep</h2>
          <div className="or" style={{ fontSize: 40, fontWeight: 800, margin: "8px 0 0" }}>
            {sens.spearman_mean != null ? sens.spearman_mean.toFixed(3) : "—"}
          </div>
          <p className="muted" style={{ marginTop: 0 }}>
            mean Spearman correlation of rankings after ±50% weight perturbation
          </p>
          <dl className="facts">
            <dt>top-50 Jaccard (mean / min)</dt>
            <dd>
              {sens.jaccard_mean?.toFixed(3) ?? "—"} / {sens.jaccard_min?.toFixed(3) ?? "—"}
            </dd>
            <dt>perturbations</dt>
            <dd>{sens.n_perturbations ?? "—"}</dd>
          </dl>
          {sens.weight_scheme && <p className="muted">{sens.weight_scheme}</p>}
          {a.warnings && a.warnings.length > 0 && (
            <details className="notes">
              <summary>{a.warnings.length} data notes</summary>
              <ul>
                {a.warnings.map((w, i) => (
                  <li key={i}>{w}</li>
                ))}
              </ul>
            </details>
          )}
        </div>

        <div>
          {Array.isArray(a.funnel) && a.funnel.length > 0 ? (
            <FunnelView steps={a.funnel as never} />
          ) : (
            <div className="card">
              <p className="muted">no funnel in payload</p>
            </div>
          )}
        </div>
      </div>

      {mb && mb.per_locus && mb.per_locus.length > 0 && (
        <div className="card">
          <span className="section-kicker">Mode B</span>
          <h2 style={{ marginTop: 6 }}>Per-locus lead-SNP ranks</h2>
          <table>
            <thead>
              <tr>
                <th>locus</th>
                <th>candidates ranked</th>
                <th>lead SNP rank</th>
                <th>in top 3?</th>
              </tr>
            </thead>
            <tbody>
              {mb.per_locus.map((p, i) => (
                <tr key={i}>
                  <td className="mono">{String(p.locus_id)}</td>
                  <td className="mono">{Number(p.n_candidates).toLocaleString()}</td>
                  <td className="mono">#{Number(p.lead_rank).toLocaleString()}</td>
                  <td>
                    {p.recovered_at_k ? (
                      <span className="badge ok">yes</span>
                    ) : (
                      <span className="badge">no</span>
                    )}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
          {mb.note && <p className="muted">{mb.note}</p>}
        </div>
      )}

      {figures.length > 0 && (
        <div className="card">
          <span className="section-kicker">Figures</span>
          <div className="figure-grid" style={{ marginTop: 12 }}>
            {figures.map((f) => (
              <figure key={f}>
                <img src={api.figureUrl(f)} alt={FRIENDLY_FIG[f] ?? f} loading="lazy" />
                <figcaption>{FRIENDLY_FIG[f] ?? f}</figcaption>
              </figure>
            ))}
          </div>
        </div>
      )}

      <div className="card">
        <span className="section-kicker">Raw payload</span>
        <h2 style={{ marginTop: 6 }}>Everything, unfiltered</h2>
        {Object.entries(data).map(([k, v]) =>
          ["figures", "enrichment_tests", "mode_a", "mode_b"].includes(k) ? null : (
            <div key={k} style={{ marginBottom: 12 }}>
              <h3>{k.replace(/_/g, " ")}</h3>
              {renderValue(v)}
            </div>
          ),
        )}
      </div>
    </>
  );
}
