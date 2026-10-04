/**
 * Variant detail: key facts, plain-English explanation, the 5-step evidence
 * chain (numbered, icon-led), gene card and the genomic locus plot.
 */
import { Suspense, lazy, useEffect, useState } from "react";
import { Link } from "react-router-dom";
import { motion } from "motion/react";
import {
  Activity,
  Copy,
  Dna,
  GitBranch,
  Layers,
  Target,
  TriangleAlert,
} from "lucide-react";
import toast from "react-hot-toast";
import {
  api,
  type GeneCardResponse,
  type LocusResponse,
  type VariantDetail,
} from "../api/client";
import { CountUp, ScoreBar } from "./ui";

const LocusPlot = lazy(() => import("./LocusPlot"));

function Section({
  n,
  icon,
  role,
  body,
  low,
}: {
  n: number;
  icon: React.ReactNode;
  role: string;
  body: React.ReactNode;
  low?: boolean;
}) {
  return (
    <motion.div
      className={`chain-step ${low ? "low" : ""}`}
      initial={{ opacity: 0, x: -8 }}
      animate={{ opacity: 1, x: 0 }}
      transition={{ duration: 0.3, delay: 0.05 * n, ease: [0.22, 1, 0.36, 1] }}
    >
      <span className="num" aria-hidden>
        {icon}
      </span>
      <div className="role">
        {n} · {role}
      </div>
      <div className="body">{body}</div>
    </motion.div>
  );
}

export default function VariantDetailCard({
  runId,
  variantId,
}: {
  runId: string;
  variantId: string;
}) {
  const [detail, setDetail] = useState<VariantDetail | null>(null);
  const [locus, setLocus] = useState<LocusResponse | null>(null);
  const [gene, setGene] = useState<GeneCardResponse | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let dead = false;
    setError(null);
    setDetail(null);
    api
      .variant(runId, variantId)
      .then((d) => !dead && setDetail(d))
      .catch((e) => !dead && setError(String(e.message ?? e)));
    api
      .locus(runId, variantId)
      .then((l) => !dead && setLocus(l))
      .catch(() => !dead && setLocus(null));
    return () => {
      dead = true;
    };
  }, [runId, variantId]);

  const rawGenes = detail?.variant?.target_genes;
  const geneList = Array.isArray(rawGenes)
    ? rawGenes
    : typeof rawGenes === "string" && rawGenes
      ? rawGenes.split(",")
      : [];
  const primaryGene =
    geneList[0] ?? (detail?.variant?.nearest_gene as string | null) ?? "";

  useEffect(() => {
    if (!primaryGene) return;
    let dead = false;
    api
      .gene(primaryGene)
      .then((g) => !dead && setGene(g))
      .catch(() => !dead && setGene(null));
    return () => {
      dead = true;
    };
  }, [primaryGene]);

  if (error)
    return (
      <div className="card">
        <p className="error">{error}</p>
        <Link to={`/runs/${runId}`}>← back to ranked list</Link>
      </div>
    );

  if (!detail)
    return (
      <div className="card">
        <div className="skel" style={{ height: 22, width: "45%", marginBottom: 14 }} />
        <div className="grid2">
          <div className="skel" style={{ height: 220 }} />
          <div className="skel plot-skel" />
        </div>
        <p className="muted" style={{ marginTop: 12 }}>
          Loading {variantId}…
        </p>
      </div>
    );

  const v = detail.variant as Record<string, unknown>;
  const genes = geneList;
  const tracks = (
    Array.isArray(v.active_tracks) ? (v.active_tracks as string[]) : String(v.active_tracks ?? "").split(",")
  ).filter(Boolean);

  return (
    <>
      <p style={{ marginTop: 0 }}>
        <Link to={`/runs/${runId}`}>← Back to ranked list</Link>
      </p>

      <div className="ctx-bar">
        <div className="ctx-title">
          <span className="mono" style={{ fontSize: 16 }}>
            {variantId}
          </span>
          <span className="badge accent">
            rank <b>#{String(v.rank)}</b>
          </span>
          <span className="badge">
            <ScoreBar value={Number(v.score)} />
          </span>
          <button
            className="ghost small"
            onClick={() => {
              void navigator.clipboard?.writeText(variantId);
              toast.success("Variant ID copied", { duration: 1400 });
            }}
          >
            <Copy size={13} /> copy
          </button>
        </div>
      </div>

      <div className="grid2">
        <div>
          <div className="card">
            <h2>Variant facts</h2>
            <dl className="facts">
              <dt>position</dt>
              <dd className="mono">
                {String(v.chrom)}:{String(v.pos)}
              </dd>
              <dt>alleles</dt>
              <dd className="mono">
                {String(v.ref)} → {String(v.alt)}
              </dd>
              <dt>rank / score</dt>
              <dd>
                #{String(v.rank)} · {Number(v.score).toFixed(4)}
              </dd>
              <dt>consequence</dt>
              <dd>{String(v.consequence_class ?? "noncoding")}</dd>
              <dt>allele freq</dt>
              <dd>
                {v.af === null || v.af === undefined
                  ? "absent from common track (rare)"
                  : Number(v.af).toFixed(4)}
              </dd>
              <dt>DP / GQ</dt>
              <dd>
                {String(v.dp ?? "—")} / {String(v.gq ?? "—")}
              </dd>
              <dt>nearest gene</dt>
              <dd>{String(v.nearest_gene ?? "—")}</dd>
            </dl>

            <h3>Why this variant ranked</h3>
            <div className="callout">
              {String(v.explanation ?? detail.explanation ?? "")}
            </div>

            {detail.warnings?.length ? (
              <div className="callout amber" style={{ marginTop: 12 }}>
                <TriangleAlert
                  size={14}
                  style={{ verticalAlign: -2, marginRight: 6, color: "var(--warn)" }}
                />
                {detail.warnings.join(" · ")}
              </div>
            ) : null}
          </div>
        </div>

        <div className="card" data-tour="chain">
          <h2>Evidence chain</h2>
          <div className="chain">
            <Section
              n={1}
              icon={<Dna />}
              role="variant"
              body={
                <span className="mono">
                  {String(v.chrom)}:{String(v.pos)} {String(v.ref)}→{String(v.alt)}
                </span>
              }
            />
            <Section
              n={2}
              icon={<Layers />}
              role="active element (MM.1S)"
              body={
                <>
                  {tracks.map((t) => (
                    <span className="pill" key={t}>
                      {t}
                    </span>
                  ))}
                  {v.in_super_enhancer ? (
                    <span className="pill gold">★ super-enhancer</span>
                  ) : null}
                  <div className="muted" style={{ marginTop: 6 }}>
                    context score {Number(v.regulatory_context_score ?? 0).toFixed(3)}
                    {v.in_encode_ccre ? " · inside ENCODE cCRE" : ""}
                  </div>
                </>
              }
            />
            <Section
              n={3}
              icon={<Target />}
              role="motif disruption"
              low={v.motif_effect === "none" || !v.motif_effect}
              body={
                v.top_motif_tf ? (
                  <>
                    <span className={`pill ${String(v.motif_effect)}`}>
                      {String(v.motif_effect)}
                    </span>{" "}
                    <b>{String(v.top_motif_tf)}</b> site:{" "}
                    <span className="mono">
                      {Number(v.motif_ref_rel).toFixed(2)} →{" "}
                      {Number(v.motif_alt_rel).toFixed(2)} (Δ{" "}
                      {Number(v.motif_delta).toFixed(2)})
                    </span>
                  </>
                ) : (
                  <span className="muted">
                    no strong motif disruption at this variant — it ranked on context,
                    linkage and gene evidence
                  </span>
                )
              }
            />
            <Section
              n={4}
              icon={<GitBranch />}
              role="target gene"
              low={v.link_confidence === "low"}
              body={
                <>
                  {(genes.length ? genes : [v.nearest_gene]).filter(Boolean).map((g) => (
                    <span className="pill" key={String(g)}>
                      {String(g)}
                    </span>
                  ))}
                  <span className={`pill ${String(v.link_confidence ?? "low")}`}>
                    {String(v.link_method ?? "none")} · {String(v.link_confidence ?? "low")}{" "}
                    confidence
                  </span>
                </>
              }
            />
            <Section
              n={5}
              icon={<Activity />}
              role="myeloma relevance"
              body={
                <>
                  <div>
                    score{" "}
                    <b>
                      <CountUp
                        value={Number(v.myeloma_relevance_score ?? 0)}
                        format={(n) => n.toFixed(3)}
                      />
                    </b>{" "}
                    {v.is_known_mm_gene ? <span className="pill created">known MM gene</span> : null}{" "}
                    {v.depmap_mm_selective ? (
                      <span className="pill created">DepMap selective</span>
                    ) : null}{" "}
                    {v.expressed_in_mm ? <span className="pill">expressed in MM</span> : null}{" "}
                    {v.in_mm_gwas_locus || v.gwas_flag ? (
                      <span className="pill gold">in MM GWAS locus</span>
                    ) : null}
                  </div>
                  {v.depmap_mm_mean_effect !== null &&
                  v.depmap_mm_mean_effect !== undefined ? (
                    <div className="muted">
                      mean CRISPR gene effect in myeloma lines:{" "}
                      {Number(v.depmap_mm_mean_effect).toFixed(3)}
                    </div>
                  ) : null}
                </>
              }
            />
          </div>

          {primaryGene && (
            <div className="card" style={{ marginTop: 14, marginBottom: 0 }}>
              <h3 style={{ marginTop: 0 }}>
                Gene card · {primaryGene}
              </h3>
              {gene ? (
                <dl className="facts">
                  <dt>DepMap selective (myeloma)</dt>
                  <dd>{gene.depmap_mm_selective ? "yes" : "no"}</dd>
                  <dt>known myeloma gene</dt>
                  <dd>{gene.is_known_mm_gene ? `yes (${gene.known_gene_source})` : "no"}</dd>
                  <dt>expressed in MM.1S lines</dt>
                  <dd>{gene.expressed_in_mm ? "yes" : "no"}</dd>
                  <dt>variants linked here</dt>
                  <dd>{gene.linked_variant_count}</dd>
                </dl>
              ) : (
                <p className="muted">no gene card data</p>
              )}
            </div>
          )}
        </div>
      </div>

      <div className="card" data-tour="locus">
        <h2>
          Genomic context{" "}
          <span className="muted" style={{ fontWeight: 400 }}>
            ({locus ? `${locus.window[0].toLocaleString()}–${locus.window[1].toLocaleString()}` : "…"})
          </span>
        </h2>
        <Suspense fallback={<div className="skel plot-skel" />}>
          {locus ? (
            <LocusPlot locus={locus} />
          ) : (
            <p className="muted">loading locus…</p>
          )}
        </Suspense>
      </div>
    </>
  );
}
