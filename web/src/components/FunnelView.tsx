/**
 * "How we got here" — the filtering funnel in plain language.
 * Real row counts per stage; bars animate in on mount (reduced-motion safe).
 */
import { motion } from "motion/react";
import type { FunnelStepDto } from "../api/client";
import { CountUp } from "./ui";

const FRIENDLY: Record<string, { h: string; p: string }> = {
  s0: {
    h: "Load input variants",
    p: "Read the VCF or locus table, normalize alleles, keep autosomes.",
  },
  s1: {
    h: "Remove noise",
    p: "Drop common population variants, low-quality calls, blacklisted and coding regions.",
  },
  s2: {
    h: "Keep myeloma-active DNA",
    p: "A variant must overlap at least one of the 18 MM.1S regulatory tracks.",
  },
  s3: {
    h: "Scan TF binding motifs",
    p: "Score ±30 bp of real reference DNA on both strands against 10 JASPAR motifs.",
  },
  s4: {
    h: "Link elements to target genes",
    p: "ABC contact model first, then distance within 100 kb, then nearest-gene fallback.",
  },
  s5: {
    h: "Score myeloma relevance",
    p: "DepMap dependency, known myeloma genes, expression and GWAS-locus support.",
  },
  s6: {
    h: "Rank & explain",
    p: "Deterministic weighted sum — every variant gets a reviewable evidence chain.",
  },
};

export default function FunnelView({
  steps,
  onPickVariant,
}: {
  steps: FunnelStepDto[];
  onPickVariant?: () => void;
}) {
  if (!steps.length)
    return (
      <div className="card">
        <p className="muted">No funnel data for this run.</p>
      </div>
    );

  const max = Math.max(1, ...steps.map((s) => s.rows_in));
  const first = steps[0];
  const last = steps[steps.length - 1];

  return (
    <div className="card">
      <span className="section-kicker">How we got here</span>
      <div className="funnel-head" style={{ marginTop: 8 }}>
        <span className="big">
          <CountUp value={first.rows_in} duration={0.7} />
        </span>
        <span className="arrow">variants in</span>
        <span className="arrow">→</span>
        <span className="big" style={{ color: "var(--accent)" }}>
          <CountUp value={last.rows_out} duration={1.1} />
        </span>
        <span className="arrow">ranked out · {steps.length} stages · every removal accounted for</span>
      </div>

      {steps.map((s, i) => {
        const pct = Math.max(0.4, (s.rows_out / max) * 100);
        const reasons = Object.entries(s.removed_reasons ?? {});
        const friendly = FRIENDLY[s.step];
        return (
          <div key={s.step}>
            <div className="funnel-row">
              <div className="funnel-label">
                <span className="s-code">{s.step}</span>
                <span>{friendly?.h ?? s.label}</span>
              </div>
              <div className="funnel-bar-wrap">
                <motion.div
                  className="funnel-bar"
                  initial={{ width: 0 }}
                  animate={{ width: `${pct}%` }}
                  transition={{
                    duration: 0.7,
                    delay: 0.08 * i,
                    ease: [0.22, 1, 0.36, 1],
                  }}
                >
                  {s.rows_out.toLocaleString()}
                </motion.div>
              </div>
              <div className={`funnel-removed ${s.removed ? "" : "none"}`}>
                {s.removed > 0 ? `−${s.removed.toLocaleString()} removed` : "kept as-is"}
              </div>
            </div>
            {friendly && <div className="funnel-note">{friendly.p}</div>}
            {reasons.length > 0 && (
              <div className="funnel-note">
                {reasons.map(([k, v]) => (
                  <span className="reason" key={k}>
                    {k.replace(/_/g, " ")}: {v.toLocaleString()}
                  </span>
                ))}
              </div>
            )}
            {s.note && <div className="funnel-note muted">{s.note}</div>}
          </div>
        );
      })}

      {onPickVariant && (
        <p style={{ marginTop: 14, marginBottom: 0 }}>
          <button className="secondary" onClick={onPickVariant}>
            Inspect the ranked list →
          </button>
        </p>
      )}
    </div>
  );
}
