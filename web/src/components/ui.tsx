/**
 * Shared UI primitives — one visual language for evidence chips, score bars,
 * animated counts, skeletons and empty states.
 */
import { useEffect, useState, type ReactNode } from "react";
import { animate, useMotionValue } from "motion/react";
import { clsx, type ClassValue } from "clsx";
import { twMerge } from "tailwind-merge";
import { SearchX } from "lucide-react";

export function cx(...v: ClassValue[]) {
  return twMerge(clsx(v));
}

/** Numeric count-up driven by motion's animate(); respects reduced motion. */
export function CountUp({
  value,
  duration = 0.9,
  format = (n: number) => Math.round(n).toLocaleString(),
}: {
  value: number;
  duration?: number;
  format?: (n: number) => string;
}) {
  const mv = useMotionValue(0);
  const [shown, setShown] = useState(0);
  useEffect(() => {
    const controls = animate(mv, value, {
      duration,
      ease: [0.22, 1, 0.36, 1],
      onUpdate: (v) => setShown(v),
    });
    return () => controls.stop();
  }, [value, duration, mv]);
  return <>{format(shown)}</>;
}

/** Score bar: value is 0..1 (the score is normalized). */
export function ScoreBar({ value }: { value: number }) {
  const pct = Math.max(0, Math.min(1, value)) * 100;
  return (
    <span className="score-cell" title={`score ${value.toFixed(4)}`}>
      <span className="num">{value.toFixed(4)}</span>
      <span className="score-bar">
        <i style={{ width: `${pct}%` }} />
      </span>
    </span>
  );
}

/** Context-score mini gauge (0..1). */
export function Gauge({ value }: { value: number }) {
  const pct = Math.max(0, Math.min(1, value || 0)) * 100;
  return (
    <span className="ctx-cell" title="regulatory context score (0–1)">
      <span className="gauge">
        <i style={{ width: `${pct}%` }} />
      </span>
      {(value || 0).toFixed(2)}
    </span>
  );
}

export function MotifChip({ effect }: { effect?: string | null }) {
  if (!effect || effect === "none")
    return <span className="muted">—</span>;
  return <span className={`pill ${effect}`}>{effect}</span>;
}

/** Three evidence dots: DepMap dependency, known myeloma gene, GWAS locus. */
export function EvidenceDots({
  depmap,
  known,
  gwas,
}: {
  depmap?: boolean;
  known?: boolean;
  gwas?: boolean;
}) {
  return (
    <span className="dots">
      <i className={cx("dot", depmap && "on-depmap")} title={depmap ? "DepMap: myeloma-selective dependency" : "no DepMap dependency"} />
      <i className={cx("dot", known && "on-known")} title={known ? "known myeloma gene" : "not in known myeloma gene set"} />
      <i className={cx("dot", gwas && "on-gwas")} title={gwas ? "in a myeloma GWAS locus" : "outside GWAS loci"} />
    </span>
  );
}

export function Skeleton({ w = "100%" }: { w?: string }) {
  return <div className="skel" style={{ width: w }} />;
}

export function SkeletonRows({ rows = 8, cols = 7 }: { rows?: number; cols?: number }) {
  return (
    <>
      {Array.from({ length: rows }).map((_, r) => (
        <tr key={r}>
          {Array.from({ length: cols }).map((__, c) => (
            <td key={c}>
              <Skeleton w={c === 1 ? "85%" : c === 0 ? "30%" : "60%"} />
            </td>
          ))}
        </tr>
      ))}
    </>
  );
}

export function EmptyState({
  title,
  body,
  action,
  actionLabel,
}: {
  title: string;
  body?: ReactNode;
  action?: () => void;
  actionLabel?: string;
}) {
  return (
    <div className="empty-state">
      <div className="ico">
        <SearchX size={30} strokeWidth={1.6} />
      </div>
      <h4>{title}</h4>
      {body && <p>{body}</p>}
      {action && (
        <button className="secondary" onClick={action}>
          {actionLabel ?? "Reset filters"}
        </button>
      )}
    </div>
  );
}

/** Question-mark tooltip chip (pure CSS hover, no runtime dependency). */
export function Help({ text }: { text: string }) {
  return (
    <span className="tip" data-tip={text} role="img" aria-label="help">
      ?
    </span>
  );
}
