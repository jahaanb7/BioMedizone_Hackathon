import { useEffect, useRef } from "react";
import Plotly from "plotly.js-dist-min";
import type { LocusResponse } from "../api/client";

const TRACK_COLORS: Record<string, string> = {
  accessibility: "#4da3ff",
  h3k27ac: "#ff6b6b",
  h3k4me1: "#ffb454",
  h3k4me3: "#7ee0a3",
  super_enhancer: "#c792ea",
  tf_peak: "#89ddff",
};

function categoryColor(category: string, name: string): string {
  const key = category.toLowerCase();
  for (const k of Object.keys(TRACK_COLORS))
    if (key.includes(k) || name.toLowerCase().includes(k)) return TRACK_COLORS[k];
  return "#5c6b8a";
}

/** Browser-like track plot: peaks as rectangles, genes as arrows, variant marked. */
export default function LocusPlot({ locus }: { locus: LocusResponse }) {
  const ref = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const el = ref.current;
    if (!el) return;
    const [start, end] = locus.window;
    const data: unknown[] = [];
    let row = 0;
    const rowLabels: string[] = [];
    const tracks = locus.tracks ?? [];
    type GeneRow = { gene: string; start: number; end: number; strand: string };
    const genes = (locus.genes ?? []) as GeneRow[];

    for (const t of tracks) {
      const color = categoryColor(t.category ?? "", t.name);
      const ivs = (t.intervals ?? []) as number[][];
      data.push({
        type: "bar",
        orientation: "h",
        x: ivs.map((iv) => ({ x0: iv[0], x1: iv[1] })),
        y: ivs.map(() => row),
        width: ivs.map((iv) => iv[1] - iv[0]),
        base: ivs.map((iv) => iv[0]),
        marker: { color },
        name: t.display_name || t.name,
        showlegend: false,
        hovertemplate: `${t.display_name || t.name}<br>%{base|%d}–%{x|%d}<extra></extra>`,
      } as unknown);
      rowLabels.push(t.display_name || t.name);
      row += 1;
    }

    for (const g of genes) {
      data.push({
        type: "bar",
        orientation: "h",
        x: [{ x0: g.start, x1: g.end }],
        y: [row],
        width: [g.end - g.start],
        base: [g.start],
        marker: { color: g.strand === "-" ? "#8b97b3" : "#dce3f2", opacity: 0.85 },
        showlegend: false,
        hovertemplate: `${g.gene} (${g.strand})<extra></extra>`,
      } as unknown);
      rowLabels.push(g.gene);
      row += 1;
    }

    data.push({
      type: "scatter",
      mode: "lines",
      x: [locus.pos, locus.pos],
      y: [-1, Math.max(1, row - 1)],
      line: { color: "#ffd166", width: 2, dash: "dot" },
      name: locus.variant_id,
      showlegend: true,
      hoverinfo: "skip",
    });

    const layout = {
      title: { text: `${locus.variant_id} @ ${locus.chrom}`, font: { size: 13 } },
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      font: { color: "#8b97b3", size: 11 },
      margin: { l: 130, r: 20, t: 34, b: 30 },
      xaxis: { range: [start, end], title: "genomic position (0-based)" },
      yaxis: {
        range: [-1, row],
        tickmode: "array",
        tickvals: rowLabels.map((_, i) => i),
        ticktext: rowLabels,
        automargin: true,
      },
      bargap: 0.2,
      showlegend: false,
      height: Math.max(240, 60 + row * 22),
      hovermode: "closest",
    };

    Plotly.react(el, data, layout, { displayModeBar: false, responsive: true });
    return () => {
      Plotly.purge(el);
    };
  }, [locus]);

  return (
    <>
      <div className="plot" ref={ref} />
      <p className="muted">
      {locus.ccre_count ?? 0} ENCODE cCREs in window ·{" "}
      {(locus.tracks ?? []).length} active tracks ·{" "}
        {locus.genes?.length ?? 0} genes
        {locus.cached ? " · precomputed demo payload (offline)" : ""}
      </p>
    </>
  );
}
