/**
 * Data manifest & provenance: every downloaded resource with its verification
 * and content hash. Status chips filter the table.
 */
import { useEffect, useMemo, useState } from "react";
import { Copy } from "lucide-react";
import toast from "react-hot-toast";
import { api, type Manifest } from "../api/client";
import { Help, Skeleton } from "./ui";

export default function ManifestView() {
  const [m, setM] = useState<Manifest | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [status, setStatus] = useState("");

  useEffect(() => {
    api
      .manifest()
      .then((x) => setM(x as unknown as Manifest))
      .catch((e) => setError(String(e.message ?? e)));
  }, []);

  const counts = useMemo(() => {
    const c: Record<string, number> = {};
    for (const r of m?.resources ?? []) c[r.status] = (c[r.status] ?? 0) + 1;
    return c;
  }, [m]);

  const rows = (m?.resources ?? []).filter((r) => !status || r.status === status);

  if (error)
    return (
      <div className="card">
        <h2>Data manifest</h2>
        <p className="error">{error}</p>
      </div>
    );

  if (!m)
    return (
      <div className="card">
        <div className="skel" style={{ height: 20, width: "40%", marginBottom: 14 }} />
        <Skeleton w="100%" />
        <p className="muted" style={{ marginTop: 12 }}>
          Loading manifest…
        </p>
      </div>
    );

  return (
    <>
      <div className="ctx-bar">
        <div className="ctx-title">
          Data manifest & provenance{" "}
          <span className="badge mono">{m.resources.length} resources</span>
          <span className="badge mono">{m.genome_build}</span>
          <span className="badge mono">
            {m.manifest_hash}
            <button
              className="ghost small"
              style={{ marginLeft: 6, padding: "0 4px" }}
              onClick={() => {
                void navigator.clipboard?.writeText(m.manifest_hash ?? "");
                toast.success("Manifest hash copied");
              }}
              aria-label="copy manifest hash"
            >
              <Copy size={11} />
            </button>
          </span>
        </div>
        <div className="spacer" />
        <span className="muted">generated {String(m.generated_at).slice(0, 19)}</span>
      </div>

      <div className="card">
        <div className="status-chips">
          <button
            className={`status-chip ${status === "" ? "on" : ""}`}
            onClick={() => setStatus("")}
          >
            All <span className="n">{m.resources.length}</span>
          </button>
          {Object.entries(counts).map(([s, n]) => (
            <button
              key={s}
              className={`status-chip ${status === s ? "on" : ""}`}
              onClick={() => setStatus(s)}
            >
              {s} <span className="n">{n}</span>
            </button>
          ))}
        </div>

        <table>
          <thead>
            <tr>
              <th>resource</th>
              <th>category</th>
              <th>accession / source</th>
              <th>size</th>
              <th>status</th>
              <th>verification</th>
            </tr>
          </thead>
          <tbody>
            {rows.map((r) => (
              <tr key={r.id + r.path}>
                <td className="mono">{r.id}</td>
                <td>{r.category}</td>
                <td>
                  <a href={r.source_url} target="_blank" rel="noreferrer">
                    {r.accession || r.source_url}
                  </a>
                </td>
                <td className="mono">{(r.size_bytes / 1e6).toFixed(1)} MB</td>
                <td>
                  <span className={`badge ${r.status === "ok" ? "ok" : "warn"}`}>
                    {r.status}
                  </span>
                </td>
                <td className="muted">{r.verification}</td>
              </tr>
            ))}
          </tbody>
        </table>
        <p className="muted" style={{ marginTop: 10 }}>
          Every result embeds this manifest hash, so any ranked list can be traced
          back to the exact inputs and configuration that produced it.
          <Help text="The hash is a content digest over the manifest entries (URLs, sizes, checks). A changed dataset produces a different hash." />
        </p>
      </div>

      <div className="card">
        <span className="section-kicker">Honest notes</span>
        <h2 style={{ marginTop: 6 }}>Limitations & substitutions</h2>
        <ul style={{ lineHeight: 1.7, fontSize: 13.5 }}>
          <li>
            <strong>gnomAD AF:</strong> gnomAD sites VCFs are ~53 GB per chromosome —
            infeasible here. Bulk AF filtering uses the UCSC dbSNP-151 common track
            (real per-allele frequencies); shortlist variants can additionally be
            checked via the Ensembl VEP/gnomAD browser APIs.
          </li>
          <li>
            <strong>ABC biosample proxy:</strong> Nasser et al. ABC predictions
            contain no plasma-cell biosample; the closest available biosample
            (B-cell/Roadmap GM12878) is used and explicitly labelled as a proxy on
            every link.
          </li>
          <li>
            <strong>DepMap via figshare:</strong> the DepMap portal is
            Cloudflare-gated; releases are fetched from the public figshare records
            instead.
          </li>
          <li>
            <strong>GWAS Catalog v2:</strong> the legacy v1 REST API returns HTTP 429
            (deprecated); the v2 API is used.
          </li>
          <li>
            <strong>TF coverage:</strong> ENCODE MM.1S TF ChIP-seq exists only for
            EZH2/CTCF/EZH2phospho — IRF4/MAF/MYC occupancy is not claimed; motif
            disruption uses JASPAR PWMs instead.
          </li>
          <li>
            <strong>GEO GSE160335:</strong> only a 2.8 GB RAW.tar was reachable and
            individual files 404'd — recorded as failed (non-critical) rather than
            silently skipped.
          </li>
          <li>
            <strong>Mode A input:</strong> GIAB HG001/NA12878 is a healthy germline
            genome — no somatic myeloma drivers are expected; the run demonstrates
            funnel specificity, not cancer discovery.
          </li>
        </ul>
      </div>
    </>
  );
}
