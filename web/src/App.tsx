import { useCallback, useEffect, useState, type ReactNode } from "react";
import {
  Link,
  Route,
  Routes,
  useLocation,
  useNavigate,
  useParams,
} from "react-router-dom";
import { AnimatePresence, MotionConfig, motion } from "motion/react";
import { Toaster, toast } from "react-hot-toast";
import {
  ArrowRight,
  BookOpen,
  CircleHelp,
  Dna,
  Download,
  GraduationCap,
  Layers,
  ListOrdered,
  PlayCircle,
  RefreshCw,
  ShieldCheck,
  Sparkles,
  FlaskConical,
} from "lucide-react";
import {
  api,
  type FunnelStepDto,
  type Health,
  type JobStatus,
  type ModeInfo,
  type RunCreateResponse,
  type ValidationResults,
} from "./api/client";
import FunnelView from "./components/FunnelView";
import VariantsTable from "./components/VariantsTable";
import VariantDetailCard from "./components/VariantDetailCard";
import ValidationView from "./components/ValidationView";
import ManifestView from "./components/ManifestView";
import Glossary from "./components/Glossary";
import Tour from "./components/Tour";
import { CountUp, Skeleton } from "./components/ui";

/* ------------------------------------------------------------------ utils */

function usePoll(runId: string | null | undefined, onTick: () => void) {
  useEffect(() => {
    if (!runId) return;
    const t = setInterval(onTick, 1500);
    return () => clearInterval(t);
  }, [runId, onTick]);
}

/** Route transition wrapper — short fade+slide, honours reduced motion. */
function Page({ children }: { children: ReactNode }) {
  return (
    <motion.div
      className="page-in"
      initial={{ opacity: 0, y: 10 }}
      animate={{ opacity: 1, y: 0 }}
      exit={{ opacity: 0, y: -6 }}
      transition={{ duration: 0.22, ease: [0.22, 1, 0.36, 1] }}
    >
      {children}
    </motion.div>
  );
}

/* ------------------------------------------------------------ analyze card */

function ModePanel({ modes, onStarted }: { modes: ModeInfo[]; onStarted: (r: RunCreateResponse) => void }) {
  const [file, setFile] = useState<File | null>(null);
  const [mode, setMode] = useState("A");
  const [busy, setBusy] = useState(false);
  const [err, setErr] = useState<string | null>(null);

  const start = async () => {
    setBusy(true);
    setErr(null);
    try {
      const r = await api.startRun(mode as "A" | "B" | "C", file ?? undefined);
      onStarted(r);
    } catch (e) {
      setErr(e instanceof Error ? e.message : String(e));
    } finally {
      setBusy(false);
    }
  };

  return (
    <div className="card">
      <h2>
        <FlaskConical size={15} style={{ marginRight: 7, verticalAlign: -2 }} />
        Analyze your own data
      </h2>
      <p className="muted" style={{ marginTop: -2 }}>
        Run the full funnel on a new VCF or on built-in loci. Results appear under
        Runs below.
      </p>
      <label>Mode</label>
      <select value={mode} onChange={(e) => setMode(e.target.value)} style={{ width: "100%" }}>
        {modes.map((m) => (
          <option key={m.id} value={m.id}>
            {m.id} — {m.name}
          </option>
        ))}
      </select>
      <p className="muted">{modes.find((m) => m.id === mode)?.description}</p>

      {mode === "A" && (
        <>
          <input
            type="file"
            accept=".vcf,.gz,.bgz"
            onChange={(e) => setFile(e.target.files?.[0] ?? null)}
          />
          <div style={{ marginTop: 10 }}>
            <button disabled={busy || !file} onClick={start}>
              Run funnel
            </button>
          </div>
          <p className="muted">
            VCF must be GRCh38/hg38. The GIAB demo runs without an upload.
          </p>
        </>
      )}
      {mode !== "A" && (
        <div style={{ marginTop: 10 }}>
          <button disabled={busy} onClick={start}>
            Run mode {mode}
          </button>
        </div>
      )}
      {err && <p className="error">{err}</p>}
    </div>
  );
}

function RunList({ refreshKey }: { refreshKey: number }) {
  const [runs, setRuns] = useState<JobStatus[]>([]);
  useEffect(() => {
    api.runs().then(setRuns).catch(() => setRuns([]));
  }, [refreshKey]);
  if (!runs.length) return null;
  return (
    <div className="card">
      <h2>Recent runs</h2>
      <table>
        <thead>
          <tr>
            <th>Run</th>
            <th>Mode</th>
            <th>Status</th>
            <th>Ranked</th>
            <th>Input</th>
          </tr>
        </thead>
        <tbody>
          {runs.map((r) => (
            <tr key={r.run_id}>
              <td className="mono">
                <Link to={`/runs/${r.run_id}`}>{r.run_id}</Link>
              </td>
              <td>{r.mode}</td>
              <td>
                <span
                  className={`badge ${r.status === "done" ? "ok" : r.status === "error" ? "warn" : ""}`}
                >
                  {r.status}
                </span>
              </td>
              <td className="mono">{r.total_ranked ? r.total_ranked.toLocaleString() : "—"}</td>
              <td className="muted">{r.input_description?.slice(0, 44)}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

/* ------------------------------------------------------------------ home */

const PIPELINE: { n: string; h: string; p: string }[] = [
  { n: "s0", h: "Load variants", p: "Read the VCF or GWAS loci; normalize alleles; autosomes only." },
  { n: "s1", h: "Remove noise", p: "Drop common population variants and low-quality calls." },
  { n: "s2", h: "Keep active DNA", p: "Retain only variants inside DNA that is active in myeloma cells." },
  { n: "s3", h: "Scan TF motifs", p: "Score ±30 bp of real reference sequence on both strands, 10 motifs." },
  { n: "s4", h: "Link to genes", p: "ABC contact model first, then distance, then nearest gene — labelled." },
  { n: "s5", h: "Ask myeloma", p: "DepMap dependency, known-gene and GWAS evidence per target gene." },
  { n: "s6", h: "Rank & explain", p: "Deterministic weighted score and a plain-English evidence chain." },
];

function HomePage({ modes, startTour }: { modes: ModeInfo[]; startTour: () => void }) {
  const nav = useNavigate();
  const [busy, setBusy] = useState<null | "A" | "B">(null);
  const [err, setErr] = useState<string | null>(null);
  const [facts, setFacts] = useState<Record<string, string>>({});

  const onStarted = (r: RunCreateResponse) => {
    if (r.status === "done" || r.status === "queued") nav(`/runs/${r.run_id}`);
  };

  /* real funnel numbers from the validation payload (no invented figures) */
  useEffect(() => {
    api
      .validation()
      .then((v) => {
        const val = v as ValidationResults & {
          mode_a?: { funnel?: FunnelStepDto[]; sensitivity?: Record<string, number> };
        };
        const f = val.mode_a?.funnel ?? [];
        const sens = val.mode_a?.sensitivity ?? {};
        const next: Record<string, string> = {};
        if (f[1]) next.s1 = `${f[1].rows_out.toLocaleString()} kept`;
        if (f[2]) next.s2 = `${f[2].rows_out.toLocaleString()} in active DNA`;
        const m = /broken=(\d+), created=(\d+)/.exec(f[3]?.note ?? "");
        if (m) next.s3 = `${m[1]} broken · ${m[2]} created`;
        const a = /ABC=(\d+)/.exec(f[4]?.note ?? "");
        if (a) next.s4 = `${Number(a[1]).toLocaleString()} ABC links`;
        if (typeof sens.spearman_mean === "number")
          next.s6 = `ρ ${sens.spearman_mean.toFixed(3)} under ±50% weight perturbation`;
        setFacts(next);
      })
      .catch(() => undefined);
  }, []);

  const openDemo = async (m: "A" | "B") => {
    setBusy(m);
    setErr(null);
    try {
      const r = await api.demo(m);
      onStarted({
        run_id: r.run_id,
        status: "done",
        mode: m,
        input_description: "precomputed demo",
        cached: true,
      });
      nav(`/runs/${r.run_id}`);
    } catch (e) {
      setErr(e instanceof Error ? e.message : String(e));
      toast.error("Could not load the precomputed demo.");
    } finally {
      setBusy(null);
    }
  };

  const card = (m: ModeInfo, variant: "a" | "b") => {
    const ready = m.demo_available;
    const total = m.demo_total_ranked;
    return (
      <div
        className={`entry-card ${variant}`}
        data-tour={variant === "a" ? "entry-card-A" : undefined}
      >
        <span className="mode-tag">Mode {m.id} · precomputed demo</span>
        <h3>{variant === "a" ? "Whole-genome scan" : "GWAS loci fine-mapping"}</h3>
        <div className="count">
          {total != null ? <CountUp value={total} /> : ready ? <Skeleton w="160px" /> : "—"}
          <small> ranked variants</small>
        </div>
        <p>{m.description}</p>
        <div className="meta">
          <span className="badge">{variant === "a" ? "3.9M input variants" : "37 GWAS loci + LD proxies"}</span>
          <span className="badge">offline</span>
        </div>
        <button
          className="open"
          disabled={!ready || busy !== null}
          onClick={() => openDemo(m.id === "A" ? "A" : "B")}
          aria-label={`Open ranked list, mode ${m.id}`}
        >
          {busy === m.id ? "Opening…" : "Open ranked list"}
          <ArrowRight size={15} />
        </button>
      </div>
    );
  };

  return (
    <>
      <section className="hero" data-tour="hero">
        <span className="section-kicker">Multiple myeloma · regulatory variant prioritization</span>
        <h1>
          Ranked variants, <span className="grad">with the evidence attached.</span>
        </h1>
        <p className="sub">
          Myelovar filters millions of non-coding variants down to the ones inside
          myeloma-active DNA, scores the transcription-factor motifs they disrupt,
          links each to a target gene and ranks everything — every result carries a
          reviewable chain, not just a score.
        </p>
        <div className="cta-row">
          <button onClick={() => openDemo("A")} disabled={busy !== null || !modes.find((m) => m.id === "A")?.demo_available}>
            {busy === "A" ? "Opening…" : "Open ranked list"}
            <ArrowRight size={15} />
          </button>
          <button className="secondary" onClick={startTour}>
            <GraduationCap size={15} />
            Take the tour
          </button>
          <span className="hint">6 steps · ~90 seconds · runs fully offline</span>
        </div>
        {err && <p className="error">{err}</p>}
      </section>

      <div className="entry-grid" data-tour="entry-cards">
        {modes.filter((m) => m.id === "A" || m.id === "B").map((m) => card(m, m.id === "A" ? "a" : "b"))}
      </div>

      <div className="link-row" style={{ marginTop: 16 }}>
        <Link className="pill-link" to="/validation">
          <ShieldCheck size={14} /> Validation statistics
        </Link>
        <Link className="pill-link" to="/manifest">
          <Layers size={14} /> Data sources & manifest
        </Link>
        <button className="pill-link" onClick={startTour} style={{ cursor: "pointer" }}>
          <PlayCircle size={14} /> Replay tutorial
        </button>
      </div>

      <div className="card" style={{ marginTop: 4 }}>
        <span className="section-kicker">How it works</span>
        <h2 style={{ marginTop: 6 }}>Seven stages from raw variants to an explained ranking</h2>
        <div className="pipeline" style={{ marginTop: 12 }}>
          {PIPELINE.map((s) => (
            <div className="pipe-step" key={s.n}>
              <span className="n">{s.n}</span>
              <h4>{s.h}</h4>
              <p>{s.p}</p>
              {facts[s.n] && <div className="fact">{facts[s.n]}</div>}
            </div>
          ))}
        </div>
        <p className="muted" style={{ marginTop: 12 }}>
          Numbers shown are from the real Mode A run — see{" "}
          <Link to="/validation">validation</Link> for enrichment statistics and
          <Link to="/manifest"> provenance</Link>.
        </p>
      </div>

      <div className="grid2">
        <ModePanel modes={modes} onStarted={onStarted} />        <div className="card">
          <h2>
            <ListOrdered size={15} style={{ marginRight: 7, verticalAlign: -2 }} />
            What you get
          </h2>
          <p>
            Every ranked variant carries an <b>evidence chain</b>:
          </p>
          <p className="mono" style={{ lineHeight: 2 }}>
            variant → active MM.1S element → motif disruption → target gene →
            myeloma relevance → explanation
          </p>
          <p className="muted">
            Mode A runs a real GIAB whole-genome VCF (a healthy donor — the funnel
            demonstrates specificity, not cancer discovery). Mode B fine-maps
            published myeloma GWAS loci with LD proxies. All numbers come from the
            datasets listed in the manifest (hash{" "}
            <span className="mono">62ead693</span>).
          </p>
          <p className="muted">
            <b>Research use only — not a diagnostic.</b>
          </p>
        </div>
      </div>

      <RunList refreshKey={0} />
    </>
  );
}

/* --------------------------------------------------------------- run page */

function RunPage() {
  const { id } = useParams<{ id: string }>();
  const nav = useNavigate();
  const [job, setJob] = useState<JobStatus | null>(null);
  const [tab, setTab] = useState<"variants" | "funnel">("variants");
  const [tick, setTick] = useState(0);

  const poll = useCallback(() => {
    if (!id) return;
    api
      .run(id)
      .then((j) => setJob(j))
      .catch((e) => {
        // Deep link to a demo on a cold server: register the precomputed
        // demo bundle first, then read the job again.
        const status = (e as { status?: number })?.status;
        if (status === 404 && /^demo-[AB]$/.test(id)) {
          api
            .demo(id.slice(5) as "A" | "B")
            .then(() => api.run(id))
            .then((j) => setJob(j))
            .catch(() => undefined);
        }
      });
  }, [id]);

  useEffect(poll, [poll, tick]);
  usePoll(
    job && (job.status === "running" || job.status === "queued") ? id : null,
    () => setTick((t) => t + 1),
  );

  if (!job)
    return (
      <div className="card">
        <div className="skel" style={{ height: 18, width: "40%", marginBottom: 12 }} />
        <div className="skel" style={{ height: 12, width: "70%" }} />
        <p className="muted" style={{ marginTop: 12 }}>
          Loading run {id}…
        </p>
      </div>
    );

  const running = job.status === "running" || job.status === "queued";

  return (
    <>
      <div className="ctx-bar">
        <Link to="/" className="ghost btn" style={{ textDecoration: "none" }}>
          ← Home
        </Link>
        <div className="ctx-title">
          <span className="mono" style={{ fontSize: 15 }}>
            {job.run_id}
          </span>
          <span className={`badge ${job.status === "done" ? "ok" : job.status === "error" ? "warn" : ""}`}>
            {job.status}
          </span>
          <span className="badge accent">mode {job.mode}</span>
          {job.total_ranked ? (
            <span className="badge mono">{job.total_ranked.toLocaleString()} ranked</span>
          ) : null}
        </div>
        <div className="spacer" />
        {job.result_available && (
          <a className="btn secondary" href={api.exportCsvUrl(job.run_id)}>
            <Download size={14} />
            Export CSV
          </a>
        )}
      </div>

      <div className="card" style={{ paddingBottom: running ? 16 : 10 }}>
        <p className="muted" style={{ margin: 0 }}>
          {job.input_description}
        </p>
        {running && (
          <div style={{ marginTop: 10 }}>
            <div className="progress">
              <div style={{ width: `${Math.round((job.fraction ?? 0) * 100)}%` }} />
            </div>
            <p className="muted" style={{ marginBottom: 0 }}>
              {job.step ?? ""} — {job.message ?? ""}
            </p>
          </div>
        )}
        {job.status === "error" && <p className="error">{job.error}</p>}
        {job.warnings && job.warnings.length > 0 && (
          <details className="notes">
            <summary>{job.warnings.length} data notes</summary>
            <ul>
              {job.warnings.map((w, i) => (
                <li key={i}>{w}</li>
              ))}
            </ul>
          </details>
        )}
      </div>

      {job.result_available && (
        <>
          <div className="tabs">
            <button
              className={`tab ${tab === "variants" ? "active" : ""}`}
              onClick={() => setTab("variants")}
            >
              Ranked variants
              <span className="count">{(job.total_ranked ?? 0).toLocaleString()}</span>
            </button>
            <button
              className={`tab ${tab === "funnel" ? "active" : ""}`}
              onClick={() => setTab("funnel")}
            >
              How we got here
            </button>
          </div>
          {tab === "variants" ? (
            <VariantsTable
              runId={job.run_id}
              mode={job.mode}
              onSelect={(vid) =>
                nav(`/runs/${job.run_id}/variants/${encodeURIComponent(vid)}`)
              }
            />
          ) : (
            <FunnelView
              steps={(job.funnel ?? []) as unknown as FunnelStepDto[]}
              onPickVariant={() => setTab("variants")}
            />
          )}
        </>
      )}
    </>
  );
}

function VariantPage() {
  const { id, vid } = useParams<{ id: string; vid: string }>();
  if (!id || !vid) return null;
  return <VariantDetailCard runId={id} variantId={decodeURIComponent(vid)} />;
}

/* ------------------------------------------------------------------- shell */

export default function App() {
  const location = useLocation();
  const nav = useNavigate();
  const [health, setHealth] = useState<Health | null>(null);
  const [modes, setModes] = useState<ModeInfo[]>([]);
  const [glossary, setGlossary] = useState(false);

  useEffect(() => {
    api.health().then(setHealth).catch(() => undefined);
    api.modes().then(setModes).catch(() => undefined);
  }, []);

  const startTour = useCallback(() => {
    window.dispatchEvent(new CustomEvent("myelovar:tour"));
  }, []);

  const openRanked = useCallback(async () => {
    try {
      const runs = await api.runs();
      const done = runs.find(
        (r) => r.status === "done" && r.result_available && (r.total_ranked ?? 0) > 0,
      );
      if (done) {
        nav(`/runs/${done.run_id}`);
        return;
      }
      const r = await api.demo("A");
      nav(`/runs/${r.run_id}`);
    } catch {
      toast.error("Ranked list is not available yet.");
    }
  }, [nav]);

  const rankedCount = modes.find((m) => m.id === "A")?.demo_total_ranked;

  return (
    <MotionConfig reducedMotion="user">
      <div className="app">
        <header className="topbar">
          <Link to="/" className="brand" aria-label="Myelovar home">
            <Dna size={17} className="helix" style={{ marginRight: 5, alignSelf: "center" }} />
            Myelo<span>var</span>
          </Link>
          <nav>
            <Link to="/" className={location.pathname === "/" ? "active" : ""}>
              Home
            </Link>
            <button className="navbtn" onClick={openRanked}>
              Ranked list
            </button>
            <Link
              to="/validation"
              className={location.pathname === "/validation" ? "active" : ""}
            >
              Validation
            </Link>
            <Link
              to="/manifest"
              className={location.pathname === "/manifest" ? "active" : ""}
            >
              Data sources
            </Link>
          </nav>
          <div className="spacer" />
          <details className="help">
            <summary>
              <CircleHelp size={15} />
              Help
            </summary>
            <div className="help-menu">
              <button
                onClick={(e) => {
                  startTour();
                  (e.currentTarget.closest("details") as HTMLElement | null)?.removeAttribute("open");
                }}
              >
                <PlayCircle size={15} /> Replay tutorial
              </button>
              <button
                onClick={(e) => {
                  setGlossary(true);
                  (e.currentTarget.closest("details") as HTMLElement | null)?.removeAttribute("open");
                }}
              >
                <BookOpen size={15} /> Glossary
              </button>
              <div className="sep" />
              <Link to="/validation">
                <ShieldCheck size={15} /> Validation statistics
              </Link>
              <div className="note">Shortcuts: <span className="kbd">/</span> search · <span className="kbd">Esc</span> clear filters</div>
            </div>
          </details>
          {health && (
            <>
              <span className={`badge ${health.data_ready ? "ok" : "warn"}`}>
                {health.data_ready ? "data ready" : "data missing"}
              </span>
              <span className="badge mono">
                v{health.version} · {health.manifest_hash.slice(0, 8) || "no hash"}
              </span>
            </>
          )}
        </header>

        <main>
          <AnimatePresence mode="wait" initial={false}>
            <Routes location={location} key={location.pathname}>
              <Route
                path="/"
                element={
                  <Page>
                    <HomePage modes={modes} startTour={startTour} />
                  </Page>
                }
              />
              <Route
                path="/runs/:id"
                element={
                  <Page>
                    <RunPage />
                  </Page>
                }
              />
              <Route
                path="/runs/:id/variants/:vid"
                element={
                  <Page>
                    <VariantPage />
                  </Page>
                }
              />
              <Route
                path="/validation"
                element={
                  <Page>
                    <ValidationView />
                  </Page>
                }
              />
              <Route
                path="/manifest"
                element={
                  <Page>
                    <ManifestView />
                  </Page>
                }
              />
              <Route
                path="/index.html"
                element={
                  <Page>
                    <HomePage modes={modes} startTour={startTour} />
                  </Page>
                }
              />
              <Route
                path="*"
                element={
                  <Page>
                    <HomePage modes={modes} startTour={startTour} />
                  </Page>
                }
              />
            </Routes>
          </AnimatePresence>
        </main>

        <footer className="app-footer">
          <span className="disclaimer">Research use only — not a diagnostic.</span>
          <span>
            Manifest <span className="mono">{health?.manifest_hash?.slice(0, 16) || "—"}</span>
          </span>
          <button
            className="ghost small"
            onClick={() => window.dispatchEvent(new CustomEvent("myelovar:tour"))}
          >
            <Sparkles size={13} /> Replay tutorial
          </button>
          <button className="ghost small" onClick={() => setGlossary(true)}>
            <BookOpen size={13} /> Glossary
          </button>
          <button className="ghost small" onClick={() => window.location.reload()}>
            <RefreshCw size={13} /> Reload data
          </button>
        </footer>

        <Glossary open={glossary} onClose={() => setGlossary(false)} />
        <Tour rankedCount={rankedCount ?? undefined} />
        <Toaster
          position="bottom-right"
          toastOptions={{
            style: {
              background: "#121826",
              color: "#e7ecf9",
              border: "1px solid #27324c",
              fontSize: "13px",
              fontFamily: "inherit",
            },
          }}
        />
      </div>
    </MotionConfig>
  );
}
