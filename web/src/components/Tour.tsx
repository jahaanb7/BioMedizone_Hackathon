/**
 * First-visit guided tutorial.
 *
 * The tour is a sequence of *segments*, each bound to one route: joyride runs
 * one instance per segment (steps on the same page), and the controller
 * navigates between segments, waits for the demo payload to be registered and
 * for every step target to exist in the DOM before starting the instance.
 *
 * Events consumed (react-joyride v3): TOUR_END (advance on action "complete",
 * abort on close/skip), STEP_BEFORE (per-step side effects, e.g. applying a
 * filter), TARGET_NOT_FOUND (graceful abort).
 *
 * Auto-starts once per browser (localStorage flag), suppressed with ?tour=0.
 * Replayable any time via the Help menu or the hero button
 * (window event "myelovar:tour").
 */
import { useCallback, useEffect, useRef, useState } from "react";
import { useLocation, useNavigate } from "react-router-dom";
import { Joyride, EVENTS, type Step } from "react-joyride";
import toast from "react-hot-toast";
import { api } from "../api/client";

const FLAG = "myelovar.onboarded";

type StepDef = {
  target: string;
  title: string;
  content: string;
  placement?: Step["placement"];
  /** side effect run when this step first shows (index within its segment) */
  onShow?: () => void;
};

type Seg = {
  key: string;
  /** where this segment lives */
  path: () => string;
  match: (p: string) => boolean;
  /** ensure the demo run is registered server-side before navigating */
  needsDemo?: boolean;
  steps: StepDef[];
};

function applyBrokenFilter() {
  window.dispatchEvent(
    new CustomEvent("myelovar:filter", { detail: { motif_effect: "broken" } }),
  );
}

function firstRowVariant(): string | null {
  const el = document.querySelector('[data-tour="first-row"]');
  return el?.getAttribute("data-vid") ?? null;
}

let demoAEnsured: Promise<void> | null = null;
function ensureDemoA(): Promise<void> {
  if (!demoAEnsured) {
    demoAEnsured = api
      .demo("A")
      .then(() => undefined)
      .catch((e) => {
        demoAEnsured = null;
        throw e;
      });
  }
  return demoAEnsured;
}

export default function Tour({ rankedCount }: { rankedCount?: number }) {
  const nav = useNavigate();
  const loc = useLocation();
  const [segIdx, setSegIdx] = useState<number | null>(null);
  const [ready, setReady] = useState(false);
  const vidRef = useRef<string | null>(null);
  const shownEffects = useRef<Set<number>>(new Set());

  const countLabel =
    rankedCount && rankedCount > 0 ? rankedCount.toLocaleString() : "the ranked";

  const segments: Seg[] = [
    {
      key: "home",
      path: () => "/",
      match: (p) => p === "/" || p === "/index.html",
      steps: [
        {
          target: '[data-tour="hero"]',
          title: "Welcome to Myelovar",
          content:
            "Myelovar ranks non-coding variants for multiple myeloma — every variant carries a reviewable evidence chain. This is where you start.",
          placement: "bottom",
        },
        {
          target: '[data-tour="entry-cards"]',
          title: "Start from a ranked list",
          content:
            "Two precomputed, fully offline analyses: a whole-genome scan and a fine-mapping of known GWAS loci. Open either in one click.",
          placement: "top",
        },
      ],
    },
    {
      key: "ranked",
      path: () => "/runs/demo-A",
      match: (p) => p === "/runs/demo-A",
      needsDemo: true,
      steps: [
        {
          target: '[data-tour="stats-strip"]',
          title: "This is the ranked list",
          content: `${countLabel} variants ordered by evidence, best first. Every number below is from a real run — click any row to open its evidence chain.`,
          placement: "bottom",
        },
        {
          target: '[data-tour="filter-motif"]',
          title: "Filter by evidence",
          content:
            "Let's keep only variants that break a transcription-factor binding site — watch the match count update.",
          placement: "bottom",
          onShow: applyBrokenFilter,
        },
        {
          target: '[data-tour="first-row"]',
          title: "Open the top row",
          content:
            "Rows are ranked by score. Open one to see why it ranked: element → motif → gene → myeloma relevance.",
          placement: "top",
        },
      ],
    },
    {
      key: "detail",
      path: () => `/runs/demo-A/variants/${encodeURIComponent(vidRef.current ?? "")}`,
      match: (p) => p.startsWith("/runs/demo-A/variants/"),
      steps: [
        {
          target: '[data-tour="chain"]',
          title: "The evidence chain",
          content:
            "Each variant shows the active element, the motif it disrupts, the target gene and how it was linked, plus its myeloma relevance — a chain you can challenge link by link.",
          placement: "left",
        },
        {
          target: '[data-tour="locus"]',
          title: "Genomic context",
          content:
            "All overlapping MM.1S tracks and nearby genes around this variant — evidence lives in regions, not single positions.",
          placement: "top",
        },
      ],
    },
    {
      key: "validation",
      path: () => "/validation",
      match: (p) => p === "/validation",
      steps: [
        {
          target: '[data-tour="validation-stats"]',
          title: "Why trust the ranking",
          content:
            "Enrichment statistics with independence flags: the cCRE test is independent of scoring, while the known-gene test is labelled partially circular and weighted down. Replay this tour anytime from Help.",
          placement: "bottom",
        },
      ],
    },
  ];

  const seg = segIdx !== null ? segments[segIdx] : null;

  /* ---- start / stop ------------------------------------------------- */
  const stop = useCallback((finished: boolean) => {
    setSegIdx(null);
    setReady(false);
    localStorage.setItem(FLAG, "1");
    if (finished)
      toast.success("Tour finished — replay anytime: Help ▸ Replay tutorial", {
        duration: 4000,
      });
  }, []);

  const start = useCallback(() => {
    setReady(false);
    setSegIdx(0);
  }, []);

  useEffect(() => {
    const onMsg = () => start();
    window.addEventListener("myelovar:tour", onMsg);
    return () => window.removeEventListener("myelovar:tour", onMsg);
  }, [start]);

  /* auto-start on first visit */
  useEffect(() => {
    if (localStorage.getItem(FLAG)) return;
    const q = new URLSearchParams(window.location.search);
    if (q.get("tour") === "0") {
      localStorage.setItem(FLAG, "1");
      return;
    }
    const t = setTimeout(() => start(), 900);
    return () => clearTimeout(t);
  }, [start]);

  /* ---- segment controller ------------------------------------------- */
  useEffect(() => {
    if (segIdx === null) return;
    let cancelled = false;
    const s = segments[segIdx];
    shownEffects.current = new Set();
    setReady(false);

    const boot = async () => {
      try {
        if (s.needsDemo) await ensureDemoA();
      } catch {
        if (!cancelled) {
          toast.error("Could not load the precomputed demo — tour stopped.");
          stop(false);
        }
        return;
      }
      if (cancelled) return;
      const want = s.path();
      if (loc.pathname !== want) {
        nav(want);
        return; // effect re-runs on location change
      }
      // wait until every target exists in the DOM
      let tries = 0;
      const iv = window.setInterval(() => {
        if (cancelled) return window.clearInterval(iv);
        const all = s.steps.every((st) => document.querySelector(st.target));
        if (all || ++tries > 50) {
          window.clearInterval(iv);
          if (!all) {
            toast("Tour target not found — skipped.", { icon: "⚠️" });
            stop(false);
          } else {
            setReady(true);
          }
        }
      }, 200);
    };
    void boot();
    return () => {
      cancelled = true;
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [segIdx, loc.pathname]);

  /* ---- joyride event handling --------------------------------------- */
  const onEvent = useCallback(
    (data: {
      type: string;
      action: string;
      index: number;
    }) => {
      if (!seg) return;
      if (data.type === EVENTS.STEP_BEFORE) {
        const def = seg.steps[data.index];
        if (def?.onShow && !shownEffects.current.has(data.index)) {
          shownEffects.current.add(data.index);
          def.onShow();
        }
        return;
      }
      if (data.type === EVENTS.TARGET_NOT_FOUND) {
        toast("Tour interrupted — replay anytime from Help.", { icon: "⚠️" });
        stop(false);
        return;
      }
      if (data.type === EVENTS.TOUR_END) {
        if (data.action !== "complete") {
          stop(false); // user closed or skipped
          return;
        }
        if (seg.key === "ranked") {
          vidRef.current = firstRowVariant();
          if (!vidRef.current) {
            stop(true);
            return;
          }
        }
        if (segIdx !== null && segIdx + 1 < segments.length) {
          setReady(false);
          setSegIdx(segIdx + 1);
        } else {
          stop(true);
        }
      }
    },
    [seg, segIdx, stop],
  );

  if (segIdx === null || !seg) return null;

  const steps: Step[] = seg.steps.map((s) => ({
    target: s.target,
    title: s.title,
    content: s.content,
    placement: s.placement ?? "bottom",
  }));

  return ready ? (
    <Joyride
      key={seg.key}
      run
      continuous
      steps={steps}
      onEvent={onEvent as never}
      options={{
        showProgress: true,
        primaryColor: "#5b9dff",
        zIndex: 10000,
        closeButtonAction: "skip",
      }}
      locale={{ back: "Back", next: "Next", close: "×", skip: "Skip tour" }}
    />
  ) : null;
}
