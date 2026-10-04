# PROMPT — Myelovar web app: semi-final product polish, ranked-list-first UX, researcher onboarding

> **Mission:** turn the working Myelovar web app from a functional prototype into a
> **semi-final, professional product** for doctors and researchers. Three pillars:
> **(1)** visual craft using external UI/graphic/motion libraries (reactbits.dev,
> uiverse.io, motion, lucide, …) so it stops reading like AI slop; **(2)** usability —
> the site's whole point is the **ranked variant list**, and today it is hard to find
> and hard to understand; **(3)** an **initial guided tutorial** that teaches a
> researcher how to use the site on first visit.
>
> **Overriding constraint, in every step: everything must stay DEMOABLE.** A demo of
> Demo A → ranked list → filter → evidence chain → locus → validation must run
> end-to-end in under 3 minutes, offline, against the precomputed demo bundle, with no
> console errors. No step of this redesign may break that path.

---

## 0. Context you must not re-discover

- **Stack:** Vite 5 + React 18 + TypeScript **strict** (`npm run build` = `tsc --noEmit && vite build`),
  React Router 6, `@tanstack/react-table` v8, `plotly.js-dist-min`, plain CSS in
  `web/src/styles.css`, typed API client in `web/src/api/client.ts` (types generated
  from FastAPI OpenAPI into `src/api/schema.d.ts`).
- **Components:** `web/src/components/{FunnelView,VariantsTable,LocusPlot,VariantDetailCard,ValidationView,ManifestView}.tsx`,
  shell/routes in `web/src/App.tsx`.
- **Backend (do not change behavior):** FastAPI at `127.0.0.1:8000`, SPA served from
  `web/dist`. Endpoints: `/api/health`, `/api/modes`, `/api/demo/{A,B}`, `/api/runs`,
  `/api/runs/{id}`, `/runs/{id}/funnel`, `/runs/{id}/variants?min_score&gene&tf&chrom&confidence&in_super_enhancer&locus_id&motif_effect&sort&order&limit≤500&offset`,
  `/runs/{id}/variants/{vid}`, `/runs/{id}/locus`, `/runs/{id}/export.csv`,
  `/api/genes/{symbol}`, `/api/validation`, `/api/validation/figures/{name}`, `/api/manifest`.
  **The redesign must consume only these existing endpoints.**
- **Real numbers to use verbatim (never fake or round differently):** Mode A demo =
  **48,113** ranked of 3,935,734 (common-pruned 1,137,368; in-active-elements 48,113;
  **169 broken / 140 created** motifs; links ABC=357 / distance=46,030 / nearest=2,083;
  18 tracks; 10 PWMs). Mode B demo = **259,211** ranked of 1,863,046 (broken=933 /
  created=1,147; 2,038 REF/FASTA skips). Validation = cCRE **OR 3.85, p=1.6e-78**
  (independent); known-gene **OR 5.46, p=2.4e-12**, flagged **partially circular**;
  weight sweep **Spearman 0.976**; Mode B lead recovery **1/8 in top-3, p=0.007**.
  Manifest hash **62ead693212ebba9**, 32 resources, 31 ok / 1 failed.
- **Audience:** clinicians and bench researchers. They do **not** upload VCFs in a
  demo; they want the ranked list, evidence chains, and a reason to trust the numbers.
  Keep the "Research use only — not a diagnostic" notice visible.
- **Environment:** macOS arm64, Node 24 / npm 11, ~4–5 GiB free disk (do NOT install
  OS-level tooling), API restart needed after rebuilding `web/dist`. Do not commit.

---

## 1. Information architecture — the ranked list is the product

**Problem:** the landing page is an "Analyze" VCF-upload form; the ranked list is
buried three interactions deep (Demo A → run page → "Ranked variants" tab), and the
funnel jargon (s1…s6) is the first thing shown.

**Changes:**

1. **New home (`/`)** — stop leading with the upload form:
   - **Hero:** one sentence of what this is ("Rare non-coding variants in multiple
     myeloma, ranked with a reviewable evidence chain"), plus **two large entry
     cards** carrying the real counts:
     - **"Genome scan — Mode A · 48,113 ranked variants"** → CTA **"Open ranked list"**
     - **"GWAS loci — Mode B · 259,211 ranked variants"** → CTA **"Open ranked list"**
   - Secondary row: "How it works" (funnel explainer), "Validation", "Data sources
     (manifest)", "Run your own VCF" (the existing Analyze panel, demoted but intact).
   - First-visit tutorial starts from here (§4).
2. **Run page (`/runs/{id}`)** — make **"Ranked variants" the default tab**; move the
   funnel to a secondary tab/section renamed **"How we got here"**, with each stage
   re-labelled in plain language (see §5 copy).
3. **Global context bar** on run pages: run name, mode badge, total ranked, active
   filter chips, CSV export, breadcrumb back home. The user must never wonder "where
   am I / what list am I looking at".
4. **Nav:** Home · Ranked · How it works · Validation · Manifest (+ Help menu with
   "Replay tutorial" and "Glossary"). Preserve deep links `/runs/{id}`,
   `/runs/{id}/variants/{vid}`, `/validation`, `/manifest` (bookmarkable).

---

## 2. Visual design — professional, not AI slop

**Direction:** a *precision scientific instrument* — restrained, high-contrast,
editorial typography, real hierarchy, one accent color. Explicitly banned AI-slop
tells: purple/blue gradient heroes, emoji as icons, glass-everything, rainbow
gradients, generic "Dashboard" card grids with no hierarchy, lorem/placeholder text,
default Plotly/Chart.js styling, three evenly spaced feature columns.

**Design tokens** (CSS custom properties, e.g. in `web/src/styles.css` under `:root`
and wired into Tailwind `@theme` per §3): ink/surface/border/elevated color ramps,
type scale (display / h1–h3 / body / small / mono), spacing scale, radii, shadows,
motion durations + easing curves, semantic colors for evidence
(broken=red-amber, created=green, changed=blue, none=grey; link confidence
high/med/low; evidence dots for DepMap/known/GWAS).

**Type:** self-hosted `Inter` (UI) + `IBM Plex Mono` (coordinates, rsIDs, scores)
via `@fontsource-variable/*` npm packages — **no Google Fonts CDN at runtime** (the
demo must work offline).

**Icons:** `lucide-react` only; consistent size/stroke; never emoji.

**Key screens to redesign:** home hero, ranked table, variant detail (evidence chain
should look like a designed stepped diagram, not a bullet list), funnel (horizontal
stage bars with before → after counts and hoverable removal reasons), locus plot
(Plotly restyled to tokens: unified track colors, no default margins/fonts), gene
card, validation (stat cards with large OR/p values + the circularity flag styled as
a deliberate callout, not an apology), manifest (status chips, hash in mono).

---

## 3. External libraries — pick, install, attribute

**Adopt Tailwind CSS v4 alongside the existing CSS** (do not rewrite `styles.css`
from scratch): `npm i -D tailwindcss @tailwindcss/vite`, `@import "tailwindcss"` +
`@theme { …tokens… }` in the stylesheet, plugin in `vite.config.ts`. This is required
so reactbits.dev and uiverse.io snippets (Tailwind-classed) drop in cleanly. Existing
CSS keeps working; new/edited components use tokens.

Install (all npm, all MIT/ISC-compatible; record attributions in
`web/THIRD_PARTY_NOTES.md`):

| Library | Use |
|---|---|
| `motion` (Framer Motion) | **primary motion engine**: route/ tab transitions, list stagger, layout animations, count-ups, shared-element open from table row → detail |
| `react-joyride` | guided tutorial (steps, targets, progress, skip/resume, event hooks) |
| `lucide-react` | icons |
| `@fontsource-variable/inter`, `@fontsource-variable/ibm-plex-mono` | offline fonts |
| `clsx` | conditional class names |
| `react-hot-toast` | feedback (CSV export started, filter reset, tour finished) |
| `tailwind-merge` | token-safe class merging when composing |

**reactbits.dev** — hand-pick 2–4 components that fit a scientific tool, e.g.
`SplitText` or `BlurText` (hero), `CountUp` (funnel/validation stats), `Shimmer`/
skeleton, `AnimatedList`. Use the npm-installable form where offered, otherwise copy
the source into `web/src/components/ui/` with an attribution comment. Max ~4 — no
component zoo.

**uiverse.io** — hand-pick 2–4 *small* elements: skeleton-shimmer loading rows,
toggle switch (dark mode or table density), tooltip/question-mark chip, stepper
dots. Vendor the markup+CSS into `web/src/components/ui/` with source URL +
license comment (**uiverse is CC BY 4.0 — attribution is mandatory**).

**Rules:** one motion engine (`motion`); respect `prefers-reduced-motion` globally
(a `useReducedMotion` guard that disables all non-essential animation); lazy-load
Plotly (`React.lazy`) so the table route stays light; no runtime CDN fetches of any
kind (fonts, icons, scripts, CSS).

---

## 4. Onboarding tutorial — first-visit guided tour

Library: **react-joyride**. Behavior:

- On first visit (flag `localStorage["myelovar.onboarded"]`), auto-start after data
  loads. Always replayable: header **Help ▸ Replay tutorial**. Auto-skip entirely when
  `?tour=0` in the URL (needed so verification scripts and tests are deterministic).
- Persistent **Glossary** panel (Help menu): cCRE, ABC model, PWM/motif, LD proxy,
  DepMap, "link confidence", "context score" — 1–2 plain sentences each.
- Skippable, back/next, progress dots, focus management, no trapping; every step's
  spotlight must actually exist on the page (verify each target in the browser).

**The six steps (exact copy intent):**

1. **Home** — "Myelovar ranks non-coding variants for myeloma. Start with a ranked
   list below." (spotlight: entry cards)
2. **Ranked list** — navigate to Mode A ranked tab; spotlight the table header +
   count: "This is the ranked list — 48,113 variants ordered by evidence, best first."
3. **Filtering** — spotlight the motif-effect filter; the tour itself applies
   **`motif_effect = broken`** and the count animates to **169**:
   "Filter to variants that break a transcription-factor binding site."
4. **Evidence chain** — click the top row; spotlight the 5-step chain: "Every variant
   shows *why* it ranked: element → motif → gene → myeloma relevance → explanation."
5. **Locus** — spotlight the locus plot: "See the variant in genomic context — all
   18 tracks and nearby genes."
6. **Trust** — route to `/validation`: "Ranks are validated: cCRE enrichment
   OR 3.85 (p=1.6e-78); circular checks are labelled." End state: toast + set the
   onboarded flag; link "Open Mode B list".

---

## 5. Usability — understand what's happening, at all times

- **Plain-language stage labels** (keep `s1…s6` only as small mono badges):
  "Remove common variants", "Keep variants in myeloma-active DNA", "Scan for motif
  damage", "Link to target genes", "Score myeloma relevance", "Rank & explain".
  Show `before → after` + *why N removed* on hover/expand (data already in
  `/funnel`).
- **Explain the score:** a "?" on the score column → tooltip with the five weighted
  components (regulatory_context 0.25, motif 0.20, link 0.20, relevance 0.25,
  gwas 0.10) and a mini bar breakdown per row on hover.
- **Evidence chips everywhere:** motif effect (color-coded), link confidence
  (high/med/low), super-enhancer ★, and three dots for DepMap/known/GWAS — in the
  table *and* the detail card, same visual language.
- **Filter bar:** search (gene/TF), dropdowns (motif effect, confidence, chrom),
  min-score slider, all existing API filters surfaced; **active filters as removable
  chips + live match count + "Reset"**; filter state mirrored to the URL query string
  (shareable/bookmarkable views).
- **States:** skeletons (uiverse shimmer) for table/detail/locus loading; friendly
  empty state with reset CTA; error state with retry; disable buttons while in
  flight; toast on CSV export.
- **Microcopy pass** across every screen — no naked jargon in headings. Example
  funnel note: "scanned 48,113 sites with 10 TF motifs on both strands; 169 broken,
  140 created (0 reference mismatches)".
- **Keyboard:** table rows focusable, Enter opens detail, `/` focuses search,
  `Esc` clears filters.

---

## 6. Ranked-list screen spec (the money screen)

Table columns: **Rank** (mono), **Variant** (`chr:pos REF>ALT`, mono, click to copy),
**Score** (number + inline bar), **Motif effect** (chip), **Top TF**, **Target
genes** (first 2 + "+n", link-confidence chip), **Context** (0–1 mini gauge),
**Evidence** (DepMap/known/GWAS dots), **Super-enhancer** ★.

- Server-side paging (`limit≤500`) and sorting (`sort=score|rank|pos`, `order`) via
  existing API — do not fetch 48k/259k rows client-side.
- Row click → `/runs/{id}/variants/{vid}` detail (5-step chain diagram, explanation
  callout, gene card, locus plot).
- Above table: **stats strip** — total ranked, broken motifs, ABC-linked,
  known-gene hits (from data already in the run/funnel payload; no new endpoints).
- Virtualized/paged with no jank at Mode B's 259,211 rows.

---

## 7. DEMOABLE — acceptance criteria (all must hold at the end)

1. `cd web && npm run build` → **0 tsc errors**, vite build green; API restarted so
   the new `web/dist` is served.
2. `.venv/bin/python -m pytest` → **33 passed** (no backend changes expected).
3. Fresh browser profile at `http://127.0.0.1:8000`: tutorial auto-starts once,
   runs through all 6 steps with valid spotlights, is replayable from Help, and
   `?tour=0` suppresses it. Flag set afterwards; reload does not re-trigger.
4. Home shows both entry cards with the real counts (48,113 / 259,211) and reaches
   the ranked list in **one click**; the ranked tab is the default on run pages.
5. Mode A ranked list: filter `motif_effect=broken` → **169 rows**; sorting by
   `rank` after `score` stays monotonic for ties; pagination, CSV export, URL filter
   state, and Reset all work.
6. Variant detail renders the 5-step evidence chain + explanation + gene card +
   locus plot for at least one top variant in **each** of Demo A and Demo B.
7. Validation page shows all 3 tests with the partially-circular flag and lead
   recovery; Manifest shows 32 resources + hash `62ead693212ebba9`.
8. **Demo A and Demo B buttons both enabled**; the scripted demo path
   (home → tour → ranked A → filter → detail → locus → Mode B → validation →
   manifest) completes in <3 min with **zero console errors** and no network
   requests outside `127.0.0.1:8000` (no CDNs).
9. `prefers-reduced-motion: reduce` disables all non-essential animation.
10. Visual pass at 1280×800 and 1512×980: screenshots of home, ranked, detail,
    locus, validation, manifest reviewed against the "AI slop" ban list in §2.
11. Third-party attributions recorded in `web/THIRD_PARTY_NOTES.md` (incl. uiverse
    CC BY 4.0); README gets a short "Using the site" section linking the tutorial.

---

## 8. Implementation order (each step leaves the demo working)

1. Tokens + Tailwind v4 bridge + fonts + lucide + `clsx` (no behavior change).
2. IA rewire: home page, default ranked tab, nav/help menu, context bar.
3. Ranked-list upgrade: columns/chips/filters/URL state/stats strip.
4. Evidence-chain + funnel restyle (motion: count-ups, bar reveals, row stagger).
5. Tutorial (joyride) + glossary + microcopy pass.
6. Locus/plotly restyle, validation/manifest polish, empty/loading/error states.
7. reactbits/uiverse component drops + route transitions + hover polish.
8. Full acceptance run (§7) + screenshots + attribution notes.

## 9. Constraints & guardrails

- **Do not touch** `src/myelovar/`, `api/` behavior, pipeline outputs, or test
  expectations. Frontend-only (plus docs).
- No new backend endpoints; if a number isn't in an existing payload, show a number
  that is — never invent.
- No runtime CDNs, no analytics, no external font/script hosts.
- Keep bundle sane: Plotly lazy-loaded; target initial-route JS gzip < ~600 KB.
- Disk is tight (~4–5 GiB): don't install OS packages or regenerate bulky data.
- Preserve existing routes and the "Research use only / not a diagnostic" notice.
- Verify by **using the site in the browser**, not just by building.
