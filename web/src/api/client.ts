/**
 * Typed API client. Endpoint types are derived from the generated OpenAPI
 * schema (src/api/schema.d.ts, regenerated with `npm run gen`), with a few
 * inline-object endpoints typed here where the spec inlines them.
 */
import type { components, paths } from "./schema";

export type Health = components["schemas"]["HealthResponse"];
export type ModeInfo = components["schemas"]["ModeInfo"];
export type JobStatus = components["schemas"]["JobStatus"];
export type RunCreateResponse = components["schemas"]["RunCreateResponse"];
export type VariantsPage = components["schemas"]["VariantsPage"];
export type VariantDetail = components["schemas"]["VariantDetail"];
export type LocusResponse = components["schemas"]["LocusResponse"];
export type LocusTrack = components["schemas"]["LocusTrack"];
export type GeneCardResponse = components["schemas"]["GeneCardResponse"];

export interface FunnelStepDto {
  step: string;
  label: string;
  rows_in: number;
  rows_out: number;
  removed: number;
  removed_reasons?: Record<string, number>;
  note?: string;
}

export interface FunnelResponse {
  run_id: string;
  mode: string;
  status: string;
  steps: FunnelStepDto[];
}

export interface ManifestEntry {
  id: string;
  category: string;
  source_url: string;
  accession?: string;
  download_date: string;
  genome_build?: string;
  cell_type?: string;
  assay?: string;
  path: string;
  size_bytes: number;
  md5?: string;
  verification?: string;
  status: string;
  error?: string;
  notes?: string;
}

export interface Manifest {
  generated_at: string;
  genome_build?: string;
  manifest_hash?: string;
  resources: ManifestEntry[];
}

export interface ValidationResults {
  [key: string]: unknown;
}

/** Resolve an OpenAPI path+method response type. */
type Op<P extends keyof paths, M extends keyof paths[P]> = paths[P][M];
type Resp<
  P extends keyof paths,
  M extends keyof paths[P],
> = Op<P, M> extends { responses: { 200: { content: { "application/json": infer S } } } }
  ? S
  : never;

export class ApiError extends Error {
  status: number;
  detail: unknown;
  constructor(status: number, detail: unknown) {
    super(
      typeof detail === "string"
        ? detail
        : JSON.stringify(detail ?? status),
    );
    this.status = status;
    this.detail = detail;
  }
}

async function req<T>(url: string, init?: RequestInit): Promise<T> {
  const res = await fetch(url, init);
  if (!res.ok) {
    let detail: unknown;
    try {
      const body = await res.json();
      detail = (body as { detail?: unknown }).detail ?? body;
    } catch {
      detail = res.statusText;
    }
    throw new ApiError(res.status, detail);
  }
  return (await res.json()) as T;
}

const get = <T>(url: string) => req<T>(url);

export const api = {
  health: () => get<Health>("/api/health"),
  manifest: () => get<Resp<"/api/manifest", "get">>("/api/manifest"),
  modes: () => get<ModeInfo[]>("/api/modes"),
  demo: (mode: "A" | "B") =>
    get<{ run_id: string; mode: string; result: unknown }>(`/api/demo/${mode}`),
  runs: () => get<JobStatus[]>("/api/runs"),
  run: (id: string) => get<JobStatus>(`/api/runs/${id}`),
  funnel: (id: string) => get<FunnelResponse>(`/api/runs/${id}/funnel`),
  variants: (id: string, q: Record<string, string | number | undefined>) => {
    const params = new URLSearchParams();
    for (const [k, v] of Object.entries(q))
      if (v !== undefined && v !== "") params.set(k, String(v));
    return get<VariantsPage>(`/api/runs/${id}/variants?${params}`);
  },
  variant: (id: string, variantId: string) =>
    get<VariantDetail>(
      `/api/runs/${id}/variants/${encodeURIComponent(variantId)}`,
    ),
  locus: (id: string, variantId: string, window = 50_000) =>
    get<LocusResponse>(
      `/api/runs/${id}/locus?variant_id=${encodeURIComponent(variantId)}&window=${window}`,
    ),
  exportCsvUrl: (id: string) => `/api/runs/${id}/export.csv`,
  gene: (symbol: string) => get<GeneCardResponse>(`/api/genes/${symbol}`),
  validation: () => get<ValidationResults>("/api/validation"),
  figureUrl: (name: string) => `/api/validation/figures/${name}`,

  startRun: async (mode: "A" | "B" | "C", file?: File) => {
    const form = new FormData();
    form.set("mode", mode);
    if (file) form.set("file", file);
    return req<RunCreateResponse>("/api/runs", { method: "POST", body: form });
  },
};
