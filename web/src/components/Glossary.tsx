import { motion } from "motion/react";
import { X, BookOpen } from "lucide-react";

const TERMS: [string, string][] = [
  [
    "Regulatory (non-coding) variant",
    "A DNA change outside protein-coding sequence. It can still matter by turning genes up or down — that is what Myelovar prioritizes.",
  ],
  [
    "MM.1S",
    "A multiple-myeloma cell line. All 18 functional tracks (open chromatin, histone marks, TF peaks, super-enhancers) were measured in MM.1S, so “active” always means active in myeloma cells.",
  ],
  [
    "cCRE",
    "candidate cis-Regulatory Element — an ENCODE-defined DNA region with regulatory activity (enhancer, promoter, or open region).",
  ],
  [
    "TF motif / PWM",
    "A transcription-factor binding site described as a position-weight matrix (JASPAR). The pipeline scores ±30 bp of real reference DNA around each variant on both strands.",
  ],
  [
    "Motif broken / created / changed",
    "broken = the alt allele loses ≥20 points from a site that scored ≥80% of the motif’s maximum; created = the mirror image; changed = a real shift below those bars.",
  ],
  [
    "ABC model (activity-by-contact)",
    "Links an enhancer to its target gene by combining element activity with 3D contact frequency. No plasma-cell contact map exists, so a B-cell dataset is used and labelled as a proxy everywhere.",
  ],
  [
    "Link confidence",
    "high = ABC contact link · medium = within 100 kb of the gene · low = nearest-gene fallback. Always shown next to the gene so you can weigh it.",
  ],
  [
    "Context score (0–1)",
    "Weighted blend of which regulatory evidence covers the variant: accessibility 0.30, H3K27ac 0.25, super-enhancer 0.15, H3K4me1/H3K4me3 0.10 each, TF peak 0.10.",
  ],
  [
    "Super-enhancer",
    "A large, exceptionally active enhancer cluster (ROSE algorithm) — often the dominant drivers of cell identity and cancer programs.",
  ],
  [
    "DepMap dependency",
    "CRISPR gene-effect scores across cell lines. “Myeloma-selective” means knocking the gene out hurts myeloma lines more than other lines (mean effect ≤ −0.5 with a selectivity margin).",
  ],
  [
    "GWAS locus / LD proxy",
    "A genomic region associated with disease risk in population studies. LD proxies are variants statistically inherited together with the lead SNP (r² ≥ 0.6, 1000G EUR).",
  ],
  [
    "Score weights",
    "regulatory context 0.25 · motif disruption 0.20 · gene link 0.20 · myeloma relevance 0.25 · GWAS 0.10. Deterministic — same inputs, same ranking (Spearman 0.976 under ±50% weight perturbation).",
  ],
  [
    "Manifest hash",
    "A content hash over every downloaded dataset (62ead693212ebba9, 32 resources). Every result embeds it so a ranked list can be traced to exact inputs.",
  ],
];

export default function Glossary({ open, onClose }: { open: boolean; onClose: () => void }) {
  if (!open) return null;
  return (
    <>
      <motion.div
        className="drawer-veil"
        initial={{ opacity: 0 }}
        animate={{ opacity: 1 }}
        exit={{ opacity: 0 }}
        onClick={onClose}
      />
      <motion.aside
        className="drawer"
        initial={{ x: 40, opacity: 0 }}
        animate={{ x: 0, opacity: 1 }}
        exit={{ x: 40, opacity: 0 }}
        transition={{ duration: 0.22, ease: [0.22, 1, 0.36, 1] }}
        role="dialog"
        aria-label="Glossary"
      >
        <header>
          <BookOpen size={17} color="var(--accent)" />
          <h2>Glossary</h2>
          <button className="ghost" onClick={onClose} aria-label="Close glossary">
            <X size={16} />
          </button>
        </header>
        <div className="body">
          <dl style={{ margin: 0 }}>
            {TERMS.map(([t, d]) => (
              <div className="term" key={t}>
                <dt>{t}</dt>
                <dd>{d}</dd>
              </div>
            ))}
          </dl>
        </div>
      </motion.aside>
    </>
  );
}
