"""
fetch_data.py — Kriegstein Lab Demo
Downloads and processes MAST differential expression data from the lab's GitHub
(cortical-arealization repo) to show PFC-enriched marker genes by cell type and stage.

Source: github.com/kriegstein-lab/cortical-arealization
File: neocortex/markers/tbls_in/markers_by_celltype_indiv_all_MAST.txt
"""

import urllib.request, csv, os
from collections import defaultdict

OUT = os.path.dirname(__file__)
URL = ("https://raw.githubusercontent.com/kriegstein-lab/cortical-arealization/"
       "master/neocortex/markers/tbls_in/markers_by_celltype_indiv_all_MAST.txt")

print("Downloading MAST marker file from Kriegstein lab GitHub …")
with urllib.request.urlopen(URL, timeout=60) as r:
    raw = r.read().decode("utf-8")
print(f"  Downloaded {len(raw):,} bytes")

rows = list(csv.DictReader(raw.splitlines(), delimiter="\t"))
print(f"  {len(rows):,} rows")

# ── Figure 1: Top PFC-enriched genes per cell type (early stage) ──────────────
# Collect max logFC per gene per cell type for early PFC markers
pfc_early = defaultdict(lambda: defaultdict(lambda: (-99, 1.0)))
for r in rows:
    if r["area"] != "pfc" or r["stage"] != "early":
        continue
    ct = r["cell_type"]
    gene = r["gene"]
    try:
        fc = float(r["avg_logfc"])
        padj = float(r["p_val_adj"])
    except (ValueError, KeyError):
        continue
    if padj < 0.05 and fc > pfc_early[ct][gene][0]:
        pfc_early[ct][gene] = (fc, padj)

CELL_TYPES = ["rg", "ipc", "neuron"]
TOP_N = 8

with open(os.path.join(OUT, "pfc_markers_early.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["cell_type", "gene", "avg_logfc", "p_val_adj"])
    for ct in CELL_TYPES:
        ranked = sorted(pfc_early[ct].items(), key=lambda x: -x[1][0])
        for gene, (fc, padj) in ranked[:TOP_N]:
            # Skip non-coding / MT genes for cleaner figure
            if gene.startswith(("MT-", "AC0", "MIR", "LINC")):
                continue
            w.writerow([ct, gene, f"{fc:.4f}", f"{padj:.2e}"])
print("Wrote pfc_markers_early.tsv")

# ── Figure 2: Stage dynamics for selected focus genes ─────────────────────────
# HOPX: outer radial glia (oRG) marker — well-known human cortex expansion gene
# CLU: Clusterin, PFC-enriched RG early, decreases over time
# NEUROG1: proneural TF, neuron-specific PFC enrichment
# STK17A: neuron-enriched PFC marker
FOCUS_GENES = {
    "rg":     ["HOPX", "CLU", "LAMP5", "ANXA2"],
    "neuron": ["NEUROG1", "STK17A", "MDK", "PPP1R17"],
}
STAGES = ["early", "mid", "late"]

# Collect best logFC per gene/ct/stage
stage_data = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: None)))
for r in rows:
    if r["area"] != "pfc":
        continue
    ct = r["cell_type"]
    gene = r["gene"]
    stage = r["stage"]
    try:
        fc = float(r["avg_logfc"])
        padj = float(r["p_val_adj"])
    except (ValueError, KeyError):
        continue
    if padj >= 0.05:
        continue
    cur = stage_data[ct][gene][stage]
    if cur is None or fc > cur:
        stage_data[ct][gene][stage] = fc

with open(os.path.join(OUT, "stage_dynamics.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["cell_type", "gene", "stage", "avg_logfc"])
    for ct, genes in FOCUS_GENES.items():
        for gene in genes:
            for stage in STAGES:
                val = stage_data[ct][gene][stage]
                if val is not None:
                    w.writerow([ct, gene, stage, f"{val:.4f}"])
print("Wrote stage_dynamics.tsv")

print("\nSummary:")
for ct, genes in FOCUS_GENES.items():
    print(f"\n  {ct.upper()}:")
    for gene in genes:
        vals = [(s, stage_data[ct][gene][s]) for s in STAGES if stage_data[ct][gene][s] is not None]
        print(f"    {gene}: " + ", ".join(f"{s}={v:.3f}" for s, v in vals))
