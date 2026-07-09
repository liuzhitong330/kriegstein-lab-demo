"""
make_figures.py — Kriegstein Lab Demo
Figure 1: Top PFC-enriched marker genes by cell type (early dev, lollipop chart)
Figure 2: Stage dynamics of key RG markers (HOPX vs CLU, early→mid→late)
"""

import csv, os, math

OUT = os.path.dirname(__file__)

# ── Load data ──────────────────────────────────────────────────────────────────
markers = []
with open(os.path.join(OUT, "pfc_markers_early.tsv")) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        markers.append((row["cell_type"], row["gene"], float(row["avg_logfc"])))

stage_rows = []
with open(os.path.join(OUT, "stage_dynamics.tsv")) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        stage_rows.append((row["cell_type"], row["gene"], row["stage"], float(row["avg_logfc"])))

# ── Figure 1: Lollipop chart — PFC markers by cell type ───────────────────────
FW, FH = 700, 440
PAD_L = 115
PAD_R = 100
PAD_T = 75
PAD_B = 55
AW = FW - PAD_L - PAD_R
AH = FH - PAD_T - PAD_B

CT_COLOR = {"rg": "#1a5c8a", "ipc": "#27ae60", "neuron": "#c0392b"}
CT_LABEL = {"rg": "Radial Glia", "ipc": "IPC", "neuron": "Neuron"}
CT_ORDER = ["rg", "ipc", "neuron"]

# Group markers by cell type
by_ct = {ct: [(g, fc) for (c, g, fc) in markers if c == ct] for ct in CT_ORDER}

# Build Y layout: group separator + gene rows
rows_y = []  # (label, logfc, color, is_header)
gap_header = 18
gap_gene = 28
gap_between = 14

for ct in CT_ORDER:
    if rows_y:
        rows_y.append(("", 0, "#999", True))  # spacer
    rows_y.append((CT_LABEL[ct], None, CT_COLOR[ct], True))
    for gene, fc in by_ct[ct]:
        rows_y.append((gene, fc, CT_COLOR[ct], False))

# Assign Y positions
ys = []
y = PAD_T
for label, fc, col, is_header in rows_y:
    if label == "" and fc == 0:
        y += gap_between
        continue
    ys.append((y, label, fc, col, is_header))
    y += gap_header if is_header else gap_gene

# Scale X
max_fc = max(fc for _, _, fc in markers) * 1.12
def xpos(fc):
    return PAD_L + fc / max_fc * AW

# Draw
lollipops = ""
for y, label, fc, col, is_header in ys:
    if is_header:
        lollipops += (f'<text x="{PAD_L - 10}" y="{y + 5}" text-anchor="end" '
                      f'font-size="11" fill="{col}" font-weight="700">{label}</text>')
        continue
    x = xpos(fc)
    # stem
    lollipops += (f'<line x1="{PAD_L}" y1="{y}" x2="{x:.1f}" y2="{y}" '
                  f'stroke="{col}" stroke-width="1.5" opacity="0.4"/>')
    # dot
    lollipops += (f'<circle cx="{x:.1f}" cy="{y}" r="6" fill="{col}" opacity="0.85"/>')
    # gene label
    lollipops += (f'<text x="{PAD_L - 8}" y="{y + 4}" text-anchor="end" '
                  f'font-size="10" fill="#333">{label}</text>')
    # fc label
    lollipops += (f'<text x="{x + 9:.1f}" y="{y + 4}" font-size="9" fill="{col}" font-weight="600">'
                  f'{fc:.2f}</text>')

# X axis
tick_vals = [0.0, 0.3, 0.6, 0.9, 1.2]
xticks = ""
for v in tick_vals:
    if v > max_fc:
        break
    tx = xpos(v)
    xticks += (f'<line x1="{tx:.1f}" y1="{PAD_T - 8}" x2="{tx:.1f}" '
               f'y2="{PAD_T + AH}" stroke="#eee" stroke-width="1"/>'
               f'<text x="{tx:.1f}" y="{PAD_T + AH + 16}" text-anchor="middle" '
               f'font-size="9" fill="#888">{v:.1f}</text>')
axis1 = (f'<line x1="{PAD_L}" y1="{PAD_T - 8}" x2="{PAD_L}" '
         f'y2="{PAD_T + AH}" stroke="#ccc" stroke-width="1"/>')

# Legend
leg = ""
lx, ly = FW - PAD_R + 12, PAD_T + 20
for i, ct in enumerate(CT_ORDER):
    leg += (f'<circle cx="{lx + 6}" cy="{ly + i*22}" r="5" fill="{CT_COLOR[ct]}"/>'
            f'<text x="{lx + 16}" y="{ly + i*22 + 4}" font-size="10" fill="{CT_COLOR[ct]}" '
            f'font-weight="600">{CT_LABEL[ct]}</text>')

svg1 = f"""<svg viewBox="0 0 {FW} {FH}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{FW//2}" y="22" text-anchor="middle" font-size="13" font-weight="600" fill="#222">
    Top PFC-Enriched Genes by Cell Type (Early Development)
  </text>
  <text x="{FW//2}" y="40" text-anchor="middle" font-size="10" fill="#666">
    log fold change (PFC vs. other areas) · p_adj &lt; 0.05 · MAST test · Kriegstein lab, cortical-arealization repo
  </text>
  <text x="{FW//2}" y="56" text-anchor="middle" font-size="10" fill="#444">
    HOPX (oRG marker) and ANXA2 enriched in PFC radial glia; NEUROG1 and STK17A in early neurons
  </text>
  <text x="{PAD_L + AW/2:.0f}" y="{PAD_T + AH + 34}" text-anchor="middle"
        font-size="10" fill="#555">Average log fold change (PFC enrichment)</text>
  {axis1}{xticks}{lollipops}{leg}
</svg>"""

with open(os.path.join(OUT, "pfc_markers.svg"), "w") as f:
    f.write(svg1)
print("Wrote pfc_markers.svg")


# ── Figure 2: Stage dynamics — RG markers over time ───────────────────────────
# Focus on RG cell type, 4 genes, early/mid/late
FW2, FH2 = 680, 380
PAD_L2 = 60
PAD_R2 = 30
PAD_T2 = 75
PAD_B2 = 70
AW2 = FW2 - PAD_L2 - PAD_R2
AH2 = FH2 - PAD_T2 - PAD_B2

RG_GENES = ["ANXA2", "HOPX", "CLU", "LAMP5"]
STAGES = ["early", "mid", "late"]
STAGE_LABELS = {"early": "Early", "mid": "Mid", "late": "Late"}
STAGE_COLORS = {"early": "#d35400", "mid": "#8e44ad", "late": "#2980b9"}
GENE_COLORS = {"ANXA2": "#c0392b", "HOPX": "#1a5c8a", "CLU": "#27ae60", "LAMP5": "#e67e22"}

# Build lookup
rg_vals = {g: {} for g in RG_GENES}
for ct, gene, stage, fc in stage_rows:
    if ct == "rg" and gene in RG_GENES:
        rg_vals[gene][stage] = fc

# Layout: 4 gene groups, 3 bars per group
n_genes = len(RG_GENES)
n_stages = 3
group_w = AW2 / n_genes
bar_w = group_w * 0.22
bar_gap = group_w * 0.05
max_fc2 = 1.25

def ypos(fc):
    return PAD_T2 + AH2 - fc / max_fc2 * AH2

bars2 = ""
for gi, gene in enumerate(RG_GENES):
    gx = PAD_L2 + gi * group_w + group_w * 0.1
    col = GENE_COLORS[gene]
    # Gene label under group
    bars2 += (f'<text x="{PAD_L2 + gi*group_w + group_w/2:.1f}" y="{PAD_T2 + AH2 + 20}" '
              f'text-anchor="middle" font-size="11" fill="{col}" font-weight="700">{gene}</text>')
    for si, stage in enumerate(STAGES):
        val = rg_vals[gene].get(stage)
        bx = gx + si * (bar_w + bar_gap)
        sc = STAGE_COLORS[stage]
        if val is not None:
            by = ypos(val)
            bh = PAD_T2 + AH2 - by
            bars2 += (f'<rect x="{bx:.1f}" y="{by:.1f}" width="{bar_w:.1f}" height="{bh:.1f}" '
                      f'fill="{sc}" opacity="0.8" rx="2"/>')
            bars2 += (f'<text x="{bx + bar_w/2:.1f}" y="{by - 5:.1f}" text-anchor="middle" '
                      f'font-size="8.5" fill="{sc}" font-weight="600">{val:.2f}</text>')
        else:
            # No data — dashed box
            bars2 += (f'<rect x="{bx:.1f}" y="{PAD_T2 + AH2 - 20:.1f}" width="{bar_w:.1f}" '
                      f'height="20" fill="none" stroke="{sc}" stroke-dasharray="3,2" opacity="0.4" rx="2"/>')

# Y axis
yticks2 = ""
ybase = PAD_T2 + AH2
for v in [0.0, 0.3, 0.6, 0.9, 1.2]:
    ty = ypos(v)
    yticks2 += (f'<line x1="{PAD_L2}" y1="{ty:.1f}" x2="{PAD_L2 + AW2}" y2="{ty:.1f}" '
                f'stroke="#eee" stroke-width="1"/>'
                f'<text x="{PAD_L2 - 8}" y="{ty + 4:.1f}" text-anchor="end" '
                f'font-size="9" fill="#888">{v:.1f}</text>')
axis2 = (f'<line x1="{PAD_L2}" y1="{PAD_T2}" x2="{PAD_L2}" y2="{ybase}" '
         f'stroke="#ccc" stroke-width="1"/>'
         f'<line x1="{PAD_L2}" y1="{ybase}" x2="{PAD_L2 + AW2}" y2="{ybase}" '
         f'stroke="#ccc" stroke-width="1"/>')

# Stage legend
leg2 = ""
lx2 = PAD_L2 + AW2 * 0.6
ly2 = PAD_T2 + 8
for i, stage in enumerate(STAGES):
    leg2 += (f'<rect x="{lx2 + i*72}" y="{ly2}" width="12" height="12" '
             f'fill="{STAGE_COLORS[stage]}" rx="2" opacity="0.8"/>'
             f'<text x="{lx2 + i*72 + 16}" y="{ly2 + 10}" font-size="10" '
             f'fill="{STAGE_COLORS[stage]}" font-weight="600">{STAGE_LABELS[stage]}</text>')

# Y label
ylabel2 = (f'<text transform="rotate(-90,18,{PAD_T2 + AH2/2:.0f})" '
           f'x="18" y="{PAD_T2 + AH2/2:.0f}" text-anchor="middle" '
           f'font-size="10" fill="#555">Avg log fold change (PFC enrichment)</text>')

svg2 = f"""<svg viewBox="0 0 {FW2} {FH2}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{FW2//2}" y="22" text-anchor="middle" font-size="13" font-weight="600" fill="#222">
    PFC Enrichment of Radial Glia Markers Across Developmental Stages
  </text>
  <text x="{FW2//2}" y="40" text-anchor="middle" font-size="10" fill="#666">
    Radial glia (RG) cell type · MAST log fold change (PFC vs. other cortical areas)
  </text>
  <text x="{FW2//2}" y="56" text-anchor="middle" font-size="10" fill="#444">
    HOPX enrichment increases mid-stage; CLU peaks early then declines — two distinct arealization dynamics
  </text>
  {axis2}{yticks2}{ylabel2}{bars2}{leg2}
</svg>"""

with open(os.path.join(OUT, "stage_dynamics.svg"), "w") as f:
    f.write(svg2)
print("Wrote stage_dynamics.svg")
