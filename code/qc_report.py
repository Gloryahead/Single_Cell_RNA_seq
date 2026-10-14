#!/usr/bin/env python3
"""
qc_report.py — Cell Ranger QC report for all samples.

Reads metrics_summary.csv from each sample's Cell Ranger output, evaluates
standard 10x PBMC thresholds, and writes:
  results/qc_report.html  — color-coded pass/fail table
  results/qc_report.csv   — machine-readable version

Run (reads only tiny CSV files — login node is fine):
  python code/qc_report.py
"""
import os, csv, sys, yaml
from pathlib import Path
from datetime import datetime

PROJECT = os.environ.get("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
DATA    = os.environ.get("SCRNA_DATA",    "/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq")
os.chdir(PROJECT)

cfg      = yaml.safe_load(open("config.yaml"))
out_dir  = Path(DATA) / "results"
samples  = [s["id"] for s in cfg["samples"]]
proj     = cfg["project"].get("name", "scRNA-seq")

# ── Thresholds ──────────────────────────────────────────────────────────────────
# "good" and "acceptable" are lower bounds; anything below "acceptable" is "poor".
# Estimated Cells uses an absolute range rather than a fold-change vs. a target.
METRICS = [
    {
        "col":   "Estimated Number of Cells",
        "label": "Estimated Cells",
        "good":        lambda x: 5_000 <= x <= 12_000,
        "acceptable":  lambda x: 2_000 <= x <= 15_000,
        "fmt":         lambda x: f"{x:,.0f}",
        "thresh_text": "Good 5 k–12 k | Acceptable 2 k–15 k | Poor outside",
    },
    {
        "col":   "Mean Reads per Cell",
        "label": "Mean Reads / Cell",
        "good":        lambda x: x >= 20_000,
        "acceptable":  lambda x: x >= 10_000,
        "fmt":         lambda x: f"{x:,.0f}",
        "thresh_text": "Good ≥20 k | Acceptable ≥10 k | Poor <10 k",
    },
    {
        "col":   "Median Genes per Cell",
        "label": "Median Genes / Cell",
        "good":        lambda x: x >= 1_500,
        "acceptable":  lambda x: x >= 1_000,
        "fmt":         lambda x: f"{x:,.0f}",
        "thresh_text": "Good ≥1,500 | Acceptable ≥1,000 | Poor <1,000",
    },
    {
        "col":   "Sequencing Saturation",
        "label": "Seq. Saturation",
        "good":        lambda x: x >= 60,
        "acceptable":  lambda x: x >= 40,
        "fmt":         lambda x: f"{x:.1f}%",
        "thresh_text": "Good ≥60% | Acceptable ≥40% | Poor <40%",
    },
    {
        "col":   "Fraction Reads in Cells",
        "label": "Reads in Cells",
        "good":        lambda x: x >= 70,
        "acceptable":  lambda x: x >= 50,
        "fmt":         lambda x: f"{x:.1f}%",
        "thresh_text": "Good ≥70% | Acceptable ≥50% | Poor <50%",
    },
    {
        "col":   "Valid Barcodes",
        "label": "Valid Barcodes",
        "good":        lambda x: x >= 90,
        "acceptable":  lambda x: x >= 80,
        "fmt":         lambda x: f"{x:.1f}%",
        "thresh_text": "Good ≥90% | Acceptable ≥80% | Poor <80%",
    },
    {
        "col":   "Q30 Bases in Barcode",
        "label": "Q30 Barcode",
        "good":        lambda x: x >= 90,
        "acceptable":  lambda x: x >= 80,
        "fmt":         lambda x: f"{x:.1f}%",
        "thresh_text": "Good ≥90% | Acceptable ≥80% | Poor <80%",
    },
    {
        "col":   "Q30 Bases in UMI",
        "label": "Q30 UMI",
        "good":        lambda x: x >= 90,
        "acceptable":  lambda x: x >= 80,
        "fmt":         lambda x: f"{x:.1f}%",
        "thresh_text": "Good ≥90% | Acceptable ≥80% | Poor <80%",
    },
    {
        "col":   "Q30 Bases in RNA Read",
        "label": "Q30 RNA Read",
        "good":        lambda x: x >= 75,
        "acceptable":  lambda x: x >= 65,
        "fmt":         lambda x: f"{x:.1f}%",
        "thresh_text": "Good ≥75% | Acceptable ≥65% | Poor <65%",
    },
    {
        "col":   "Reads Mapped Confidently to Transcriptome",
        "label": "Mapped to Transcriptome",
        "good":        lambda x: x >= 70,
        "acceptable":  lambda x: x >= 50,
        "fmt":         lambda x: f"{x:.1f}%",
        "thresh_text": "Good ≥70% | Acceptable ≥50% | Poor <50%",
    },
]

def parse_val(raw):
    return float(raw.strip().strip('"').replace(",", "").replace("%", ""))

def rate(val, m):
    if m["good"](val):       return "good"
    if m["acceptable"](val): return "acceptable"
    return "poor"

ICON = {"good": "✅", "acceptable": "⚠️", "poor": "❌", "missing": "—"}
BG   = {"good": "#d1e7dd", "acceptable": "#fff3cd", "poor": "#f8d7da", "missing": "#f8f9fa"}
FG   = {"good": "#0a3622", "acceptable": "#664d03", "poor": "#58151c", "missing": "#6c757d"}

# ── Read metrics_summary.csv for each sample ────────────────────────────────────
records = []
for sid in samples:
    csv_path = out_dir / "counts" / sid / "outs" / "metrics_summary.csv"
    if not csv_path.exists():
        print(f"WARNING: {csv_path} not found — skipping {sid}", file=sys.stderr)
        continue
    with open(csv_path, newline="") as f:
        row = next(csv.DictReader(f))
    rec = {"sample": sid, "metrics": {}}
    for m in METRICS:
        raw = row.get(m["col"], "").strip()
        if not raw:
            rec["metrics"][m["col"]] = {"val": None, "fmt": "—", "rating": "missing"}
        else:
            val = parse_val(raw)
            rec["metrics"][m["col"]] = {"val": val, "fmt": m["fmt"](val), "rating": rate(val, m)}
    ratings = [v["rating"] for v in rec["metrics"].values()]
    rec["overall"] = (
        "poor"       if "poor"       in ratings else
        "acceptable" if "acceptable" in ratings else
        "good"
    )
    records.append(rec)

if not records:
    sys.exit("No metrics_summary.csv files found. Run Cell Ranger first.")

# ── CSV ────────────────────────────────────────────────────────────────────────
csv_out = out_dir / "qc_report.csv"
with open(csv_out, "w", newline="") as f:
    cols = ["sample"] + [m["col"] for m in METRICS] + ["overall"]
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    for rec in records:
        row = {"sample": rec["sample"], "overall": rec["overall"]}
        for m in METRICS:
            e = rec["metrics"][m["col"]]
            row[m["col"]] = e["fmt"]
        w.writerow(row)
print(f"CSV  -> {csv_out}")

# ── HTML ───────────────────────────────────────────────────────────────────────
n_good = sum(1 for r in records if r["overall"] == "good")
n_acc  = sum(1 for r in records if r["overall"] == "acceptable")
n_poor = sum(1 for r in records if r["overall"] == "poor")

# Rotate sample header IDs slightly for readability in narrow cells
sample_headers = "".join(
    f'<th class="sample-header"><div><span>{r["sample"]}</span></div></th>'
    for r in records
)

metric_rows = []
for m in METRICS:
    cells = []
    for rec in records:
        e   = rec["metrics"][m["col"]]
        bg  = BG[e["rating"]]
        fg  = FG[e["rating"]]
        ico = ICON[e["rating"]]
        cells.append(
            f'<td style="background:{bg};color:{fg};text-align:center;'
            f'white-space:nowrap;padding:6px 10px;">'
            f'{ico} {e["fmt"]}</td>'
        )
    metric_rows.append(
        f'<tr>'
        f'<td class="metric-name" title="{m["thresh_text"]}">{m["label"]}</td>'
        f'<td class="thresh-cell">{m["thresh_text"]}</td>'
        + "".join(cells) +
        f'</tr>'
    )

overall_cells = "".join(
    f'<td style="background:{BG[r["overall"]]};color:{FG[r["overall"]]};'
    f'text-align:center;font-weight:600;padding:6px 10px;">'
    f'{ICON[r["overall"]]} {r["overall"].capitalize()}</td>'
    for r in records
)

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{proj} — Cell Ranger QC</title>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
          font-size: 14px; background: #f5f5f5; color: #212529; padding: 24px; }}
  h1   {{ font-size: 1.4rem; font-weight: 600; margin-bottom: 4px; }}
  .meta {{ color: #6c757d; font-size: 0.85rem; margin-bottom: 20px; }}
  .summary {{ display: flex; gap: 12px; margin-bottom: 20px; flex-wrap: wrap; }}
  .badge {{ padding: 6px 14px; border-radius: 6px; font-weight: 600; font-size: 0.85rem; }}
  .badge-good {{ background: #d1e7dd; color: #0a3622; }}
  .badge-acc  {{ background: #fff3cd; color: #664d03; }}
  .badge-poor {{ background: #f8d7da; color: #58151c; }}
  .table-wrap {{ overflow-x: auto; border-radius: 8px;
                 box-shadow: 0 1px 4px rgba(0,0,0,.12); }}
  table {{ border-collapse: collapse; background: #fff; min-width: 100%; }}
  th, td {{ border: 1px solid #dee2e6; font-size: 13px; }}
  thead th {{ background: #343a40; color: #fff; padding: 8px 10px;
              font-weight: 600; white-space: nowrap; }}
  .metric-name {{ padding: 7px 12px; font-weight: 500; white-space: nowrap;
                  cursor: help; text-decoration: underline dotted #aaa; }}
  .thresh-cell {{ padding: 7px 12px; color: #6c757d; font-size: 11px;
                  white-space: nowrap; min-width: 240px; }}
  .sample-header div {{ writing-mode: vertical-rl; transform: rotate(180deg);
                        padding: 8px 4px; min-height: 90px;
                        display: flex; align-items: center; }}
  tr:hover td {{ filter: brightness(0.96); }}
  .overall-row td {{ font-weight: 600; border-top: 2px solid #343a40; }}
  .legend {{ margin-top: 16px; display: flex; gap: 16px; font-size: 12px;
             color: #495057; flex-wrap: wrap; }}
  .legend span {{ display: flex; align-items: center; gap: 4px; }}
</style>
</head>
<body>
<h1>{proj} — Cell Ranger QC Report</h1>
<p class="meta">Generated {datetime.now().strftime("%Y-%m-%d %H:%M")} &nbsp;|&nbsp;
{len(records)} samples &nbsp;|&nbsp; 10 metrics</p>

<div class="summary">
  <div class="badge badge-good">✅ Good: {n_good} sample{"s" if n_good != 1 else ""}</div>
  <div class="badge badge-acc">⚠️ Acceptable: {n_acc} sample{"s" if n_acc != 1 else ""}</div>
  <div class="badge badge-poor">❌ Poor: {n_poor} sample{"s" if n_poor != 1 else ""}</div>
</div>

<div class="table-wrap">
<table>
  <thead>
    <tr>
      <th>Metric</th>
      <th>Thresholds</th>
      {sample_headers}
    </tr>
  </thead>
  <tbody>
    {"".join(metric_rows)}
    <tr class="overall-row">
      <td class="metric-name">Overall</td>
      <td class="thresh-cell">All metrics good / any acceptable / any poor</td>
      {overall_cells}
    </tr>
  </tbody>
</table>
</div>

<div class="legend">
  <span>✅ Good</span>
  <span>⚠️ Acceptable</span>
  <span>❌ Poor — investigate before proceeding</span>
  <span style="color:#aaa">Hover metric names to see thresholds</span>
</div>
</body>
</html>
"""

html_out = out_dir / "qc_report.html"
html_out.write_text(html, encoding="utf-8")
print(f"HTML -> {html_out}")
