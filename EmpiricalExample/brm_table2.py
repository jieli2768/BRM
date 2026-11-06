"""
brm_table2.py -- every entry of Table 2, from the four run folders.

Use:  python3 brm_table2.py [results_dir]      (default: results)
Also prints the related statements in the text: OP items whose exposure rate the 0.30
target could not hold, and final ML estimates at the scale boundary. The last column
is the simulator's run of Condition 2 (Figure 3), which records no timing or position
checks. OP items are I001-I200; FT items are excluded from the exposure statistics.
"""
import csv, sys, statistics
from pathlib import Path

RES = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results")
RUNS = [("1", "2026-10-30_012151_brm1"), ("2", "2026-10-30_014920_brm2"),
        ("3", "2026-10-30_021646_brm3"), ("4", "2026-10-30_162500_brm4cs"),
        ("2 (simulator)", "brm2s_20261030_192335")]

def rows(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))

def summarize(run):
    d = RES / run
    s = [r for r in rows(d / "scores.csv") if r["domain_name"] == "Overall"]
    th = [float(r["true_theta"]) for r in s]; est = [float(r["theta"]) for r in s]
    se = [float(r["se"]) for r in s]; err = [e - t for e, t in zip(est, th)]
    def bias_at(v):
        e = [x for x, t in zip(err, th) if round(t * 2) / 2 == v]
        return sum(e) / len(e), len(e)
    ex = [r for r in rows(d / "exposure.csv") if r["type"] == "item" and int(r["id"][1:]) <= 200]
    rate = [float(r["exposure_rate"]) for r in ex]
    con = rows(d / "constraints.csv")
    pos = rows(d / "positions_check.csv") if (d / "positions_check.csv").exists() else None
    t = [float(r["ms"]) for r in rows(d / "timing.csv") if int(r["call"]) > 1] if (d / "timing.csv").exists() else None
    return {
        "N": len(s),
        "Bias, pooled": sum(err) / len(err),
        "Bias at theta ~ -1.5": bias_at(-1.5),
        "Bias at theta ~ +1.5": bias_at(1.5),
        "RMSE": (sum(e * e for e in err) / len(err)) ** 0.5,
        "Mean SE": sum(se) / len(se),
        "Coverage, +-1.96 SE": sum(abs(e) <= 1.96 * v for e, v in zip(err, se)) / len(err),
        "Exposure rate, max": max(rate),
        "Exposure rate, mean": sum(rate) / len(rate),
        "Exposure rate, SD": statistics.pstdev(rate),
        "OP pool administered": sum(r > 0 for r in rate) / len(rate),

        "OP items above 0.31": sum(r > 0.31 for r in rate),
        "Constraints met (n, all, positions)": (len(con), all(float(r["on_target"]) == 1 for r in con),
                            "" if pos is None else all(float(r["all_members_in_position"]) == 1 for r in pos)),
        "Median response time (ms)": "" if t is None else statistics.median(t),
        "Final estimates at -5 (raw score)": [(r["student"], r["raw_score"]) for r in s if float(r["theta"]) <= -4.99],
    }

out = {c: summarize(r) for c, r in RUNS}
for k in out["1"]:
    print("%-34s" % k, " | ".join(
        ("%.3f" % v) if isinstance(v, float) else
        ("%.3f (n=%d)" % v) if isinstance(v, tuple) and isinstance(v[0], float) else str(v)
        for v in (out[c][k] for c, _ in RUNS)))
