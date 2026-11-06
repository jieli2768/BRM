"""
make_brm_figures.py -- Figures 4, 5 and 6 of the Empirical Example.

Use:  python3 make_brm_figures.py [results_dir] [output_dir]     (defaults: . and figures)
  Figure 4: mean reported SE at each of the 36 true abilities, Conditions 1 and 2, with
            the theoretical limit.
  Figure 5: posterior density of theta for Students A and B after 0, 10, 20, 30, 40 items.
  Figure 6: posterior mean +- 1 SD by item position, Students A and B.
Requires numpy, scipy and matplotlib.
"""
import csv, math, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

RES = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("figures")
OUT.mkdir(parents=True, exist_ok=True)
RUN1 = RES / "2026-10-30_012151_brm1" / "scores.csv"
RUN2 = RES / "2026-10-30_014920_brm2" / "scores.csv"
LIMIT = RES / "2026-10-30_011719_brmlim_limit" / "theoretical_limit.csv"
STU = RES / "2026-10-30_191311_brm4cs_students"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": "#dddddd",
                     "grid.linewidth": 0.6, "legend.frameon": False})
INK, MID, LIGHT = "#1f1f1f", "#5a5a5a", "#9a9a9a"

def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / f"{name}.png")

def mean_se_by_theta(path):
    acc = {}
    for r in csv.DictReader(open(path)):
        if r["domain_name"] != "Overall":
            continue
        acc.setdefault(float(r["true_theta"]), []).append(float(r["se"]))
    th = sorted(acc)
    return np.array(th), np.array([np.mean(acc[t]) for t in th]), np.array([np.median(acc[t]) for t in th])

th1, se1, md1 = mean_se_by_theta(RUN1)
th2, se2, md2 = mean_se_by_theta(RUN2)
lim = [(float(r["theta"]), float(r["limit_sem"])) for r in csv.DictReader(open(LIMIT))]
lt, ls = np.array([x for x, _ in lim]), np.array([y for _, y in lim])
keep = (lt >= -2.5) & (lt <= 2.5)
fig, ax = plt.subplots(figsize=(6.5, 4.0))
ax.plot(lt[keep], ls[keep], color=INK, lw=2.2, label="Theoretical limit")
ax.plot(th1, se1, color=MID, lw=1.4, ls="--", marker="o", ms=3, label="Condition 1: no exposure control")
ax.plot(th2, se2, color=INK, lw=1.4, ls=":", marker="s", ms=3, label="Condition 2: exposure-control rate 0.30")
ax.set_xlabel(r"True ability $\theta$")
ax.set_ylabel("Mean reported SE")
ax.set_xlim(-2.5, 2.5)
ax.set_ylim(0.2, 0.9)
ax.legend(loc="upper right")
save(fig, "fig_sem_limit")

print("cell -2.20: mean SE cond1 %.3f (median %.3f), cond2 %.3f (median %.3f)" % (se1[0], md1[0], se2[0], md2[0]))

draws = {}
for r in csv.DictReader(open(STU / "students_draws.csv")):
    draws[(r["student"], int(r["position"]))] = np.array([float(r["v%d" % i]) for i in range(1, 501)])
traj = {}
for r in csv.DictReader(open(STU / "students_trajectory.csv")):
    traj.setdefault(r["student"], []).append((int(r["position"]), float(r["post_mean"]), float(r["post_sd"]),
                                              float(r["true_theta"])))
LAB = {"A": "Student A (true $\\theta$ = 1.5)", "B": "Student B (true $\\theta$ = −1.5)"}

fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.4), sharey=True)
grid = np.linspace(-4, 4, 400)
styles = {0: (LIGHT, ":", 1.2), 10: (LIGHT, "-", 1.2), 20: (MID, "--", 1.3), 30: (MID, "-", 1.3), 40: (INK, "-", 2.0)}
for ax, s in zip(axes, ("A", "B")):
    for pos, (c, lst, lw) in styles.items():
        v = draws[(s, pos)]
        ax.plot(grid, gaussian_kde(v)(grid), color=c, ls=lst, lw=lw, label=f"after {pos} items" if pos else "prior (0 items)")
    t = traj[s][0][3]
    ax.axvline(t, color=INK, lw=0.8, ls=(0, (2, 2)))
    ax.set_title(LAB[s], fontsize=10)
    ax.set_xlabel(r"$\theta$")
axes[0].set_ylabel("Posterior density")
axes[1].legend(loc="upper right", fontsize=8)
save(fig, "fig_post_density")

fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.4), sharey=True)
for ax, s in zip(axes, ("A", "B")):
    p = np.array([x[0] for x in traj[s]]); m = np.array([x[1] for x in traj[s]]); sd = np.array([x[2] for x in traj[s]])
    ax.plot(p, m, color=INK, lw=1.6, label="Posterior mean")
    ax.plot(p, m + sd, color=MID, lw=1.0, ls=":", label="Mean ± 1 SD")
    ax.plot(p, m - sd, color=MID, lw=1.0, ls=":")
    ax.axhline(traj[s][0][3], color=INK, lw=0.8, ls=(0, (4, 3)), label="True $\\theta$")
    ax.set_title(LAB[s], fontsize=10)
    ax.set_xlabel("Item position")
    ax.set_xlim(0, 40)
axes[0].set_ylabel(r"$\theta$")
axes[0].legend(loc="lower right", fontsize=8)
save(fig, "fig_post_path")
for s in ("A", "B"):
    last = traj[s][-1]
    print("student %s: final mean %.3f, SD %.3f (true %.1f)" % (s, last[1], last[2], last[3]))
