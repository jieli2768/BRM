"""
make_brm_calib_draws.py -- item-parameter draws for Condition 4.

For each OP item, 500 draws from the large-sample posterior that a calibration of N
examinees from N(0, 1) would give: a normal distribution centred on the bank's parameter
values, with covariance equal to the inverse of N times the Fisher information averaged
over N(0, 1). 2PL items are parameterized as (log a, b). GPC items are parameterized as
(log a, delta_1..delta_m), the category thresholds delta_k = b - d_k, which are
identified; each draw is written with B equal to the mean threshold and steps summing
to zero, which gives the same response probabilities.

Use:  python3 make_brm_calib_draws.py --steps [N ...]      (default N = 1000)
Reads input/math_bayesian_normprior.csv (the N(0, 1) prior row is kept) and the item
bank from results/2026-10-30_162500_brm4cs/testdef_brm4cs.json. Writes
input/math_bayesian_normprior_calibN<N>_steps.csv and copies it to config/.
"""
import csv, math, random, shutil, statistics as st, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "input" / "math_bayesian_normprior.csv"
TESTDEF = HERE / "results" / "2026-10-30_162500_brm4cs" / "testdef_brm4cs.json"
SEED = 20261030
Q = [-6 + 12 * i / 400 for i in range(401)]
W = [math.exp(-x * x / 2) for x in Q]; _s = sum(W); W = [w / _s for w in W]

def info_loga_b(it):
    a, b, sc = it["a"], it["b"], it["sc"]
    steps = [0.0] + list(it["d"]) if it["model"] == "GPCM" else [0.0, 0.0]
    Iaa = Iab = Ibb = 0.0
    for th, w in zip(Q, W):
        T, X = [0.0], [0]
        for k in range(1, len(steps)):
            T.append(T[-1] + (th - b + steps[k])); X.append(k)
        z = [sc * a * t for t in T]; m = max(z); e = [math.exp(v - m) for v in z]; S = sum(e)
        p = [v / S for v in e]
        ET = sum(pi * t for pi, t in zip(p, T)); EX = sum(pi * x for pi, x in zip(p, X))
        VT = sum(pi * (t - ET) ** 2 for pi, t in zip(p, T))
        VX = sum(pi * (x - EX) ** 2 for pi, x in zip(p, X))
        C = sum(pi * (t - ET) * (x - EX) for pi, t, x in zip(p, T, X))

        Iaa += w * (sc * a) ** 2 * VT
        Ibb += w * (sc * a) ** 2 * VX
        Iab += w * -(sc * a) ** 2 * C
    return Iaa, Iab, Ibb

def info_delta(it):
    a, sc = it["a"], it["sc"]
    delta = [it["b"] - dk for dk in it["d"]] if it["model"] == "GPCM" else [it["b"]]
    mm = len(delta); n = mm + 1
    I = [[0.0] * n for _ in range(n)]
    for th, w in zip(Q, W):
        T = [0.0]
        for k in range(mm):
            T.append(T[-1] + th - delta[k])
        z = [sc * a * t for t in T]; mx = max(z); e = [math.exp(v - mx) for v in z]; S = sum(e)
        p = [v / S for v in e]

        F = [[sc * a * T[x]] + [-sc * a * (1.0 if x >= k + 1 else 0.0) for k in range(mm)]
             for x in range(mm + 1)]
        mu = [sum(p[x] * F[x][j] for x in range(mm + 1)) for j in range(n)]
        for i in range(n):
            for j in range(n):
                I[i][j] += w * sum(p[x] * (F[x][i] - mu[i]) * (F[x][j] - mu[j]) for x in range(mm + 1))
    return I, delta

def _inv(M):
    n = len(M); A = [r[:] + [float(i == j) for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(A[r][c])); A[c], A[piv] = A[piv], A[c]
        pv = A[c][c]; A[c] = [v / pv for v in A[c]]
        for r in range(n):
            if r != c:
                f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [r[n:] for r in A]

def _chol(V):
    n = len(V); L = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = V[i][j] - sum(L[i][k] * L[j][k] for k in range(j))
            L[i][j] = math.sqrt(s) if i == j else s / L[j][j]
    return L

def main():
    import json
    steps = "--steps" in sys.argv
    Ns = [int(x) for x in sys.argv[1:] if x != "--steps"] or [1000]
    bank = {it["itm_id"]: it for it in json.load(open(TESTDEF))["item_bank"] if it["status"] == "OP"}
    src = list(csv.reader(open(SRC, newline="")))
    S = len(src[0]) - 2
    for N in Ns:
        rng = random.Random(SEED + N)
        rows, sdb, sdla, rho = [], [], [], []
        draws = {}
        for iid in sorted(bank):
            it = bank[iid]
            I, delta = info_delta(it)
            V = [[v / N for v in r] for r in _inv(I)]
            L = _chol(V); n = len(V); mm = n - 1
            mdel = sum(delta) / mm
            la, bb, dl = [], [], []
            for _ in range(S):
                zz = [rng.gauss(0, 1) for _ in range(n)]
                x = [sum(L[i][k] * zz[k] for k in range(i + 1)) for i in range(n)]
                la.append(math.log(it["a"]) + x[0])
                bb.append(it["b"] + sum(x[1:]) / mm)
                dl.append([delta[k] + x[1 + k] for k in range(mm)])
            vla = V[0][0]
            vb = sum(V[i][j] for i in range(1, n) for j in range(1, n)) / mm ** 2
            cov = sum(V[0][j] for j in range(1, n)) / mm
            draws[iid] = (la, bb, dl, it["model"] == "GPCM")
            sdb.append(math.sqrt(vb)); sdla.append(math.sqrt(vla)); rho.append(cov / math.sqrt(vla * vb))
        for r in src:
            if r[0] in ("entry", "STUDENT"):
                rows.append(r); continue
            la, bb, dl, gpcm = draws[r[0]]
            if steps and gpcm:
                bs = [sum(dd) / len(dd) for dd in dl]
                if r[1] == "A":
                    rows.append([r[0], "A"] + ["%.6f" % math.exp(v) for v in la])
                else:
                    rows.append([r[0], "B"] + ["%.6f" % v for v in bs])
                    for k in range(len(dl[0])):
                        rows.append([r[0], "D%d" % (k + 1)] + ["%.6f" % (b_ - dd[k]) for b_, dd in zip(bs, dl)])
                continue
            vals = [math.exp(v) for v in la] if r[1] == "A" else bb
            rows.append([r[0], r[1]] + ["%.6f" % v for v in vals])
        out = HERE / "input" / ("math_bayesian_normprior_calibN%d%s.csv" % (N, "_steps" if steps else ""))
        with open(out, "w", newline="") as f:
            csv.writer(f).writerows(rows)
        chk = list(csv.reader(open(out, newline="")))
        assert chk[0] == src[0] and chk[1] == src[1]
        assert [r[:2] for r in chk if r[1] in ("A", "B", "THETA") or r[0] == "entry"] == [r[:2] for r in src]
        if steps:
            byk = {(r[0], r[1]): [float(v) for v in r[2:]] for r in chk[2:]}
            for iid, (la, bb, dl, gpcm) in draws.items():
                if not gpcm: continue
                for k in range(len(dl[0])):
                    assert max(abs(b_ - d_ - dd[k]) for b_, d_, dd in zip(byk[(iid, "B")], byk[(iid, "D%d" % (k + 1))], dl)) < 2e-6
                assert max(abs(sum(byk[(iid, "D%d" % (k + 1))][s_] for k in range(len(dl[0])))) for s_ in range(S)) < 1e-5
            print("  --steps: %d D rows written; b - d_k reproduces every drawn threshold; steps sum to zero in every draw"
                  % sum(1 for r in chk if r[1].startswith("D")))
        shutil.copy(out, HERE / "config" / out.name)
        q = lambda v, p: sorted(v)[int(p * (len(v) - 1))]
        bvals = {i: bank[i]["b"] for i in bank}
        print("%s: %d item rows, prior row unchanged, copied to config/" % (out.name, len(chk) - 2))
        print("  theoretical sd(b): median %.3f IQR %.3f-%.3f range %.3f-%.3f"
              % (st.median(sdb), q(sdb, .25), q(sdb, .75), min(sdb), max(sdb)))
        print("  theoretical sd(log a): median %.3f IQR %.3f-%.3f range %.3f-%.3f"
              % (st.median(sdla), q(sdla, .25), q(sdla, .75), min(sdla), max(sdla)))
        print("  corr(log a, b): median %.2f range %.2f-%.2f" % (st.median(rho), min(rho), max(rho)))
        ids = sorted(bank)
        for lo, hi in [(-9, -0.5), (-0.5, 0.5), (0.5, 1.5), (1.5, 2.5), (2.5, 9)]:
            s = [sdb[j] for j, i in enumerate(ids) if lo <= bvals[i] < hi]
            if s: print("  b in [%g,%g): n %3d, median sd(b) %.3f" % (lo, hi, len(s), st.median(s)))
        emp = [st.pstdev(draws[i][1]) / sdb[j] for j, i in enumerate(ids)]
        print("  empirical/theoretical sd(b) over 500 draws: %.3f-%.3f" % (min(emp), max(emp)))

if __name__ == "__main__":
    main()
