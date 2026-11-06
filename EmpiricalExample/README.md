# Empirical Example

This folder contains the materials for the Empirical Example section of

> Li, J., Patz, R. J., & van der Linden, W. J. OptimalCAT: Bayesian adaptive testing with optimal test assembly. *Behavior Research Methods*.

It holds everything needed to repeat the simulations reported there and to recompute Table 2, Figures 3–6 and the other results stated in the text: the engine configurations, the simulated examinees and their item responses, the item-parameter draws, the engine's outputs, and the scripts.

## The example

The example is a 40-item adaptive mathematics test: 35 operational (OP) items and 5 field-test (FT) items, drawn from a bank of 200 OP and 200 FT items. Every test must meet 231 constraints. The test was administered to 3,600 simulated examinees, 100 at each of 36 quantiles of N(0, 1), under four conditions:

| Condition | Test ID | Ability estimation | Item parameters | Exposure control |
|---|---|---|---|---|
| 1 | `brm1` | maximum likelihood (engine default) | fixed | none |
| 2 | `brm2` | maximum likelihood | fixed | 0.30 |
| 3 | `brm3` | EAP, N(0, 1) prior | fixed | 0.30 |
| 4 | `brm4cs` | fully Bayesian (MCMC) | 500 posterior draws per item | 0.30 |

There are two further test IDs:
- `brm2s` is a copy of Condition 2, used to run the simulator bundled with OptimalCAT (Figure 3).
- `brmlim` is Condition 3 without exposure control. It was used for the theoretical limit on the standard error (Figure 4) and for pre-flight checks, so that no probe entered a reported condition's exposure history.

Each examinee's response to every OP item in the bank was generated once, at the true ability and the bank's item parameters, and every condition reads that one file. The four conditions are therefore paired by examinee.

## Requirements

- **OptimalCAT 2.4**, the free version for research use (www.optimal-cat.com). It needs FICO Xpress, and the free Xpress Community License is sufficient for these problems.
- **R** with the `httr` and `jsonlite` packages. The R client drives the engine through its API (`POST /cat`).
- **Python 3** with `numpy`, `scipy` and `matplotlib`, for Table 2, Figures 4–6 and the item-parameter draws.

The reported runs used OptimalCAT 2.4 with Xpress on a Windows laptop. The R client ran on a separate computer on the same local network.

## Contents

| Path | What it is |
|---|---|
| `config/` | Engine configuration workbooks for the six tests, and the item-parameter draws file read by `brm4cs.xlsx`. |
| `input/responses_seed20261030_n3600_gpcmplus.csv` | True θ and the response of each of the 3,600 examinees to every item. |
| `input/math_bayesian_normprior.csv` | Template for the draws file: the N(0, 1) prior row (500 draws) and the item rows to be replaced. |
| `input/math_bayesian_normprior_calibN1000_steps.csv` | The item-parameter draws used in Condition 4 (same file as in `config/`). |
| `input/brm_constraint_sets.csv` | The 231 constraints in a form the client uses to write `constraints.csv`. |
| `input/s_brm2s.csv`, `input/simulation_brm2s.txt` | Simulator input for Figure 3: the same examinees and responses, in the same order. |
| `results/` | One folder per run (below). |
| `figures/` | Figures 3–6 as produced (PNG and PDF). |
| `brm_driver.R`, `ocat_common.R` | The R client that administers the simulated tests through the engine's API. |
| `make_brm_calib_draws.py` | Generates the Condition 4 item-parameter draws. |
| `brm_table2.py` | Recomputes every entry of Table 2 from the run folders, and the related statements in the text. |
| `exact_posterior_4cs.py` | Computes the exact posterior for every Condition 4 examinee on the items they were given, as a reference for the engine's sampled posterior. |
| `make_brm_figures.py` | Draws Figures 4–6 from the run folders. |

### Results, by table and figure

| Folder in `results/` | Used for |
|---|---|
| `2026-10-30_012151_brm1` | Table 2, Condition 1; Figure 4 |
| `2026-10-30_014920_brm2` | Table 2, Condition 2; Figure 4 |
| `2026-10-30_021646_brm3` | Table 2, Condition 3 |
| `2026-10-30_162500_brm4cs` | Table 2, Condition 4 |
| `2026-10-30_011719_brmlim_limit` | Theoretical limit on the SE (Figure 4) |
| `brm2s_20261030_192335` | Figure 3: output of the bundled simulator, including `statistics.png` |
| `2026-10-30_191311_brm4cs_students` | Figures 5 and 6: posterior draws after each item for Students A (θ = 1.5) and B (θ = −1.5) |
| `2026-11-01_122222_brm4cs_updatecheck` | Condition 4: the engine's posterior compared with the exact posterior after a single update and as updates accumulate |
| `2026-11-02_091255_brm4cs_bigS` | Condition 4: every test rescored with 5,000 rather than 500 posterior draws retained at each update |

Each run folder written by the R client contains the following files.

Results:
- `scores.csv`, `constraints.csv` and `exposure.csv`, in the same layout the simulator uses;
- `items_long.csv`: the item administered at each position, the response, and the interim estimate;
- `timing.csv`: the round-trip time of every API call;
- `positions_check.csv`.

Provenance:
- `run_meta.json`;
- `sessionInfo.txt`;
- a checksum of the response file;
- the engine's exposure-control state before and after the run.

## Reproducing the results

**From the outputs provided** (no engine needed):

```
python3 brm_table2.py results                 # every entry of Table 2
python3 make_brm_figures.py results figures   # Figures 4, 5 and 6
python3 exact_posterior_4cs.py                # exact posterior for Condition 4
```

The comparison behind the statement that Condition 4's slightly lower coverage comes from Monte Carlo error accumulated over the sequential updates is in two results folders:
- `2026-11-01_122222_brm4cs_updatecheck`: a single update agrees with the exact posterior, and the error grows as updates accumulate;
- `2026-11-02_091255_brm4cs_bigS`: retaining 5,000 draws at each update brings RMSE, mean SD and coverage close to the exact posterior. See `rescore_bigS_summary.csv`, to be read against the output of `exact_posterior_4cs.py`.

**Repeating the simulations:**

1. **Load the tests.** Copy the six workbooks in `config/` and the draws file into the engine's `data` folder. Build and check them with the engine's `buildjson` and `checkjson` utilities, then restart the engine. `GET /` on the engine should list the six test IDs.
2. **Run the conditions from R,** with this folder as the working directory:
   ```r
   options(ocat.host = "localhost")   # or the address of the computer running the engine
   source("brm_driver.R")
   brm_selftest()          # offline checks
   brm_preflight()         # engine reachable, tests loaded, scoring checked (uses brmlim)
   brm_probe_limit()       # theoretical limit on the SE
   brm_run("brm1"); brm_run("brm2"); brm_run("brm3"); brm_run("brm4cs")
   brm_capture_students()  # Students A and B, Condition 4
   brm_update_check()      # Condition 4: single update vs accumulated updates
   ```
   For the 5,000-draw rescoring, set `bayesian.post_burn_in = 50000` in the engine's `application.properties` (thinning 10, so 5,000 draws are retained at each update), restart the engine, run `brm_rescore_bigS()`, then restore the default and restart.
   Each call writes a new time-stamped folder in `results/`. Run each condition once after an engine restart: exposure-control rates accumulate per test from the time the engine starts.
3. **Figure 3:** on the computer running the engine, from the simulator folder with `s_brm2s.csv` beside it, run `java -jar ocat-simulator-2.4.jar simulation_brm2s.txt`.

**What to expect from a repeated run.** It should reproduce the reported statistics within simulation error, but not necessarily item for item, for two reasons:
- each shadow test is solved to a relative MIP gap of 0.005 within a time limit, so near-tied solutions can differ between runs;
- Condition 4's MCMC updates are stochastic.

The response file itself is fixed. `brm_responses()` in `brm_driver.R` reads it and checks it against the design. If the file is absent, `brm_make_responses()` regenerates it from its seed.

## Item-parameter draws for Condition 4

`make_brm_calib_draws.py --steps` generates the draws to represent the posterior that a calibration of 1,000 examinees from N(0, 1) would produce.

- **How each item's draws are generated.** For each OP item, the draws come from a normal approximation centred on the bank's parameter values. Its covariance is the inverse of the Fisher information of 1,000 examinees.
- **GPC items.** Their location and step parameters are determined only up to a common shift. So the approximation is for the discrimination (on the log scale) and the four category thresholds, which are identified. Each draw is stored with B equal to the mean threshold and steps that sum to zero, which gives the same response probabilities.
- **The prior row.** Row 1 holds the 500 draws that represent the N(0, 1) ability prior.

## Configuration workbooks

The six workbooks share the item bank, the 231 constraints and the test design of `Application Example/math_bayesian.xlsx` in this repository. They differ from it, and from each other, only in these entries of the `Test Settings` sheet:

| Entry | brm1 | brm2, brm2s | brm3 | brm4cs | brmlim |
|---|---|---|---|---|---|
| Test Interim Theta Estimation | (default) | (default) | EAP NORMAL 0 1 | (not used) | EAP NORMAL 0 1 |
| Test Final Theta Estimation | (default) | (default) | EAP NORMAL 0 1 | (not used) | EAP NORMAL 0 1 |
| Item Exposure Control Rate | (none) | 0.3 | 0.3 | 0.3 | (none) |
| Bayesian Draws | (none) | (none) | (none) | math_bayesian_normprior_calibN1000_steps.csv | (none) |
| Pre-Generate Shadow Test Forms at Theta | −2, −1, 0, 1, 2 | −2, −1, 0, 1, 2 | −2, −1, 0, 1, 2 | −2, −1, 0, 1, 2 | −2, −1, 0, 1, 2 |

The engine's defaults for a blank estimation entry are EAP under N(0, 1) during the test until the responses include both correct and incorrect answers, then maximum likelihood, with a maximum-likelihood final score.

## Licence

The scripts are released under the MIT License (`LICENSE`, at the top of this repository). The data and other materials are released under CC BY 4.0.
