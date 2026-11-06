# brm_driver.R -- R client for the Empirical Example.
#
# Administers the simulated tests through the OptimalCAT API (POST /cat), one item per
# call, reading every examinee's responses from input/responses_seed20261030_n3600_gpcmplus.csv,
# and writes one time-stamped folder per run in results/.
#
# Use: with this folder as the working directory, set options(ocat.host = ...) to the
# computer running the engine and source("brm_driver.R"). The steps are listed in the
# RUN block at the end. Requires the R packages httr and jsonlite, and ocat_common.R.


# Folders: this script runs from the folder that holds config/, input/ and results/.
CELL     <- "BRM"
CELL_DIR <- normalizePath(getwd(), winslash = "/")
ROOT     <- CELL_DIR
options(ocat.cell_dir = CELL_DIR)
source(file.path(CELL_DIR, "ocat_common.R"))
DRIVER_PATH <- file.path(CELL_DIR, "brm_driver.R")

# Design: 36 quantiles of N(0, 1) at probabilities (i - 0.5)/36, 100 examinees at each;
# the six tests loaded in the engine; the four conditions in order.
SEED       <- 20261030L
N_PER_Q    <- 100L
QUANTILES  <- round(qnorm(((1:36) - 0.5) / 36), 4)
CONDITIONS <- c("brm1", "brm2", "brm3", "brm4cs")
PROBE_ID   <- "brmlim"
EXPECTED_TESTS <- c("brm1", "brm2", "brm3", "brm4cs", "brm2s", "brmlim")

# GPC response function as in the OptimalCAT manual:
# P(x) proportional to exp(sum_{k<=x} sc * a * (theta - b + d_k)).
GPCM_SIGN  <- +1

# The 231 constraints, used to summarize constraint conformance.
CONSTRAINT_SETS <- file.path(CELL_DIR, "input", "brm_constraint_sets.csv")

# True abilities of the 3,600 examinees, in the fixed order in which they are
# administered (the same order in every condition).
brm_thetas <- function() {
  set.seed(SEED)
  th <- rep(QUANTILES, each = N_PER_Q)
  th[sample.int(length(th))]
}

# GPC step parameters d_1..d_m for every item, from the test definition.
.steps_list <- function(ib) {
  d <- ib$d
  m <- as.integer(ib$max_points)
  out <- vector("list", nrow(ib))
  for (i in seq_len(nrow(ib))) {
    di <- if (is.list(d)) d[[i]] else if (is.matrix(d) || is.data.frame(d)) unlist(d[i, ]) else d[i]
    di <- suppressWarnings(as.numeric(unlist(di)))
    di <- di[is.finite(di)]
    if (!is.na(m[i]) && m[i] > 1) {
      if (length(di) == m[i] + 1) di <- di[-1]
      if (length(di) != m[i])
        stop("item ", ib$itm_id[i], ": max_points ", m[i], " but ", length(di), " step values")
    } else di <- numeric(0)
    out[[i]] <- di
  }
  out
}

# Item bank (OP items) as loaded in the engine, from GET /<testID>.
brm_bank <- function(testID = CONDITIONS[1], base_url = ocat_base_url()) {
  td <- ocat_get(paste0("/", testID), base_url)
  ib <- td$item_bank
  steps <- .steps_list(ib)
  op <- ib$status == "OP"
  col <- function(x, default) {
    if (is.null(x)) return(rep(default, sum(op)))
    v <- suppressWarnings(as.numeric(unlist(lapply(x, function(e) if (is.null(e) || !length(e)) NA else e[1]))))[op]
    v[is.na(v)] <- default; v
  }
  bank <- data.frame(itm_id = ib$itm_id[op], model = toupper(ib$model[op]),
                     max_points = as.integer(ib$max_points[op]),
                     sc = col(ib$sc, 1), a = col(ib$a, NA), b = col(ib$b, NA), c = col(ib$c, 0),
                     stringsAsFactors = FALSE)
  if (anyNA(bank$a) || anyNA(bank$b)) stop("operational items with missing a or b in GET /", testID)
  bank$steps <- steps[op]
  attr(bank, "all_ids") <- ib$itm_id
  attr(bank, "status")  <- setNames(ib$status, ib$itm_id)
  bank
}

# Category probabilities for one item at each theta.
item_probs <- function(item, theta, sign = GPCM_SIGN) {
  if (item$model %in% c("GPCM", "GPC")) {
    d <- item$steps
    z <- sapply(d, function(dk) item$sc * item$a * (theta - item$b + sign * dk))
    z <- matrix(z, nrow = length(theta))
    num <- cbind(0, t(apply(z, 1, cumsum)))
    if (length(d) == 1) num <- cbind(0, z)
    num <- num - apply(num, 1, max)
    p <- exp(num)
    p / rowSums(p)
  } else {
    p1 <- item$c + (1 - item$c) / (1 + exp(-item$sc * item$a * (theta - item$b)))
    cbind(1 - p1, p1)
  }
}

# Fisher information for one item at each theta.
item_info <- function(item, theta, sign = GPCM_SIGN) {
  p <- item_probs(item, theta, sign)
  if (item$model %in% c("GPCM", "GPC")) {
    x <- 0:(ncol(p) - 1)
    ex  <- p %*% x
    ex2 <- p %*% (x^2)
    as.numeric((item$sc * item$a)^2 * (ex2 - ex^2))
  } else {
    P <- p[, 2]; Q <- 1 - P
    as.numeric((item$sc * item$a)^2 * (Q / P) * ((P - item$c) / (1 - item$c))^2)
  }
}

# One item's parameters, in the form item_probs() and item_info() take.
bank_row <- function(bank, id) {
  i <- match(id, bank$itm_id)
  if (is.na(i)) return(NULL)
  c(as.list(bank[i, setdiff(names(bank), "steps")]), list(steps = bank$steps[[i]]))
}

# Generates every examinee's response to every item at the true theta.
brm_make_responses <- function(bank, thetas, seed = SEED + 1L, sign = GPCM_SIGN) {
  set.seed(seed)
  m <- matrix(NA_integer_, length(thetas), nrow(bank), dimnames = list(NULL, bank$itm_id))
  for (j in seq_len(nrow(bank))) {
    p <- item_probs(bank_row(bank, bank$itm_id[j]), thetas, sign)
    u <- runif(length(thetas))
    m[, j] <- rowSums(u > t(apply(p, 1, cumsum))[, -ncol(p), drop = FALSE])
  }
  data.frame(simulee = seq_along(thetas), true_theta = thetas, m, check.names = FALSE)
}

# Reads the response file (or generates it if absent) and checks it against the design.
brm_responses <- function(bank, thetas = brm_thetas()) {
  f <- file.path(CELL_DIR, "input",
                 sprintf("responses_seed%d_n%d_gpcm%s.csv", SEED, length(thetas),
                         if (GPCM_SIGN > 0) "plus" else "minus"))
  if (file.exists(f)) {
    resp <- read.csv(f, check.names = FALSE)
    if (nrow(resp) != length(thetas) || !isTRUE(all.equal(resp$true_theta, thetas)))
      stop("saved response matrix ", basename(f), " does not match the design; do not mix them")
    if (!all(bank$itm_id %in% names(resp))) stop("saved response matrix does not cover the bank")
    message("responses: read ", basename(f))
  } else {
    resp <- brm_make_responses(bank, thetas)
    write.csv(resp, f, row.names = FALSE)
    message("responses: generated and saved ", basename(f))
  }
  attr(resp, "path") <- f
  resp
}

# One POST /cat call.
brm_post <- function(testID, administered, sessionData = NULL, start_theta = NULL,
                     base_url = ocat_base_url()) {
  payload <- list(testID = testID, administeredItems = administered, sessionData = sessionData)
  if (!is.null(start_theta))
    payload$testStartTheta <- list(distribution = "FIXED", param1 = start_theta, param2 = 0)
  body <- toJSON(payload, auto_unbox = TRUE, digits = NA)
  t0 <- proc.time()[["elapsed"]]
  r <- POST(paste0(base_url, "/cat"), body = body, encode = "raw",
            add_headers(accept = "application/json"), content_type("application/json"))
  ms <- 1000 * (proc.time()[["elapsed"]] - t0)
  stop_for_status(r)
  txt <- content(r, as = "text", encoding = "UTF-8")
  out <- fromJSON(txt)
  attr(out, "raw") <- txt
  attr(out, "ms")  <- ms
  out
}

# Administers one complete test, answering each item from the examinee's responses.
brm_run_one <- function(testID, answer, max_calls = 60, base_url = ocat_base_url()) {
  administered <- data.frame()
  steps <- list(); ms <- numeric(0)
  d <- brm_post(testID, administered, NULL, base_url = base_url)
  ms <- c(ms, attr(d, "ms")); k <- 0
  while (identical(d$status, "OK") && !isTRUE(d$testCompleted) && k < max_calls) {
    k <- k + 1
    nxt <- d$nextItemsToAdminister
    if (is.null(nxt) || !NROW(nxt)) break
    nxt$itemScore <- as.integer(answer(nxt$itemID))
    administered <- rbind(administered, nxt)
    d <- brm_post(testID, administered, d$sessionData, base_url = base_url)
    ms <- c(ms, attr(d, "ms"))
    steps[[k]] <- data.frame(step = k,
                             slot = if (!is.null(nxt$slot)) nxt$slot else NA,
                             itemID = nxt$itemID, itemScore = nxt$itemScore,
                             theta_hat = num1(d$sessionData$test_theta),
                             se = num1(d$sessionData$test_se),
                             stringsAsFactors = FALSE)
  }
  if (!identical(d$status, "OK")) stop("engine status ", d$status, ": ", attr(d, "raw"))
  if (!isTRUE(d$testCompleted)) stop("test did not complete within ", max_calls, " calls")
  list(steps = do.call(rbind, steps), final = d, raw = attr(d, "raw"), ms = ms,
       administered = administered)
}

# Answer function for one examinee: the response to whichever item is administered.
brm_answer_fn <- function(resp_row, bank) {
  fixed <- unlist(resp_row[bank$itm_id]); names(fixed) <- bank$itm_id
  function(ids) { s <- unname(fixed[ids]); s[is.na(s)] <- -1L; s }
}

# EAP estimate and posterior SD on a grid (used by the self-test and pre-flight checks).
eap <- function(bank, ids, scores, sign = GPCM_SIGN, lb = -5, ub = 5, step = 0.005) {
  grid <- seq(lb, ub, by = step)
  ll <- dnorm(grid, log = TRUE)
  for (k in seq_along(ids)) {
    it <- bank_row(bank, ids[k]); if (is.null(it) || scores[k] < 0) next
    p <- item_probs(it, grid, sign)
    ll <- ll + log(p[, scores[k] + 1])
  }
  w <- exp(ll - max(ll)); w <- w / sum(w)
  m <- sum(w * grid)
  c(theta = m, se = sqrt(sum(w * (grid - m)^2)))
}

# Writes scores.csv in the simulator's layout.
.write_scores <- function(run_dir, res) {
  write.csv(res$scores, file.path(run_dir, "scores.csv"), row.names = FALSE, na = "")
}

# Writes constraints.csv: for each constraint, the distribution of counts across tests.
.write_constraints <- function(run_dir, adm_long) {
  cs <- read.csv(CONSTRAINT_SETS, stringsAsFactors = FALSE)
  sims <- sort(unique(adm_long$simulee))
  ids_by_sim <- split(adm_long$itemID, adm_long$simulee)
  rows <- lapply(unique(cs$constraint_id), function(cid) {
    s <- cs[cs$constraint_id == cid, ]
    cnt <- vapply(ids_by_sim, function(v) sum(v %in% s$itm_id), numeric(1))
    lb <- s$lb[1]; ub <- s$ub[1]
    dist <- tabulate(cnt + 1, nbins = 36); names(dist) <- 0:35
    data.frame(section = "A", constraint_id = cid, apply = "ITM", stm_id = "",
               type = s$type[1], lb = lb, ub = ub, violation_cost = "",
               N = length(cnt), on_target = round(mean(cnt >= lb & cnt <= ub), 3),
               min = min(cnt), max = max(cnt), mean = round(mean(cnt), 3),
               sd = round(sqrt(mean((cnt - mean(cnt))^2)), 3),
               t(dist), check.names = FALSE, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  write.csv(out, file.path(run_dir, "constraints.csv"), row.names = FALSE)
  pos <- unique(cs[nzchar(cs$positions) & !is.na(cs$positions), c("constraint_id", "positions")])
  if (nrow(pos) && "slot" %in% names(adm_long)) {
    pc <- do.call(rbind, lapply(seq_len(nrow(pos)), function(k) {
      cid <- pos$constraint_id[k]; members <- cs$itm_id[cs$constraint_id == cid]
      ps <- unlist(lapply(strsplit(pos$positions[k], ",")[[1]], function(x) {
        x <- trimws(x); if (grepl("-", x)) { r <- as.integer(strsplit(x, "-")[[1]]); r[1]:r[2] } else as.integer(x) }))
      ok <- vapply(split(adm_long, adm_long$simulee), function(a) {
        inm <- a$itemID %in% members
        all(a$slot[inm] %in% ps)
      }, logical(1))
      data.frame(constraint_id = cid, positions = pos$positions[k], N = length(ok),
                 all_members_in_position = round(mean(ok), 4))
    }))
    write.csv(pc, file.path(run_dir, "positions_check.csv"), row.names = FALSE)
  }
  out
}

# Writes exposure.csv: administrations and exposure rate for every item.
.write_exposure <- function(run_dir, adm_long, all_ids) {
  n <- length(unique(adm_long$simulee))
  cnt <- table(factor(adm_long$itemID, levels = all_ids))
  ex <- rbind(data.frame(type = "section", id = "A", count = n, exposure_rate = 1.0),
              data.frame(type = "item", id = all_ids, count = as.integer(cnt),
                         exposure_rate = round(as.integer(cnt) / n, 4)))
  write.csv(ex, file.path(run_dir, "exposure.csv"), row.names = FALSE)
  ex
}

# Offline checks of the design, the response functions, the response generator and
# the inputs. No contact with the engine.
brm_selftest <- function() {
  ok <- TRUE; say <- function(flag, msg) { cat(if (flag) "  OK   " else "  FAIL ", msg, "\n"); ok <<- ok && flag }
  th <- brm_thetas()
  say(length(th) == 3600 && all(table(th) == 100), "design: 36 quantiles x 100")
  say(abs(var(QUANTILES) * 35 / 36 - 0.9654) < 5e-4, sprintf("design variance %.4f (expect 0.9654)", var(QUANTILES) * 35 / 36))
  g <- list(model = "GPCM", sc = 1, a = 0.734, b = 1.381, c = 0, steps = c(1.597, -0.473, -0.776, -0.151))
  tt <- seq(-3, 3, by = 0.5)
  p <- item_probs(g, tt)
  say(max(abs(rowSums(p) - 1)) < 1e-12, "GPCM probabilities sum to 1")
  ex <- function(t) sum(item_probs(g, t) * 0:4)
  h <- 1e-5; num <- sapply(tt, function(t) (ex(t + h) - ex(t - h)) / (2 * h))
  say(max(abs(num * g$sc * g$a - item_info(g, tt))) < 1e-6, "GPCM information = sc*a * dE[X]/dtheta")
  two <- list(model = "2PL", sc = 1, a = 0.842, b = 0.768, c = 0, steps = numeric(0))
  P <- item_probs(two, tt)[, 2]
  say(max(abs(item_info(two, tt) - two$a^2 * P * (1 - P))) < 1e-12, "2PL information = a^2 PQ")
  bank <- data.frame(itm_id = c("X1", "X2", "X3"), model = c("2PL", "2PL", "GPCM"),
                     max_points = c(1L, 1L, 4L), sc = 1, a = c(1, 0.8, 0.734), b = c(0, 0.5, 1.381), c = 0)
  bank$steps <- list(numeric(0), numeric(0), g$steps)
  set.seed(1); big <- brm_make_responses(bank, rep(0.5, 20000), seed = 7)
  say(abs(mean(big$X1) - 1 / (1 + exp(-0.5))) < 0.01, "2PL responses match P(theta)")
  pg <- item_probs(g, 0.5)
  say(max(abs(as.numeric(prop.table(table(factor(big$X3, levels = 0:4)))) - as.numeric(pg))) < 0.01, "GPCM responses match category probabilities")
  e <- eap(bank, c("X1", "X2", "X3"), c(1, 0, 2))
  say(is.finite(e[["theta"]]) && e[["se"]] > 0 && e[["se"]] < 1, sprintf("EAP runs (%.3f, %.3f)", e[["theta"]], e[["se"]]))
  say(file.exists(CONSTRAINT_SETS), "constraint sets present (input/brm_constraint_sets.csv)")
  if (file.exists(CONSTRAINT_SETS)) {
    cs <- read.csv(CONSTRAINT_SETS, stringsAsFactors = FALSE)
    say(length(unique(cs$constraint_id)) == 231, sprintf("%d constraints (expect 231)", length(unique(cs$constraint_id))))
  }
  for (p in c("httr", "jsonlite")) say(requireNamespace(p, quietly = TRUE), paste("package", p))
  cat(if (ok) "\nself-test passed\n" else "\nSELF-TEST FAILED -- do not contact the engine yet\n")
  invisible(ok)
}

# Checks before the runs, on brmlim: the engine is reachable, the six tests are loaded,
# and the engine's scores for 36 examinees agree with an independent EAP computation.
brm_preflight <- function(n = 36, base_url = ocat_base_url()) {
  cat("engine:", base_url, "\n")
  root <- tryCatch(ocat_get("/", base_url), error = function(e)
    stop("cannot reach the engine at ", base_url, " -- is it running, is port 8080 open ",
         "on the Private profile, and is the address right? Try it in a browser.\n", conditionMessage(e)))
  txt <- toJSON(root, auto_unbox = TRUE)
  missing <- EXPECTED_TESTS[!vapply(EXPECTED_TESTS, function(t) grepl(sprintf('"%s"', t), txt), logical(1))]
  cat("loaded tests: ", txt, "\n")
  if (length(missing)) stop("not loaded: ", paste(missing, collapse = ", "))
  bank <- brm_bank(PROBE_ID, base_url)
  st <- attr(bank, "status")
  cat(sprintf("bank: %d OP (%s), %d FT\n", nrow(bank),
              paste(names(table(bank$model)), table(bank$model), collapse = ", "), sum(st == "FT")))
  stopifnot(nrow(bank) == 200, sum(bank$model == "GPCM") == 8)
  pf_dir <- new_run_dir(CELL, "preflight")
  capture_provenance(pf_dir, PROBE_ID, script_path = DRIVER_PATH, seed = SEED + 99L, base_url = base_url,
                     notes = "pre-flight: GPCM sign, FT handling, per-call time")
  set.seed(SEED + 99L)
  th <- QUANTILES[seq(1, 36, length.out = n)]
  resp <- brm_make_responses(bank, th, seed = SEED + 99L)
  rows <- list(); ms <- numeric(0)
  for (i in seq_len(n)) {
    r <- brm_run_one(PROBE_ID, brm_answer_fn(resp[i, ], bank), base_url = base_url)
    ms <- c(ms, r$ms)
    a <- r$administered
    ov <- score_overall(r$final)
    ep <- eap(bank, a$itemID, a$itemScore, sign = +1)
    em <- eap(bank, a$itemID, a$itemScore, sign = -1)
    n_ft <- sum(!(a$itemID %in% bank$itm_id)); n_g <- sum(a$itemID %in% bank$itm_id[bank$model == "GPCM"])
    rows[[i]] <- data.frame(true = th[i], engine = ov$theta, engine_se = ov$se,
                            eap_plus = ep[["theta"]], eap_minus = em[["theta"]],
                            se_plus = ep[["se"]], n_items = nrow(a), n_ft = n_ft, n_gpcm = n_g)
    if (i == 1) writeLines(r$raw, file.path(pf_dir, "first_final_response.json"))
  }
  x <- do.call(rbind, rows)
  dp <- max(abs(x$engine - x$eap_plus)); dm <- max(abs(x$engine - x$eap_minus))
  dse <- max(abs(x$engine_se - x$se_plus))
  write.csv(x, file.path(pf_dir, "preflight.csv"), row.names = FALSE)
  cat(sprintf("\nitems per test: %s; FT per test: %s; GPCM per test: %s\n",
              paste(range(x$n_items), collapse = "-"), paste(range(x$n_ft), collapse = "-"),
              paste(range(x$n_gpcm), collapse = "-")))
  cat(sprintf("max |engine - EAP|, GPCM sign +: %.5f    sign -: %.5f    (SE, sign +: %.5f)\n", dp, dm, dse))
  cat(sprintf("POST /cat round trip: median %.1f ms, 90th pct %.1f ms, over %d calls\n",
              median(ms), quantile(ms, 0.9), length(ms)))
  per_run_h <- median(ms) * mean(x$n_items + 1) * 3600 / 3.6e6
  cat(sprintf("=> roughly %.1f h per 3,600-simulee condition at this speed (condition 4 will be slower)\n", per_run_h))
  best <- if (dp <= dm) +1 else -1
  if (abs(dp - dm) < 0.002)
    warning("the two GPCM signs give nearly the same theta-hat here (", round(dp, 4), " vs ", round(dm, 4),
            "); the sign test is inconclusive -- report this before running the conditions")
  if (min(dp, dm) > 0.01)
    warning("neither sign reproduces the engine's theta-hat to 0.01 -- check FT handling and the ",
            "estimation settings before running anything")
  if (best != GPCM_SIGN)
    stop("the engine matches GPCM sign ", best, " but GPCM_SIGN is ", GPCM_SIGN,
         ". Set GPCM_SIGN <- ", best, " at the top of this file and re-source it.")
  close_provenance(pf_dir, PROBE_ID, base_url)
  cat("\npre-flight passed: GPCM sign", GPCM_SIGN, "confirmed; FT items excluded from theta.\n")
  invisible(x)
}

# Theoretical limit on the SE at each theta: the information of the first shadow
# test assembled at that theta (brmlim, no exposure control).
brm_probe_limit <- function(grid = seq(-3, 3, by = 0.25), base_url = ocat_base_url()) {
  bank <- brm_bank(PROBE_ID, base_url)
  run_dir <- new_run_dir(CELL, "brmlim_limit")
  capture_provenance(run_dir, PROBE_ID, script_path = DRIVER_PATH, seed = NA, base_url = base_url,
                     notes = "theoretical-limit probes: first shadow test at a fixed start theta")
  out <- do.call(rbind, lapply(grid, function(t) {
    d <- brm_post(PROBE_ID, data.frame(), NULL, start_theta = t, base_url = base_url)
    ids <- d$sessionData$shadowTest$itemID
    op <- ids[ids %in% bank$itm_id]
    info <- sum(vapply(op, function(id) item_info(bank_row(bank, id), t), numeric(1)))
    data.frame(theta = t, n_items = length(ids), n_op = length(op), info = info,
               limit_sem = 1 / sqrt(info), ms = attr(d, "ms"))
  }))
  write.csv(out, file.path(run_dir, "theoretical_limit.csv"), row.names = FALSE)
  close_provenance(run_dir, PROBE_ID, base_url)
  cat("written:", run_dir, "\n")
  options(brm.limit_file = file.path(run_dir, "theoretical_limit.csv"))
  invisible(out)
}

# Runs one condition for all 3,600 examinees and writes its results folder.
brm_run <- function(testID, n = NULL, resume = NULL, base_url = ocat_base_url()) {
  stopifnot(testID %in% c(CONDITIONS, "brm4w4", "brm4w8", "brm4c", "brm4cs"))
  bank <- brm_bank(testID, base_url)
  thetas <- brm_thetas()
  resp <- brm_responses(bank, thetas)
  N <- if (is.null(n)) nrow(resp) else min(n, nrow(resp))
  if (is.null(resume)) {
    run_dir <- new_run_dir(CELL, testID)
    capture_provenance(run_dir, testID, script_path = DRIVER_PATH, seed = SEED, base_url = base_url,
                       notes = sprintf("R driver; %d simulees; responses %s; GPCM sign %+d",
                                       N, basename(attr(resp, "path")), GPCM_SIGN))
    writeLines(unname(tools::md5sum(attr(resp, "path"))), file.path(run_dir, "responses_used.md5"))
    acc <- list(scores = list(), items = list(), timing = list(), done = 0L)
  } else {
    run_dir <- resume
    acc <- readRDS(file.path(run_dir, "progress.rds"))
    message("resuming ", basename(run_dir), " after simulee ", acc$done)
  }
  lim <- getOption("brm.limit_file")
  if (!is.null(lim) && file.exists(lim)) file.copy(lim, file.path(run_dir, "theoretical_limit.csv"), overwrite = TRUE)
  t_start <- Sys.time()
  for (i in seq.int(acc$done + 1L, N)) {
    r <- tryCatch(brm_run_one(testID, brm_answer_fn(resp[i, ], bank), base_url = base_url),
                  error = function(e) { saveRDS(acc, file.path(run_dir, "progress.rds"))
                    stop("stopped at simulee ", i, ": ", conditionMessage(e),
                         "\nResume with brm_run(\"", testID, "\", resume = \"", run_dir, "\")", call. = FALSE) })
    ov <- score_overall(r$final)
    if (!NROW(ov)) {
      writeLines(r$raw, file.path(run_dir, sprintf("unreadable_score_simulee%d.json", i)))
      saveRDS(acc, file.path(run_dir, "progress.rds"))
      stop("no overall score block in the final response for simulee ", i, "; saved it as ",
           sprintf("unreadable_score_simulee%d.json", i), call. = FALSE)
    }
    if (!is.finite(ov$theta) || abs(ov$theta) >= 4.9 || ov$raw %in% c(0, ov$max_pts))
      writeLines(r$raw, file.path(run_dir, sprintf("extreme_score_simulee%d.json", i)))
    bth <- r$final$score$bayesian_theta
    bt    <- if (is.list(bth)) num1(bth$mean) else num1(bth)
    bt_sd <- if (is.list(bth)) num1(bth$sd)   else NA_real_
    acc$scores[[i]] <- data.frame(thread = 1L, student = i, domain_name = "Overall", domain_value = "Overall",
                                  item_count = ov$n_items, raw_score = ov$raw, max_points = ov$max_pts,
                                  true_theta = thetas[i], theta = ov$theta, se = ov$se,
                                  scale_score = ov$scale, scale_score_se = ov$scale_se,
                                  performance_level = NA, bayesian_theta = bt,
                                  bayesian_sd = bt_sd)
    st <- r$steps; st$simulee <- i
    acc$items[[i]] <- st
    acc$timing[[i]] <- data.frame(simulee = i, call = seq_along(r$ms), ms = r$ms)
    if (i == 1) writeLines(r$raw, file.path(run_dir, "first_final_response.json"))
    acc$done <- i
    if (i %% 100 == 0 || i == N) {
      saveRDS(acc, file.path(run_dir, "progress.rds"))
      el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      message(sprintf("%s  %d/%d  %.1f min elapsed", testID, i, N, el))
    }
  }
  scores <- do.call(rbind, acc$scores)
  items  <- do.call(rbind, acc$items)
  timing <- do.call(rbind, acc$timing)
  .write_scores(run_dir, list(scores = scores))
  write.csv(items, file.path(run_dir, "items_long.csv"), row.names = FALSE)
  write.csv(timing, file.path(run_dir, "timing.csv"), row.names = FALSE)
  med <- median(timing$ms[timing$call > 1])
  writeLines(sprintf("%.1f ms median POST /cat round trip (calls after the first; %d calls), R on %s -> engine %s",
                     med, sum(timing$call > 1), Sys.info()[["nodename"]], base_url),
             file.path(run_dir, "runtime.txt"))
  .write_constraints(run_dir, items)
  .write_exposure(run_dir, items, attr(bank, "all_ids"))
  close_provenance(run_dir, testID, base_url)
  unlink(file.path(run_dir, "progress.rds"))
  if (any(is.finite(scores$bayesian_theta))) {
    cat(sprintf("Bayesian mode: max |overall theta - posterior mean| = %.5f; max |overall se - posterior sd| = %.5f\n",
                max(abs(scores$theta - scores$bayesian_theta), na.rm = TRUE),
                max(abs(scores$se - scores$bayesian_sd), na.rm = TRUE)))
  }
  ext <- with(scores, !is.finite(theta) | abs(theta) >= 4.9 | raw_score %in% 0 | raw_score == max_points)
  cat(sprintf("\n%s done: %d simulees, %s\n", testID, N, run_dir))
  cat(sprintf("extreme final scores (raw 0 or max, |theta-hat| >= 4.9, or non-finite): %d of %d%s\n",
              sum(ext, na.rm = TRUE), N,
              if (any(ext, na.rm = TRUE)) " -- see extreme_score_simulee*.json" else ""))
  invisible(run_dir)
}

# Scores a given set of administered items and responses (POST /rescore).
brm_rescore_post <- function(testID, adm, base_url = ocat_base_url()) {
  body <- toJSON(list(testID = testID, administeredItems = adm, sessionData = NULL),
                 auto_unbox = TRUE, digits = NA)
  t0 <- proc.time()[["elapsed"]]
  r <- POST(paste0(base_url, "/rescore"), body = body, encode = "raw",
            add_headers(accept = "application/json"), content_type("application/json"))
  ms <- 1000 * (proc.time()[["elapsed"]] - t0)
  stop_for_status(r)
  out <- fromJSON(content(r, as = "text", encoding = "UTF-8"))
  attr(out, "ms") <- ms
  out
}

# Most recent completed results folder for a test.
.latest_run <- function(testID) {
  d <- list.dirs(file.path(CELL_DIR, "results"), recursive = FALSE)
  d <- d[grepl(paste0("_", testID, "$"), d) & file.exists(file.path(d, "COMPLETED.txt"))]
  if (!length(d)) stop("no completed run folder for ", testID)
  sort(d, decreasing = TRUE)[1]
}

# Students A (theta = 1.5) and B (theta = -1.5) on Condition 4: the posterior mean,
# SD and draws after each item (Figures 5 and 6).
brm_capture_students <- function(testID = "brm4cs", thetas = c(A = 1.5, B = -1.5),
                                 seed = 20261031L, base_url = ocat_base_url()) {
  bank <- brm_bank(testID, base_url)
  resp <- brm_make_responses(bank, unname(thetas), seed = seed)
  out_dir <- new_run_dir(CELL, paste0(testID, "_students"))
  capture_provenance(out_dir, testID, script_path = DRIVER_PATH, seed = seed, base_url = base_url,
                     notes = sprintf("posterior trajectories: %s", paste(names(thetas), thetas, collapse = ", ")))
  write.csv(resp, file.path(out_dir, "students_responses.csv"), row.names = FALSE)
  post_of <- function(d) {
    b <- d$sessionData$bayesian_theta
    if (is.null(b) || !length(b)) b <- d$score$bayesian_theta
    v <- if (is.null(b)) numeric(0) else as.numeric(unlist(b$vector))
    list(mean = if (is.null(b)) NA_real_ else num1(b$mean),
         sd = if (is.null(b)) NA_real_ else num1(b$sd), v = v)
  }
  summ <- list(); draws <- list()
  for (j in seq_along(thetas)) {
    lab <- names(thetas)[j]; ans <- brm_answer_fn(resp[j, ], bank)
    administered <- data.frame(); k <- 0
    d <- brm_post(testID, administered, NULL, base_url = base_url)
    repeat {
      writeLines(attr(d, "raw"), file.path(out_dir, sprintf("raw_%s_%02d.json", lab, k)))
      po <- post_of(d)
      n_scored <- if (NROW(administered)) sum(administered$itemScore >= 0) else 0L
      summ[[length(summ) + 1]] <- data.frame(student = lab, true_theta = thetas[[j]], position = k,
        n_scored = n_scored, item = if (k > 0) administered$itemID[k] else "",
        score = if (k > 0) administered$itemScore[k] else NA, post_mean = po$mean, post_sd = po$sd,
        n_draws = length(po$v), stringsAsFactors = FALSE)
      if (length(po$v)) draws[[length(draws) + 1]] <- data.frame(student = lab, position = k,
        matrix(po$v, nrow = 1, dimnames = list(NULL, paste0("v", seq_along(po$v)))), stringsAsFactors = FALSE)
      if (!identical(d$status, "OK") || isTRUE(d$testCompleted) || k >= 60) break
      nxt <- d$nextItemsToAdminister
      if (is.null(nxt) || !NROW(nxt)) break
      nxt$itemScore <- as.integer(ans(nxt$itemID))
      administered <- rbind(administered, nxt); k <- k + NROW(nxt)
      d <- brm_post(testID, administered, d$sessionData, base_url = base_url)
    }
  }
  x <- do.call(rbind, summ)
  write.csv(x, file.path(out_dir, "students_trajectory.csv"), row.names = FALSE)
  if (length(draws)) {
    dr <- do.call(rbind, draws)
    write.csv(dr, file.path(out_dir, "students_draws.csv"), row.names = FALSE)
  }
  close_provenance(out_dir, testID, base_url)
  cat("\nFinal posterior (mean, SD) and draws captured per position:\n")
  print(x[x$position == ave(x$position, x$student, FUN = max), c("student", "true_theta", "position", "n_scored", "post_mean", "post_sd")], row.names = FALSE)
  cat(sprintf("positions with draws: %d of %d\n", sum(x$n_draws > 0), nrow(x)))
  cat("written:", out_dir, "\n")
  invisible(x)
}

# Category probabilities averaged over each item's parameter draws, on a grid.
.item_avg_probs <- function(testID, grid, base_url = ocat_base_url()) {
  td <- ocat_get(paste0("/", testID), base_url)
  ib <- td$item_bank; bay <- ib$bayesian
  op <- which(ib$status == "OP")
  get <- function(f, i) { x <- bay[[f]]; if (is.list(x)) x[[i]] else x[i, ] }
  out <- list()
  for (i in op) {
    a <- as.numeric(unlist(get("a_vector", i))); b <- as.numeric(unlist(get("b_vector", i)))
    sc <- if (is.null(ib$sc)) 1 else as.numeric(ib$sc[i]); if (!is.finite(sc)) sc <- 1
    if (toupper(ib$model[i]) %in% c("GPCM", "GPC")) {
      dv <- get("d_vector", i)
      D <- if (is.list(dv)) do.call(cbind, lapply(dv, function(v) as.numeric(unlist(v)))) else as.matrix(dv)
      if (nrow(D) != length(a)) D <- t(D)
      m <- ncol(D); acc <- matrix(0, m + 1, length(grid))
      for (s in seq_along(a)) {
        z <- rbind(0, apply(sc * a[s] * outer(D[s, ], grid, function(dk, th) th - b[s] + dk), 2, cumsum))
        z <- sweep(z, 2, apply(z, 2, max)); e <- exp(z); acc <- acc + sweep(e, 2, colSums(e), "/")
      }
      out[[ib$itm_id[i]]] <- acc / length(a)
    } else {
      P <- 1 / (1 + exp(-sc * a * outer(-b, grid, "+")))
      p1 <- colMeans(P); out[[ib$itm_id[i]]] <- rbind(1 - p1, p1)
    }
  }
  attr(out, "prior") <- as.numeric(unlist(td$test_settings$bayesian_theta$vector))
  out
}

# Mean and SD of a posterior given on a grid.
.moments <- function(grid, logw) { w <- exp(logw - max(logw)); w <- w / sum(w)
  m <- sum(w * grid); c(mean = m, sd = sqrt(sum(w * (grid - m)^2))) }

# Equally spaced quantiles of a posterior given on a grid, as draws.
.quantile_draws <- function(grid, logw, n = 500) {
  w <- exp(logw - max(logw)); cw <- cumsum(w) / sum(w)
  keep <- !duplicated(cw)
  approx(cw[keep], grid[keep], xout = ((1:n) - 0.5) / n, rule = 2)$y
}

# Rescores with a given set of draws as the prior.
brm_rescore_prior <- function(testID, adm, prior_draws, base_url = ocat_base_url()) {
  body <- toJSON(list(testID = testID, administeredItems = adm,
                      sessionData = list(bayesian_theta = list(vector = prior_draws))),
                 auto_unbox = TRUE, digits = NA)
  r <- POST(paste0(base_url, "/rescore"), body = body, encode = "raw",
            add_headers(accept = "application/json"), content_type("application/json"))
  stop_for_status(r)
  fromJSON(content(r, as = "text", encoding = "UTF-8"))
}

# Posterior mean and SD from an engine response.
.post_of <- function(d) {
  b <- d$score$bayesian_theta; if (is.null(b)) b <- d$sessionData$bayesian_theta
  v <- as.numeric(unlist(b$vector)); c(mean = mean(v), sd = sd(v))
}

# Condition 4: the engine's posterior against the exact posterior, after a single
# update from an exact prior and after the full sequence of updates.
brm_update_check <- function(testID = "brm4cs", run_dir = .latest_run("brm4cs"),
                             per_q = 2, ks_A = c(1, 2, 3, 5, 8, 12, 17, 23, 30, 35), reps_A = 2,
                             ks_B = c(5, 15, 30), reps_B = 5, seed = 20261101L,
                             base_url = ocat_base_url()) {
  set.seed(seed)
  grid <- seq(-7, 7, by = 0.01)
  AP <- .item_avg_probs(testID, grid, base_url)
  pd <- attr(AP, "prior"); m0 <- mean(pd); s0 <- sd(pd)
  lprior <- dnorm(grid, m0, s0, log = TRUE)
  s  <- read.csv(file.path(run_dir, "scores.csv"))
  it <- read.csv(file.path(run_dir, "items_long.csv"), stringsAsFactors = FALSE)
  pick <- unlist(lapply(split(s$student, s$true_theta), function(v) head(sort(v), per_q)))
  out_dir <- new_run_dir(CELL, paste0(testID, "_updatecheck"))
  capture_provenance(out_dir, testID, script_path = DRIVER_PATH, seed = seed, base_url = base_url,
                     notes = sprintf("update check on %s: %d simulees; A ks %s x%d; B ks %s x%d",
                                     basename(run_dir), length(pick), paste(ks_A, collapse = ","), reps_A,
                                     paste(ks_B, collapse = ","), reps_B))
  t0 <- Sys.time(); A <- list(); B <- list()
  for (k_sim in pick) {
    a <- it[it$simulee == k_sim, ]; a <- a[order(a$step), ]
    adm_all <- data.frame(slot = a$slot, sectionID = "A", itemID = a$itemID, stimulusID = "",
                          itemScore = as.integer(a$itemScore), stringsAsFactors = FALSE)
    sc_rows <- which(adm_all$itemScore >= 0 & adm_all$itemID %in% names(AP))
    tt <- s$true_theta[s$student == k_sim]
    lex <- lprior; adf <- c(m0, s0); EX <- list(); AD <- list(); LW <- list()
    LW[[1]] <- lex
    for (j in seq_along(sc_rows)) {
      r <- sc_rows[j]; ll <- log(AP[[adm_all$itemID[r]]][adm_all$itemScore[r] + 1, ])
      lex <- lex + ll; EX[[j]] <- .moments(grid, lex); LW[[j + 1]] <- lex
      lad <- dnorm(grid, adf[1], adf[2], log = TRUE) + ll; adf <- .moments(grid, lad); AD[[j]] <- adf
    }
    for (k in ks_A[ks_A <= length(sc_rows)]) {
      adm <- adm_all[seq_len(sc_rows[k]), ]
      for (rp in seq_len(reps_A)) {
        e <- tryCatch(.post_of(brm_rescore_post(testID, adm, base_url)), error = function(err) c(mean = NA, sd = NA))
        A[[length(A) + 1]] <- data.frame(student = k_sim, true_theta = tt, k = k, rep = rp,
          eng_mean = e[["mean"]], eng_sd = e[["sd"]], exact_mean = EX[[k]][["mean"]], exact_sd = EX[[k]][["sd"]],
          adf_mean = AD[[k]][["mean"]], adf_sd = AD[[k]][["sd"]])
      }
    }
    for (k in ks_B[ks_B <= length(sc_rows)]) {
      draws <- .quantile_draws(grid, LW[[k]])
      r <- sc_rows[k]; ll <- log(AP[[adm_all$itemID[r]]][adm_all$itemScore[r] + 1, ])
      nm <- dnorm(grid, mean(draws), sd(draws), log = TRUE) + ll
      ek <- EX[[k]]; nk <- .moments(grid, nm); pk <- .moments(grid, LW[[k]])
      for (rp in seq_len(reps_B)) {
        e <- tryCatch(.post_of(brm_rescore_prior(testID, adm_all[r, , drop = FALSE], draws, base_url)),
                      error = function(err) c(mean = NA, sd = NA))
        B[[length(B) + 1]] <- data.frame(student = k_sim, true_theta = tt, k = k, rep = rp,
          prior_mean = pk[["mean"]], prior_sd = pk[["sd"]], eng_mean = e[["mean"]], eng_sd = e[["sd"]],
          exact_mean = ek[["mean"]], exact_sd = ek[["sd"]], normal_mean = nk[["mean"]], normal_sd = nk[["sd"]])
      }
    }
  }
  A <- do.call(rbind, A); B <- do.call(rbind, B)
  write.csv(A, file.path(out_dir, "partA_tracking.csv"), row.names = FALSE)
  write.csv(B, file.path(out_dir, "partB_onestep.csv"), row.names = FALSE)
  close_provenance(out_dir, testID, base_url)
  reg <- function(t) cut(t, c(-Inf, -1.6, 1.6, Inf), labels = c("low", "mid", "high"))
  sA <- do.call(rbind, lapply(split(A, list(A$k, reg(A$true_theta)), drop = TRUE), function(g) data.frame(
    k = g$k[1], region = as.character(reg(g$true_theta[1])), n = nrow(g),
    eng_exact_mean = mean(g$eng_mean - g$exact_mean, na.rm = TRUE), eng_exact_sd = sd(g$eng_mean - g$exact_mean, na.rm = TRUE),
    sd_ratio = mean(g$eng_sd / g$exact_sd, na.rm = TRUE),
    adf_exact_mean = mean(g$adf_mean - g$exact_mean), adf_sd_ratio = mean(g$adf_sd / g$exact_sd))))
  sA <- sA[order(sA$region, sA$k), ]; sA[-(1:2)] <- round(sA[-(1:2)], 4)
  sB <- do.call(rbind, lapply(split(B, list(B$k, reg(B$true_theta)), drop = TRUE), function(g) data.frame(
    k = g$k[1], region = as.character(reg(g$true_theta[1])), n = nrow(g),
    eng_exact_mean = mean(g$eng_mean - g$exact_mean, na.rm = TRUE), eng_exact_sd = sd(g$eng_mean - g$exact_mean, na.rm = TRUE),
    rep_sd = mean(tapply(g$eng_mean, g$student, sd), na.rm = TRUE),
    sd_ratio_exact = mean(g$eng_sd / g$exact_sd, na.rm = TRUE), sd_ratio_normal = mean(g$eng_sd / g$normal_sd, na.rm = TRUE),
    eng_normal_mean = mean(g$eng_mean - g$normal_mean, na.rm = TRUE))))
  sB <- sB[order(sB$region, sB$k), ]; sB[-(1:2)] <- round(sB[-(1:2)], 4)
  write.csv(sA, file.path(out_dir, "partA_summary.csv"), row.names = FALSE)
  write.csv(sB, file.path(out_dir, "partB_summary.csv"), row.names = FALSE)
  cat(sprintf("\n%d simulees; %.1f min. NA rescoring results: A %d, B %d\n", length(pick),
              as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(is.na(A$eng_mean)), sum(is.na(B$eng_mean))))
  cat("\nPart A, engine (default prior, first k items) vs exact and vs deterministic normal re-approximation:\n"); print(sA, row.names = FALSE)
  cat("\nPart B, one update from the exact posterior after k-1 items:\n"); print(sB, row.names = FALSE)
  cat("written:", out_dir, "\n")
  invisible(list(A = A, B = B))
}

# Condition 4 rescored with a larger posterior sample retained at each update.
# Requires bayesian.post_burn_in = 50000 in the engine's application.properties
# (thinning 10: 5,000 draws per update) and an engine restart.
brm_rescore_bigS <- function(testID = "brm4cs", run_dir = .latest_run("brm4cs"), expected_S = 5000L,
                             label = "bigS", resume = FALSE, base_url = ocat_base_url()) {
  s  <- read.csv(file.path(run_dir, "scores.csv"))
  s  <- s[s$domain_name == "Overall", ]
  it <- read.csv(file.path(run_dir, "items_long.csv"), stringsAsFactors = FALSE)
  adm_of <- function(k) { a <- it[it$simulee == k, ]; a <- a[order(a$step), ]
    data.frame(slot = a$slot, sectionID = "A", itemID = a$itemID, stimulusID = "",
               itemScore = as.integer(a$itemScore), stringsAsFactors = FALSE) }
  t1 <- proc.time()[["elapsed"]]
  d <- brm_rescore_post(testID, adm_of(s$student[1]), base_url)
  v <- as.numeric(unlist(d$score$bayesian_theta$vector))
  cat(sprintf("Step 1: returned posterior sample size %d (expected %d); one rescoring took %.2f s\n",
              length(v), expected_S, proc.time()[["elapsed"]] - t1))
  if (length(v) != expected_S)
    stop("sample size is not ", expected_S, ": check application.properties and the restart; if both are right, ",
         "the sample size is tied to the draws file and the fallback is a file with repeated item draws")
  out_dir <- if (resume) sort(list.dirs(file.path(CELL_DIR, "results"), recursive = FALSE)[
    grepl(paste0("_", testID, "_", label, "$"), list.dirs(file.path(CELL_DIR, "results"), recursive = FALSE))],
    decreasing = TRUE)[1] else new_run_dir(CELL, paste0(testID, "_", label))
  f <- file.path(out_dir, "rescore_bigS.csv")
  if (!resume) {
    capture_provenance(out_dir, testID, script_path = DRIVER_PATH, base_url = base_url,
                       notes = sprintf("rescoring all %s tests with %d posterior draws per update", basename(run_dir), expected_S))
    write.csv(data.frame(student = integer(0), true_theta = numeric(0), orig_mean = numeric(0), orig_sd = numeric(0),
                         new_mean = numeric(0), new_sd = numeric(0), n_draws = integer(0), ms = numeric(0)), f, row.names = FALSE)
  }
  done <- read.csv(f)$student
  todo <- setdiff(s$student, done)
  t0 <- Sys.time(); buf <- list()
  for (j in seq_along(todo)) {
    k <- todo[j]; t1 <- proc.time()[["elapsed"]]
    d <- brm_rescore_post(testID, adm_of(k), base_url)
    v <- as.numeric(unlist(d$score$bayesian_theta$vector))
    buf[[length(buf) + 1]] <- data.frame(student = k, true_theta = s$true_theta[s$student == k],
      orig_mean = s$theta[s$student == k], orig_sd = s$se[s$student == k],
      new_mean = mean(v), new_sd = sd(v), n_draws = length(v), ms = 1000 * (proc.time()[["elapsed"]] - t1))
    if (length(buf) == 100 || j == length(todo)) {
      write.table(do.call(rbind, buf), f, sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE); buf <- list()
      cat(sprintf("  %d/%d  %.1f min\n", length(done) + j, nrow(s), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      flush.console()
    }
  }
  x <- read.csv(f)
  cov <- function(m, sd, t) mean(abs(m - t) <= 1.96 * sd)
  reg <- cut(x$true_theta, c(-Inf, -1.6, -0.6, 0.6, 1.6, Inf))
  summ <- function(g) data.frame(n = nrow(g),
    rmse_orig = sqrt(mean((g$orig_mean - g$true_theta)^2)), rmse_new = sqrt(mean((g$new_mean - g$true_theta)^2)),
    sd_orig = mean(g$orig_sd), sd_new = mean(g$new_sd),
    cov_orig = cov(g$orig_mean, g$orig_sd, g$true_theta), cov_new = cov(g$new_mean, g$new_sd, g$true_theta),
    bias_orig = mean(g$orig_mean - g$true_theta), bias_new = mean(g$new_mean - g$true_theta))
  tab <- rbind(cbind(region = "all", summ(x)),
               do.call(rbind, lapply(split(x, reg), function(g) cbind(region = as.character(cut(g$true_theta[1], c(-Inf, -1.6, -0.6, 0.6, 1.6, Inf))), summ(g)))))
  tab[-1] <- round(tab[-1], 4)
  write.csv(tab, file.path(out_dir, "rescore_bigS_summary.csv"), row.names = FALSE)
  close_provenance(out_dir, testID, base_url)
  cat(sprintf("\n%d tests, posterior draws per test %s, median %.0f ms per rescoring\n", nrow(x),
              paste(unique(x$n_draws), collapse = ","), median(x$ms)))
  print(tab, row.names = FALSE)
  cat("Reference (exact posterior on the same tests): RMSE 0.390, mean SD 0.396, coverage 0.959\n")
  cat("written:", out_dir, "\n")
  invisible(x)
}

# =====================================================================
# RUN -- one step at a time, with this folder as the working directory
# =====================================================================
if (FALSE) {
  options(ocat.host = "localhost")   # or the address of the computer running the engine
  source("brm_driver.R")
  brm_selftest()                 # offline checks
  brm_preflight()                # engine reachable, tests loaded, scoring checked (brmlim)
  brm_probe_limit()              # theoretical limit on the SE (Figure 4)
  brm_run("brm1")                # the four conditions, in this order, once each (Table 2)
  brm_run("brm2")
  brm_run("brm3")
  brm_run("brm4cs")
  brm_capture_students()         # Students A and B (Figures 5 and 6)
  brm_update_check()             # Condition 4: single update vs accumulated updates
  # set bayesian.post_burn_in = 50000 in application.properties, restart the engine, then:
  brm_rescore_bigS()             # Condition 4 with 5,000 draws per update
  # ...then restore the default and restart.
}
# Exposure-control rates accumulate per test from the time the engine starts: run each
# condition once, after a restart, and make no trial runs on brm1-brm4cs.
