# ocat_common.R -- functions for calling the OptimalCAT API from R, used by brm_driver.R.

library(httr)
library(jsonlite)

# Folder handling: results are written to <folder>/results/<time stamp>[_label]/.
ocat_project_root <- function(start = getwd()) {
  p <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(p, "OptimalCAT.Rproj"))) return(p)
    parent <- dirname(p)
    if (parent == p) stop("OptimalCAT.Rproj not found above: ", start)
    p <- parent
  }
}
ocat_cell_dir <- function(cell) getOption("ocat.cell_dir", file.path(ocat_project_root(), cell))
iso_stamp <- function(t = Sys.time()) format(t, "%Y-%m-%d_%H%M%S")
new_run_dir <- function(cell, label = NULL) {
  stamp <- iso_stamp()
  if (!is.null(label) && nzchar(label)) stamp <- paste0(stamp, "_", label)
  d <- file.path(ocat_cell_dir(cell), "results", stamp)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# Engine address: set options(ocat.host = ..., ocat.port = ...) to override.
OCAT_DEFAULT_HOST <- "localhost"
OCAT_DEFAULT_PORT <- 8080
ocat_base_url <- function(host = getOption("ocat.host", OCAT_DEFAULT_HOST),
                          port = getOption("ocat.port", OCAT_DEFAULT_PORT)) {
  paste0("http://", host, ":", port)
}

# GET requests to the engine.
ocat_get_text <- function(path, base_url = ocat_base_url()) {
  r <- GET(paste0(base_url, path), add_headers(accept = "application/json"))
  stop_for_status(r)
  content(r, as = "text", encoding = "UTF-8")
}
ocat_get <- function(path, base_url = ocat_base_url()) {
  fromJSON(ocat_get_text(path, base_url))
}

# Score fields from an engine response.
num1 <- function(x) if (is.null(x) || !length(x)) NA_real_ else as.numeric(x)[1]
score_overall <- function(d) {
  ov <- d$score$overall
  if (is.null(ov)) return(data.frame())
  data.frame(n_items  = num1(ov$itemCount),
             raw      = num1(ov$rawScore),
             max_pts  = num1(ov$maxPoints),
             theta    = num1(ov$theta),
             se       = num1(ov$se),
             scale    = num1(ov$scaleScore),
             scale_se = num1(ov$scaleScore_se))
}

# Run records: run_meta.json (engine, R and package versions, seed), sessionInfo.txt,
# a copy of the script, and the engine's exposure-control state before and after the run.
capture_provenance <- function(run_dir, testID, script_path = NULL,
                               seed = NA_integer_, notes = "",
                               base_url = ocat_base_url()) {
  w <- function(txt, f) writeLines(txt, file.path(run_dir, f), useBytes = TRUE)
  root <- ocat_get_text("/", base_url)
  w(root, "engine_root.json")
  tdef <- ocat_get_text(paste0("/", testID), base_url)
  w(tdef, paste0("testdef_", testID, ".json"))
  iec <- tryCatch(ocat_get_text(paste0("/IEC/", testID), base_url),
                  error = function(e) paste0('{"error":"', conditionMessage(e), '"}'))
  w(iec, "IEC_pre.json")
  meta <- list(run_started = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
               testID = testID, engine_url = base_url, seed = seed,
               notes = notes,
               engine_root = fromJSON(root),
               r_version = R.version.string,
               platform = R.version$platform,
               host = Sys.info()[["nodename"]],
               packages = sapply(c("httr", "jsonlite"),
                                 function(p) as.character(packageVersion(p))))
  write_json(meta, file.path(run_dir, "run_meta.json"),
             auto_unbox = TRUE, pretty = TRUE)
  capture.output(sessionInfo(), file = file.path(run_dir, "sessionInfo.txt"))
  if (!is.null(script_path) && file.exists(script_path))
    file.copy(script_path, file.path(run_dir, basename(script_path)), overwrite = TRUE)
  invisible(run_dir)
}
close_provenance <- function(run_dir, testID, base_url = ocat_base_url()) {
  iec <- tryCatch(ocat_get_text(paste0("/IEC/", testID), base_url),
                  error = function(e) paste0('{"error":"', conditionMessage(e), '"}'))
  writeLines(iec, file.path(run_dir, "IEC_post.json"), useBytes = TRUE)
  writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
             file.path(run_dir, "COMPLETED.txt"))
  invisible(run_dir)
}
