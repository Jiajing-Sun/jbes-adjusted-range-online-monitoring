# ============================================================
# 00_utils.R  -- small utilities used across scripts
# ============================================================

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(..., "\n")
}

ensure_pkgs <- function(pkgs, auto_install = TRUE) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    if (!auto_install) stop("Missing packages: ", paste(missing, collapse = ", "))
    msg("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
  invisible(TRUE)
}

# Try to locate the directory of the currently running script.
script_dir <- function() {
  # Case 1: running via Rscript --file=...
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  m <- grep(file_arg, args)
  if (length(m) > 0) {
    path <- sub(file_arg, "", args[m[1]])
    return(dirname(normalizePath(path)))
  }
  # Case 2: sourced in interactive sessions
  ofile <- tryCatch(attr(sys.frame(1), "ofile"), error = function(e) NULL)
  if (!is.null(ofile)) {
    return(dirname(normalizePath(ofile)))
  }
  # Fallback
  return(getwd())
}

# Project root: R source files live directly at the root
project_root <- function() {
  d <- script_dir()
  base <- basename(d)
  if (base %in% c("scripts", "R")) {
    return(normalizePath(file.path(d, "..")))
  }
  return(normalizePath(d))
}

# Safe matrix inverse with ridge regularization
inv_pd <- function(A, ridge = 1e-10) {
  A2 <- A + diag(ridge, nrow(A))
  solve(A2)
}

# First index where x > c, else returns length(x)+1
first_exceed <- function(x, cval) {
  idx <- which(x > cval)
  if (length(idx) == 0) return(length(x) + 1L)
  return(min(idx))
}

# Detection summaries for a single vector of stopping times. Here k_star is
# the last pre-change monitoring index; the post-change sample starts at
# k_star + 1 in the simulation scripts.
detection_metrics <- function(stops, k_star, mT) {
  stops <- as.integer(stops)
  pre_alarm <- stops <= k_star
  horizon_reject <- stops <= mT
  post_detect <- stops > k_star & stops <= mT

  post_delay <- stops[post_detect] - k_star

  c(
    DP_horizon = mean(horizon_reject),
    PFA_pre = mean(pre_alarm),
    DP_post = mean(post_detect),
    CADD = if (any(post_detect)) mean(post_delay) else NA_real_,
    ARL_excess = mean(pmax(0, stops - k_star)),
    ADD_excess = if (any(horizon_reject)) mean(pmax(0, stops[horizon_reject] - k_star)) else NA_real_
  )
}

# Build a compact worker export bundle that includes the local closure state
# together with any globals recursively referenced by those functions.
build_worker_exports <- function(fun, export_env) {
  ensure_pkgs("codetools", auto_install = TRUE)

  bundle_env <- new.env(parent = emptyenv())

  add_binding <- function(name, value) {
    if (!exists(name, envir = bundle_env, inherits = FALSE)) {
      assign(name, value, envir = bundle_env)
    }
  }

  local_names <- ls(envir = export_env, all.names = TRUE)
  for (nm in local_names) {
    add_binding(nm, get(nm, envir = export_env, inherits = FALSE))
  }

  add_binding(".jb_fun_main", fun)
  queue <- c(
    local_names[vapply(local_names, function(nm) {
      is.function(get(nm, envir = export_env, inherits = FALSE))
    }, logical(1))],
    ".jb_fun_main"
  )
  seen <- character(0)

  while (length(queue) > 0L) {
    nm <- queue[[1]]
    queue <- queue[-1]
    if (nm %in% seen || !exists(nm, envir = bundle_env, inherits = FALSE)) next

    obj <- get(nm, envir = bundle_env, inherits = FALSE)
    seen <- c(seen, nm)
    if (!is.function(obj)) next

    fg <- codetools::findGlobals(obj, merge = FALSE)
    refs <- unique(c(fg$functions, fg$variables))
    refs <- setdiff(refs, c("...", ".Generic", ".Method"))

    for (ref in refs) {
      if (exists(ref, envir = bundle_env, inherits = FALSE)) next

      if (exists(ref, envir = export_env, inherits = FALSE)) {
        add_binding(ref, get(ref, envir = export_env, inherits = FALSE))
      } else if (exists(ref, envir = globalenv(), inherits = FALSE)) {
        add_binding(ref, get(ref, envir = globalenv(), inherits = FALSE))
      } else {
        next
      }

      if (is.function(get(ref, envir = bundle_env, inherits = FALSE))) {
        queue <- c(queue, ref)
      }
    }
  }

  list(
    env = bundle_env,
    vars = setdiff(ls(envir = bundle_env, all.names = TRUE), ".jb_fun_main")
  )
}

# Parallel apply that works on Windows/macOS/Linux using PSOCK clusters.
# By default, exports objects from environment(fun) plus recursively referenced
# globals needed by nested helper closures. You can override via export_vars /
# export_env when a script needs stricter control.
par_lapply_safe <- function(X, fun, n_cores = 1L, seed = NULL,
                            export_env = environment(fun),
                            export_vars = NULL,
                            packages = NULL, ...) {
  if (n_cores <= 1L) {
    if (!is.null(seed)) {
      RNGkind("L'Ecuyer-CMRG")
      set.seed(seed)
    }
    return(lapply(X, fun, ...))
  }
  ensure_pkgs(c("parallel"))
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # Reproducible RNG streams
  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }

  # Load packages on workers (optional)
  if (!is.null(packages) && length(packages) > 0) {
    parallel::clusterCall(cl, function(pkgs) {
      for (p in pkgs) {
        suppressMessages(require(p, character.only = TRUE))
      }
      NULL
    }, packages)
  }

  # Export needed objects to workers
  export_source_env <- export_env
  if (is.null(export_vars)) {
    bundle <- build_worker_exports(fun = fun, export_env = export_env)
    export_source_env <- bundle$env
    export_vars <- bundle$vars
  }
  # Preserve clusterSetRNGStream: a caller's global RNG state must never
  # overwrite the distinct worker streams during automatic export.
  export_vars <- setdiff(export_vars, ".Random.seed")
  parallel::clusterExport(cl, varlist = export_vars, envir = export_source_env)

  parallel::parLapply(cl, X, fun, ...)
}



# ============================================================
# Project-root helpers (avoid fragile relative paths)
# ============================================================
jbes_set_root <- function(root) {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  options(JBES_ROOT = root)
  invisible(root)
}

jbes_root <- function() {
  getOption("JBES_ROOT", default = normalizePath(getwd(), winslash = "/", mustWork = FALSE))
}

jbes_path <- function(...) {
  file.path(jbes_root(), ...)
}

# ============================================================
# Environment-variable helpers (useful for smoke tests / overrides)
# ============================================================
env_string <- function(name, default = "") {
  val <- Sys.getenv(name, unset = "")
  if (!nzchar(val)) return(default)
  val
}

env_int <- function(name, default) {
  val <- env_string(name, "")
  if (!nzchar(val)) return(as.integer(default))
  as.integer(val)
}

env_num <- function(name, default) {
  val <- env_string(name, "")
  if (!nzchar(val)) return(as.numeric(default))
  as.numeric(val)
}

env_chr_vec <- function(name, default, sep = ",") {
  val <- env_string(name, "")
  if (!nzchar(val)) return(default)
  trimws(strsplit(val, split = sep, fixed = TRUE)[[1]])
}

env_num_vec <- function(name, default, sep = ",") {
  val <- env_string(name, "")
  if (!nzchar(val)) return(default)
  as.numeric(trimws(strsplit(val, split = sep, fixed = TRUE)[[1]]))
}
