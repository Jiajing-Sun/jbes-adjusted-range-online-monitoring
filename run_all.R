# Compatibility entry point for the documented production pipeline.
# Use reproduce_reported.R --list to inspect it without running simulations.
root <- normalizePath(Sys.getenv("JBES_PROJECT_ROOT", unset=getwd()))
status <- system2(file.path(R.home("bin"), "Rscript"),
  c(shQuote(file.path(root, "reproduce_reported.R")), "--run"))
if(status != 0L) stop("The production pipeline failed; inspect logs/.")
