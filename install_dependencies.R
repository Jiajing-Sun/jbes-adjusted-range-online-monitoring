# Install public CRAN dependencies; no data or simulation output is downloaded.
packages <- c('MASS','Matrix','data.table','dplyr','tidyr','readr','lubridate',
 'ggplot2','patchwork','showtext','fda','xtable','codetools','jsonlite','lpSolve','sandwich')
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))install.packages(missing,repos='https://cloud.r-project.org')
# Exact versions define the recorded numerical implementation.
for(spec in list(c('aersn','0.2.3'),c('highs','1.14.0.2'))) {
 name<-spec[1];version<-spec[2]
 if(!requireNamespace(name,quietly=TRUE)||as.character(packageVersion(name))!=version) {
  if(!requireNamespace('remotes',quietly=TRUE))install.packages('remotes',repos='https://cloud.r-project.org')
  remotes::install_version(name,version=version,repos='https://cloud.r-project.org',upgrade='never')
 }
}
cat('Dependencies installed. Run Rscript verify_backend.R before production.\n')
