args <- commandArgs(FALSE)
file <- sub('^--file=', '', args[grepl('^--file=',args)][1])
ROOT <- normalizePath(Sys.getenv('JBES_PROJECT_ROOT',dirname(normalizePath(file))))
OUT <- file.path(ROOT,'results','automatic_bandwidth');dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
source(file.path(ROOT,'00_utils.R')); jbes_set_root(ROOT)
for(f in c('01_critical_values.R','02_methods.R','10_dgps_mean.R','30_construct_psi.R')) source(file.path(ROOT,f))
stopifnot(requireNamespace('sandwich',quietly=TRUE))
N <- as.integer(Sys.getenv('JBES_N_MC','1000')); cores <- as.integer(Sys.getenv('JBES_N_CORES','4'))
m <- 250L; gamma <- .15
spec <- expand.grid(kernel=c('Bartlett','Parzen','Quadratic Spectral'),selector=c('Andrews AR1','Newey-West'),stringsAsFactors=FALSE)
names_method <- c('RSMS','CSMS default',paste(spec$kernel,spec$selector,sep=' / '))
# Full-lag kernel sum; in particular QS is not truncated at the bandwidth.
lrv_from_acf <- function(ac,bw,kernel) {
 w <- sandwich::kweights((0:(m-1))/bw,kernel=kernel)
 A <- apply(ac, c(2,3),function(x) sum(w*x))
 A+t(A)-ac[1,,]
}
set.seed(7101); u <- matrix(rnorm(m*2),m,2);u <- scale(u,TRUE,FALSE)
ac <- acf(u,lag.max=m-1,type='covariance',plot=FALSE,demean=FALSE)$acf
for(k in unique(spec$kernel)) {
 bw <- 7.4; a <- lrv_from_acf(ac,bw,k); b <- crossprod(u)/m
 for(l in 1:(m-1)){v <- crossprod(u[(l+1):m,,drop=FALSE],u[1:(m-l),,drop=FALSE])/m;b <- b+sandwich::kweights(l/bw,kernel=k)*(v+t(v))}
 stopifnot(max(abs(a-b))<1e-12)
 ref <- aersn::aersn_hac_lrv(aersn_training_fit(u),kernel=k,bandwidth=bw,prewhite=FALSE,adjust=FALSE,center=TRUE)$matrix
 stopifnot(max(abs(a-ref))<1e-10)
}
stopifnot(max(abs(lrv_from_acf(ac,5,'Bartlett')-hac_lrv(u,L=4)))<1e-12)
writeLines(c('Verified FFT covariance sum against direct all-lag sums for all three kernels.',paste('sandwich version:',packageVersion('sandwich'))),file.path(OUT,'implementation_check.txt'))
run_rep <- function(i,T,phase,cv=NULL) {
 set.seed(20260914L + T*100000L + phase*10000L+i)
 h <- m*T; X <- sim_var1_dgp('DGP1',m+h,m,2,burnin=250)
 kstar <- 0L
 if(phase %in% c(3,4)) {
  kstar <- floor(if(phase==3) .1*h else .8*h)
  X <- apply_mean_change(X,m,T,kstar,if(phase==3) .5 else 1.5,rep(1,2),'abrupt')
 }
 psi <- make_mean_psi(X,m,T)$psi
 train <- scale(psi[1:m,,drop=FALSE],TRUE,FALSE)
 ac <- acf(train,lag.max=m-1,type='covariance',plot=FALSE,demean=FALSE)$acf
 bws <- numeric(nrow(spec)); gammas <- list()
 for(j in seq_len(nrow(spec))) {
  fn <- if(spec$selector[j]=='Andrews AR1') sandwich::bwAndrews else sandwich::bwNeweyWest
  bws[j] <- fn(train,kernel=spec$kernel[j],prewhite=0,weights=rep(1,2))
  stopifnot(is.finite(bws[j]),bws[j]>0)
  gammas[[j]] <- aersn::aersn_hac_lrv(aersn_training_fit(train),kernel=spec$kernel[j],
    bandwidth=bws[j],prewhite=FALSE,adjust=FALSE,center=TRUE)$matrix
 }
 S <- apply(psi[m+seq_len(h),,drop=FALSE],2,cumsum)
 k <- seq_len(h);sc <- 1/m/(1+k/m)^2/(k/(m+k))^(2*gamma)
 paths <- cbind(rsms_path(psi,m,T,gamma),csms_path(psi,m,T,gamma),sapply(gammas,function(g) rowSums((S%*%inv_pd(g,ridge=1e-8))*S)*sc))
 stopifnot(all(is.finite(paths)))
 if(is.null(cv)) return(list(maxima=apply(paths,2,max),bw=bws))
 stops <- sapply(seq_len(ncol(cv)),function(j) vapply(seq_len(ncol(paths)),function(k)first_exceed(paths[,k],cv[k,j]),numeric(1)))
 list(stops=stops,bw=bws)
}
rows <- list();bwrows <- list();cvrows <- list()
for(T in c(1L,5L)) {
 cat('T',T,'calibration',N,'\n');flush.console()
 cal <- parallel::mclapply(seq_len(N),run_rep,T=T,phase=1,mc.cores=cores,mc.set.seed=FALSE)
 stopifnot(!any(vapply(cal,inherits,logical(1),'try-error')))
 maxima <- do.call(rbind,lapply(cal,`[[`,'maxima'))
 chat <- apply(maxima,2,quantile,probs=.95,type=8)
 raw <- c(get_critical_value_rsms(2,T,gamma,.05),rep(get_critical_value_csms(2,T,gamma,.05),7))
 cv <- cbind(asymptotic=raw,independent_calibration=chat)
 cvrows[[length(cvrows)+1]] <- data.frame(T,method=names_method,raw,chat)
 saveRDS(cal,file.path(OUT,paste0('calibration_T',T,'.rds')))
 for(phase in 2:4) {
  cat('T',T,'phase',phase,N,'\n');flush.console()
  z <- parallel::mclapply(seq_len(N),run_rep,T=T,phase=phase,cv=cv,mc.cores=cores,mc.set.seed=FALSE)
  stopifnot(!any(vapply(z,inherits,logical(1),'try-error')))
  saveRDS(z,file.path(OUT,paste0('evaluation_T',T,'_phase',phase,'.rds')))
  bw <- do.call(rbind,lapply(z,`[[`,'bw'))
  for(j in 1:6) bwrows[[length(bwrows)+1]] <- data.frame(T,phase,method=names_method[j+2],mean=mean(bw[,j]),median=median(bw[,j]),q05=quantile(bw[,j],.05),q95=quantile(bw[,j],.95))
  for(b in 1:2) for(j in 1:8) {
   stops <- vapply(z,function(x)x$stops[j,b],numeric(1)); h <- m*T
   kstar <- if(phase==2) 0 else floor(if(phase==3) .1*h else .8*h)
   met <- detection_metrics(stops,kstar,h)
   rows[[length(rows)+1]] <- data.frame(T,phase=c('null','early','late')[phase-1],boundary=colnames(cv)[b],method=names_method[j],N,DP_horizon=met['DP_horizon'],PFA_pre=met['PFA_pre'],DP_post=met['DP_post'],CADD=met['CADD'],MCSE=sqrt(met['DP_post']*(1-met['DP_post'])/N))
  }
  write.csv(do.call(rbind,rows),file.path(OUT,'selector_results.csv'),row.names=FALSE)
 }
}
write.csv(do.call(rbind,bwrows),file.path(OUT,'selected_bandwidths.csv'),row.names=FALSE)
write.csv(do.call(rbind,cvrows),file.path(OUT,'selector_critical_values.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(OUT,'sessionInfo.txt'))
cat('DONE\n')

for(pair in list(c("selector_results.csv","csms_automatic_bandwidth_results.csv"),c("selected_bandwidths.csv","csms_automatic_selected_bandwidths.csv"),c("selector_critical_values.csv","csms_automatic_critical_values.csv"))) file.copy(file.path(OUT,pair[1]),file.path(ROOT,"results",pair[2]),overwrite=TRUE)
