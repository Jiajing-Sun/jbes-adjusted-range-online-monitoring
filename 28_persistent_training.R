# Persistent training break: fixed design, independent calibration/evaluation,
# common random numbers across B, gamma, directions and competing methods.
# Run from replication/. The original data generator and all old results are retained.
source('00_utils.R'); source('01_critical_values.R')
source('02_methods.R'); source('10_dgps_mean.R'); source('06_training_geometry.R')
jbes_set_root(normalizePath(getwd()))
N <- env_int('JBES_PERSIST_N',1000L)
CORES <- env_int('JBES_N_CORES',8L)
MS <- c(250L,500L,1000L); BS <- c(0,.25,.5,.75,1)
GS <- c(0,.15); METHODS <- c('RSMS','SSMS','CSMS')
MODES <- c('clean','oracle'); TT <- 2; PI <- .2; DELTA <- .5
SCENES <- data.frame(name=c('none','parallel_early','parallel_late',
  'reverse_early','reverse_late','transverse_early','transverse_late'),
  bf=c(0,.1,.8,.1,.8,.1,.8),dx=c(0,.5,.5,-.5,-.5,.5,.5),
  dy=c(0,.5,.5,-.5,-.5,-.5,-.5))
SEED <- 91427000L
OUT <- Sys.getenv('JBES_PERSIST_OUT','results/persistent_training')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
clean <- array(NA_real_,c(2,3))
for(g in 1:2)for(j in 1:3) clean[g,j] <- switch(METHODS[j],
 RSMS=get_critical_value_rsms(2,TT,GS[g],.05),
 SSMS=get_critical_value_ssms(2,TT,GS[g],.05),
 CSMS=get_critical_value_csms(2,TT,GS[g],.05))
summary_rows <- list(); cv_rows <- list(); shape_rows <- list(); pair_rows <- list()
for(mi in seq_along(MS)) {
 m <- MS[mi]; nmon <- m*TT; k <- seq_len(nmon)
 weights <- vapply(GS,function(g)1/((1+k/m)^2*(k/(m+k))^(2*g)),numeric(nmon))
 oracle <- array(NA_real_,c(5,2,3))
 for(phase in 1:3) {
   # phase 1: calibration no new break; phase 2: independent no-new-break check;
   # phase 3: independent alternatives (all six scenarios, with paired noises).
   scenarios <- if(phase<3)1L else 2:7
   one <- function(i) {
    RNGkind("L'Ecuyer-CMRG");set.seed(SEED+mi*1000000L+phase*100000L+i)
    X <- sim_var1_dgp('DGP2',m+nmon,m,2,burnin=300)
    mx <- array(NA_real_,c(5,2,3)); stops <- array(nmon+1L,c(5,2,length(scenarios),3,2))
    kap <- numeric(5)
    for(bi in 1:5) {
     Y <- X; ix <- (m-floor(PI*m)+1L):(m+nmon)
     Y[ix,] <- Y[ix,,drop=FALSE]+BS[bi]
     ob <- persistent_training_objects(Y,m);kap[bi] <- ob$kappa
     Z0 <- apply(ob$psi[m+seq_len(nmon),,drop=FALSE],2,cumsum)/sqrt(m)
     for(si in seq_along(scenarios)) {
      sc <- SCENES[scenarios[si],]; kk <- floor(sc$bf*nmon)
      Z <- Z0 + (pmax(k-kk,0)/sqrt(m)) %o% c(sc$dx,sc$dy)
      raw <- persistent_base_statistics(ob,Z)
      for(gi in 1:2) {
       paths <- raw*weights[,gi]
       if(phase==1)mx[bi,gi,] <- apply(paths,2,max) else
        for(j in 1:3)for(ca in 1:2) {
         cv <- if(ca==1)clean[gi,j] else oracle[bi,gi,j]
         stops[bi,gi,si,j,ca] <- first_exceed(paths[,j],cv)
        }
      }
     }
    }
    list(maxima=if(phase==1)mx else NULL,stops=if(phase>1)stops else NULL,kappa=kap)
   }
   file <- file.path(OUT,sprintf('m%d_phase%d_N%d.rds',m,phase,N))
   if(file.exists(file)) {
    ans <- readRDS(file);stopifnot(length(ans)==N)
   } else {
    cat('Running m=',m,' phase=',phase,' N=',N,'\n');flush.console()
    ans <- par_lapply_safe(seq_len(N),one,CORES,seed=NULL)
    saveRDS(ans,file,compress=FALSE)
   }
   if(phase==1) {
    for(bi in 1:5)for(gi in 1:2)for(j in 1:3) {
     oracle[bi,gi,j] <- quantile(vapply(ans,function(z)z$maxima[bi,gi,j],0),.95,type=8)
     cv_rows[[length(cv_rows)+1L]] <- data.frame(m=m,B=BS[bi],gamma=GS[gi],method=METHODS[j],
       clean_cv=clean[gi,j],oracle_cv=oracle[bi,gi,j],N_CAL=N)
    }
   } else {
    for(bi in 1:5)for(gi in 1:2)for(ca in 1:2) {
     ratio <- if(ca==1)clean[gi,2]/clean[gi,1] else oracle[bi,gi,2]/oracle[bi,gi,1]
     kap <- vapply(ans,function(z)z$kappa[bi],0); event <- kap<=ratio
     if(phase==2)shape_rows[[length(shape_rows)+1L]] <- data.frame(m=m,B=BS[bi],gamma=GS[gi],calibration=MODES[ca],
       ratio=ratio,containment_rate=mean(event),se=sqrt(mean(event)*(1-mean(event))/N),
       kappa_median=median(kap),kappa_q95=quantile(kap,.95),N=N)
     for(si in seq_along(scenarios)) {
      sc <- SCENES[scenarios[si],]; ks <- floor(sc$bf*nmon)
      st <- sapply(1:3,function(j)vapply(ans,function(z)z$stops[bi,gi,si,j,ca],0))
      for(j in 1:3) {
       alarm <- st[,j]<=nmon; post <- st[,j]>ks & alarm; delay <- st[post,j]-ks
       summary_rows[[length(summary_rows)+1L]] <- data.frame(m=m,B=BS[bi],pi=PI,T=TT,gamma=GS[gi],
         scenario=sc$name,break_frac=sc$bf,dx=sc$dx,dy=sc$dy,method=METHODS[j],calibration=MODES[ca],
         alarm=mean(alarm),pre=if(phase==2)NA_real_ else mean(st[,j]<=ks),
         post=if(phase==2)NA_real_ else mean(post),alarm_se=sqrt(mean(alarm)*(1-mean(alarm))/N),
         post_se=if(phase==2)NA_real_ else sqrt(mean(post)*(1-mean(post))/N),
         CADD=if(phase==2||!length(delay))NA_real_ else mean(delay),
         CADD_se=if(phase==2||length(delay)<2)NA_real_ else sd(delay)/sqrt(length(delay)),
         post_count=if(phase==2)NA_integer_ else sum(post),N=N)
      }
      if(phase==3) {
       pr <- st[,1]>ks & st[,1]<=nmon; ps <- st[,2]>ks & st[,2]<=nmon
       gap <- as.integer(pr)-as.integer(ps)
       paired <- data.frame(m=m,B=BS[bi],gamma=GS[gi],scenario=sc$name,calibration=MODES[ca],
         post_gap_R_minus_S=mean(gap),paired_se=sd(gap)/sqrt(N),
         containment_rate=mean(event),R_pre=mean(st[,1]<=ks),
         bound_lower=-mean(!event)-mean(st[,1]<=ks),
         containment_violations=sum(event & st[,1]>st[,2]),N=N)
       common <- event & st[,1]>ks
       paired$common_survival_count <- sum(common)
       paired$common_survival_probability <- mean(common)
       paired$R_post_given_common <- if(any(common))mean(pr[common]) else NA_real_
       paired$S_post_given_common <- if(any(common))mean(ps[common]) else NA_real_
       paired$R_earlier_given_common <- if(any(common))mean(st[common,1]<st[common,2]) else NA_real_
       paired$R_pre_count <- sum(st[,1]<=ks)
       paired$S_pre_count <- sum(st[,2]<=ks)
       paired$R_post_count <- sum(pr)
       paired$S_post_count <- sum(ps)
       paired$R_no_crossing_count <- sum(st[,1]>nmon)
       paired$S_no_crossing_count <- sum(st[,2]>nmon)
       # The theory uses unregularized SSMS; the retained 1e-8 ridge only reduces its path.
       stopifnot(paired$containment_violations==0L,mean(gap)+1e-12>=paired$bound_lower)
       pair_rows[[length(pair_rows)+1L]] <- paired
      }
     }
    }
   }
   cat('Completed m=',m,' phase=',phase,'\n');flush.console()
 }
}
write.csv(do.call(rbind,summary_rows),file.path(OUT,'persistent_summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,cv_rows),file.path(OUT,'persistent_critical_values.csv'),row.names=FALSE)
write.csv(do.call(rbind,shape_rows),file.path(OUT,'persistent_geometry.csv'),row.names=FALSE)
write.csv(do.call(rbind,pair_rows),file.path(OUT,'persistent_paired_comparison.csv'),row.names=FALSE)
cat('All prespecified persistent-training results saved.\n')
