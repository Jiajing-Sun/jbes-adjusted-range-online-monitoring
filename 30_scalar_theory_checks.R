# Direct checks of the manuscript's scalar persistent-training propositions.
# Independent calibration, null evaluation and alternative evaluation streams.
# The variance-profile experiment is deliberately outside this script's scope.
root<-normalizePath(Sys.getenv('JBES_PROJECT_ROOT',getwd()))
for(f in c('00_utils.R','01_critical_values.R','02_methods.R'))source(file.path(root,f))
jbes_set_root(root)
out<-file.path(root,Sys.getenv('JBES_THEORY_OUT','results/scalar_theory'));dir.create(out,recursive=TRUE,showWarnings=FALSE)
N_CAL<-env_int('JBES_THEORY_CAL',10000L);N_EVAL<-env_int('JBES_THEORY_EVAL',5000L)
CORES<-env_int('JBES_N_CORES',4L);MS<-env_num_vec('JBES_THEORY_M',c(250,1000,4000))
PIS<-c(.25,.5);NOISES<-c('iid_gaussian','ar_t5');T0<-2;TAU<-.4;ALPHA<-.05;B<-1
H<-c(-4,-2,0,1,2,4);METHODS<-c('RSMS','SSMS','CSMS')
# A stationary AR(1) with standardized t5 innovations and long-run variance one.
noise_sample<-function(n,kind) {
 if(kind=='iid_gaussian')return(rnorm(n))
 a<-.4;e<-rt(n+500,5)/sqrt(5/3)*(1-a)
 as.numeric(stats::filter(e,filter=a,method='recursive'))[501:(n+500)]
}
select_scalar_training<-function(x,dependent) {
 m<-length(x);r<-(1:m)/m;lo<-ceiling(.1*m);hi<-floor(.9*m);idx<-lo:hi
 cpath<-cumsum(x-mean(x))/m;Q<-cpath[idx]^2/(r[idx]*(1-r[idx]))
 j<-idx[which.max(Q)];selected<-max(Q)/mean((x-mean(x))^2)>m^(-.5)
 if(!selected)return(list(selected=FALSE,split=NA_integer_,pi=NA_real_,b=NA_real_,sigma=NA_real_))
 bhat<-mean(x[(j+1):m])-mean(x[1:j]);res<-c(x[1:j]-mean(x[1:j]),x[(j+1):m]-mean(x[(j+1):m]))
 v<-if(dependent)hac_lrv(matrix(res),L=floor(m^.2))[1,1] else mean(res^2)
 list(selected=TRUE,split=j,pi=1-j/m,b=bhat,sigma=sqrt(v))
}
# Exact first crossings on the full calendar grid; finite-horizon noncrossing=N+1.
scalar_record<-function(i,phase,m,pi0,noise,oracle=NULL,clean=NULL) {
 config<-(match(noise,NOISES)-1L)*length(PIS)*length(MS)+(match(pi0,PIS)-1L)*length(MS)+match(m,MS)
 set.seed(20261014L+config*100000L+phase*20000L+i)
 N<-as.integer(m*T0);k<-seq_len(N);ks<-floor(m*TAU);rho<-1-pi0
 eps<-noise_sample(m+N,noise);x<-eps+B*c(rep(0,floor(m*rho)),rep(1,m+N-floor(m*rho)))
 theta<-mean(x[1:m]);train<-matrix(x[1:m]-theta)
 fit<-aersn_training_fit(train)
 den<-c(diff(range(fit$path$G))^2,
        aersn::aersn_shao_normalizer(fit,integration='calendar')$matrix[1,1],
        hac_lrv(train)[1,1])+c(0,1e-8,1e-8)
 base<-cumsum(x[m+k]-theta)/sqrt(m);w<-(1+k/m)^2
 paths<-function(shift) outer((base+pmax(k-ks,0)*shift/sqrt(m))^2/w,1/den)
 if(phase==1)return(apply(paths(0),2,max))
 est<-select_scalar_training(x[1:m],noise!='iid_gaussian')
 fitted<-clean
 if(est$selected) {
  pp<-est$pi;rr<-1-pp;beta<-est$sigma/abs(est$b);ell<-T0^2/(pp^2*(1+T0)^2)
  cc<-2*T0/(pp^2*rr*(1+T0)^2);vv<-c(T0+T0^2/pp,T0+T0^2*(3*pp^2+2*pp+1)/(5*pp^2))
  fitted[1:2]<-c(1,3)*(ell+cc*beta*sqrt(vv/m)*qnorm(1-ALPHA))
 }
 # Only R/S have the stated fitted scalar boundary. CSMS fitted entries are NA.
 fitted[3]<-NA_real_;cv<-cbind(clean=clean,oracle=oracle,fitted=fitted)
 hh<-if(phase==2)0 else H
 stops<-array(N+1L,c(length(hh)+as.integer(phase==3),3,3))
 for(j in seq_along(hh)) {
  pp<-paths(B*hh[j]/sqrt(m))
  for(a in 1:3)for(ca in 1:3)stops[j,a,ca]<-if(is.finite(cv[a,ca]))first_exceed(pp[,a],cv[a,ca]) else NA_real_
 }
 if(phase==3) {
  pp<-paths(.75*B)
  for(a in 1:3)for(ca in 1:3)stops[length(hh)+1,a,ca]<-if(is.finite(cv[a,ca]))first_exceed(pp[,a],cv[a,ca]) else NA_real_
 }
 list(stops=stops,selected=est$selected,split_error=est$split-floor(m*rho))
}
summary_rows<-pair_rows<-cal_rows<-select_rows<-curve_rows<-list()
for(noise in NOISES)for(pi0 in PIS)for(m in MS) {
 tag<-sprintf('%s_pi%.2f_m%d',noise,pi0,m);cat('Scalar theory:',tag,'\n');flush.console()
 clean<-c(get_critical_value_rsms(1,T0,0,ALPHA),get_critical_value_ssms(1,T0,0,ALPHA),get_critical_value_csms(1,T0,0,ALPHA))
 cal<-par_lapply_safe(seq_len(N_CAL),function(i)scalar_record(i,1,m,pi0,noise),CORES)
 maxima<-do.call(rbind,cal);oracle<-apply(maxima,2,quantile,1-ALPHA,type=8)
 saveRDS(maxima,file.path(out,paste0(tag,'_calibration.rds')))
 # Quantile brackets isolate calibration uncertainty from evaluation uncertainty.
 prob<-pmin(1,pmax(0,(1-ALPHA)+c(-1,1)*qnorm(.975)*sqrt(ALPHA*(1-ALPHA)/N_CAL)))
 for(j in 1:3)cal_rows[[length(cal_rows)+1]]<-data.frame(noise,pi=pi0,m,method=METHODS[j],critical=oracle[j],lo=quantile(maxima[,j],prob[1]),hi=quantile(maxima[,j],prob[2]),n_cal=N_CAL)
 for(phase in 2:3) {
  ans<-par_lapply_safe(seq_len(N_EVAL),function(i)scalar_record(i,phase,m,pi0,noise,oracle,clean),CORES)
  saveRDS(ans,file.path(out,paste0(tag,'_phase',phase,'.rds')))
  select_rows[[length(select_rows)+1]]<-data.frame(noise,pi=pi0,m,phase,selection=mean(vapply(ans,`[[`,FALSE,'selected')),mean_split_error=mean(vapply(ans,`[[`,0,'split_error'),na.rm=TRUE))
  scenes<-if(phase==2)'null' else c(paste0('local_h',H),'fixed_0.75')
  ks<-floor(m*TAU);N<-m*T0;v<-c(T0+T0^2/pi0,T0+T0^2*(3*pi0^2+2*pi0+1)/(5*pi0^2))
  for(si in seq_along(scenes))for(ca in 1:3) {
   st<-vapply(ans,function(z)z$stops[si,,ca],numeric(3));st<-t(st)
   for(j in 1:3) {
    if(all(is.na(st[,j])))next
    pre<-st[,j]<=ks;post<-st[,j]>ks&st[,j]<=N;d<-st[post,j]-ks
    theo<-if(phase==3&&si<=length(H)&&j<=2)pnorm(qnorm(1-ALPHA)-B*H[si]*(T0-TAU)/sqrt(v[j]),lower.tail=FALSE) else NA_real_
    summary_rows[[length(summary_rows)+1]]<-data.frame(noise,pi=pi0,m,scenario=scenes[si],calibration=colnames(cbind(clean,oracle,fitted=clean))[ca],method=METHODS[j],alarm=mean(st[,j]<=N),pre=mean(pre),post=mean(post),post_count=sum(post),survival_count=sum(!pre),conditional_post=if(any(!pre))mean(post[!pre]) else NA_real_,post_se=sqrt(mean(post)*(1-mean(post))/N_EVAL),cadd=if(length(d))mean(d) else NA_real_,cadd_se=if(length(d)>1)sd(d)/sqrt(length(d)) else NA_real_,theory_power=theo,n_eval=N_EVAL)
   }
   pr<-st[,1]>ks&st[,1]<=N;ps<-st[,2]>ks&st[,2]<=N;gap<-as.integer(pr)-as.integer(ps)
   dr<-ifelse(pr,st[,1]-ks,0);ds<-ifelse(ps,st[,2]-ks,0)
   mr<-if(any(pr))mean(dr[pr]) else NA_real_;ms<-if(any(ps))mean(ds[ps]) else NA_real_
   if(scenes[si]=='fixed_0.75') {
    rho<-1-pi0;hh<-pi0*rho;cc<-rho*T0/(1+T0)
    t0<-(cc+.75*TAU)/(rho+.75-cc);ff<-(rho+.75)*t0-.75*TAU
    derivative<-2*ff*(rho+.75+.75*TAU)/(hh^2*(1+t0)^3)
    delay_limit<-2*rho*T0/(hh^2*(1+T0)^2)*qnorm(1-ALPHA)*(sqrt(v[2])-sqrt(v[1]))/derivative
   } else delay_limit<-NA_real_
   for(uj in c(.1,.2,.4,.6,.8,.9,1)*T0)for(j in 1:3)if(!all(is.na(st[,j])))
    curve_rows[[length(curve_rows)+1]]<-data.frame(noise,pi=pi0,m,scenario=scenes[si],calibration=c('clean','oracle','fitted')[ca],method=METHODS[j],u=uj,first_crossing_probability=mean(st[,j]<=floor(m*uj)),n_eval=N_EVAL)
   # Joint ratio-estimator influence functions preserve separate detection events.
   infl<-if(any(pr)&&any(ps))ps*(ds-ms)/mean(ps)-pr*(dr-mr)/mean(pr) else rep(NA_real_,N_EVAL)
   pair_rows[[length(pair_rows)+1]]<-data.frame(noise,pi=pi0,m,scenario=scenes[si],calibration=c('clean','oracle','fitted')[ca],post_gap_R_minus_S=mean(gap),paired_se=sd(gap)/sqrt(N_EVAL),cadd_gap_S_minus_R=ms-mr,scaled_gap=(ms-mr)/sqrt(m),theory_scaled_gap=delay_limit,cadd_gap_se=if(all(is.finite(infl)))sd(infl)/sqrt(N_EVAL) else NA_real_,n_eval=N_EVAL)
  }
 }
 write.csv(do.call(rbind,summary_rows),file.path(out,'scalar_summary.csv'),row.names=FALSE)
 write.csv(do.call(rbind,pair_rows),file.path(out,'scalar_paired.csv'),row.names=FALSE)
 write.csv(do.call(rbind,cal_rows),file.path(out,'scalar_calibration.csv'),row.names=FALSE)
 write.csv(do.call(rbind,curve_rows),file.path(out,'scalar_first_crossing.csv'),row.names=FALSE)
 write.csv(do.call(rbind,select_rows),file.path(out,'scalar_selection.csv'),row.names=FALSE)
}
writeLines(capture.output(sessionInfo()),file.path(out,'session.txt'))
