# Paired comparison using the original mean-null designs and RNG streams.
source('00_utils.R');source('01_critical_values.R');source('02_methods.R')
source('10_dgps_mean.R');source('30_construct_psi.R');jbes_set_root(getwd())
mg<-read.csv('results/affine_matched_grid_critical_values.csv')
ref<-read.csv('results/affine_finite_critical_values.csv')
N<-1000L;cnt<-0L;rows<-list()
for(dgp in c('DGP1','DGP2','DGP3','DGP4'))for(m in c(100,250,500))for(T in c(1,2,5,10))for(gamma in c(0,.15)) {
 cnt<-cnt+1L;mT<-as.integer(m*T)
 msg('Matched-grid size comparison: ',dgp,', m=',m,', T=',T,', gamma=',gamma)
 one<-function(i){
  X<-sim_var1_dgp(dgp=dgp,n=m+mT,m_train=m,q=2,burnin=250)
  psi<-make_mean_psi(X,m,T)$psi
  max(rsms_path(psi,m,T,gamma))
 }
 vals<-unlist(par_lapply_safe(1:N,one,n_cores=6L,seed=12345+1000L*cnt))
 for(a in c(.05,.10)){
  c400<-ref$critical_value[ref$method=='RSMS'&ref$d==2&ref$T==T&ref$gamma==gamma&ref$alpha==a]
  cmg<-mg$critical_value[mg$method=='RSMS'&mg$d==2&mg$T==T&mg$gamma==gamma&mg$alpha==a&mg$m==m]
  stopifnot(length(c400)==1,length(cmg)==1)
  r400<-vals>c400;rmg<-vals>cmg
  rows[[length(rows)+1L]]<-data.frame(dgp=dgp,m=m,T=T,gamma=gamma,alpha=a,N=N,
   critical_grid400=c400,critical_matched=cmg,size_grid400=mean(r400),size_matched=mean(rmg),
   paired_change=mean(rmg-r400),paired_se=sd(as.numeric(rmg)-as.numeric(r400))/sqrt(N),
   absolute_error_grid400=abs(mean(r400)-a),absolute_error_matched=abs(mean(rmg)-a))
 }
 write.csv(do.call(rbind,rows),'results/matched_grid_size_comparison.csv',row.names=FALSE)
}
out<-do.call(rbind,rows)
old<-read.csv('results/mean_null_size.csv');old<-subset(old,method=='RSMS')
chk<-merge(old,out,by=c('dgp','m','T','gamma','alpha'))
stopifnot(nrow(chk)==192L,max(abs(chk$size-chk$size_grid400))<1e-12)
cat('Paired comparison complete; original 400-grid rejection rates reproduced exactly.\n')
