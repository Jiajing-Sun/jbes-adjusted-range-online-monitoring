# Public API and independent-matrix checks; these are not power experiments.
root<-normalizePath(getwd());Sys.setenv(JBES_PROJECT_ROOT=root)
for(f in c('00_utils.R','02_methods.R','03_affine_extrema.R','05_affine_stopping.R','06_training_geometry.R'))source(file.path(root,f))
jbes_set_root(root);set.seed(20261014);checks<-0L
for(q in c(1L,2L,5L))for(m in c(40L,100L))for(gamma in c(0,.15)) {
 X<-matrix(rnorm(3*m*q),3*m,q);psi<-sweep(X,2,colMeans(X[1:m,,drop=FALSE]))
 fit<-aersn::aersn(rep(0,q),psi[1:m,,drop=FALSE]);G<-fit$path$G
 Z<-matrix(apply(psi[m+seq_len(2*m),,drop=FALSE],2,cumsum),ncol=q)/sqrt(m)
 k<-seq_len(2*m);w<-(1+k/m)^2*(k/(m+k))^(2*gamma);ix<-unique(round(seq(1,2*m,length.out=9)))
 ref<-as.numeric(aersn::aersn_gauge(fit$hull,Z[ix,,drop=FALSE]))^2/w[ix]
 stopifnot(max(abs(rsms_path(psi,m,2,gamma)[ix]-ref)/(1+abs(ref)))<1e-7)
 D<-crossprod(G)/m;S<-rowSums((Z%*%inv_pd(D,ridge=1e-8))*Z)/w
 stopifnot(max(abs(ssms_path(psi,m,2,gamma)-S)/(1+abs(S)))<1e-10)
 tr<-scale(psi[1:m,,drop=FALSE],TRUE,FALSE);L<-ceiling(1.06*m^.2);V<-crossprod(tr)/m
 for(l in 1:L){ac<-crossprod(tr[(l+1):m,,drop=FALSE],tr[1:(m-l),,drop=FALSE])/m;V<-V+(1-l/(L+1))*(ac+t(ac))}
 C<-rowSums((Z%*%inv_pd(V,ridge=1e-8))*Z)/w
 stopifnot(max(abs(csms_path(psi,m,2,gamma)-C)/(1+abs(C)))<1e-10)
 if(q==5L&&m==40L&&gamma==0){p<-rsms_path(psi,m,2,gamma);stopifnot(abs(rsms_maximum(psi,m,2,gamma)-max(p))<1e-7)
  for(cv in as.numeric(quantile(p,c(.25,.5,.9))))stopifnot(rsms_first_crossing(psi,m,2,gamma,cv)==first_exceed(p,cv))}
 A<-diag(q)+matrix(rnorm(q*q,sd=.1),q);ga<-aersn_hull_from_G(G%*%A)
 stopifnot(max(abs(aersn_gauge_evaluator(ga)(Z[ix,,drop=FALSE]%*%A)-sqrt(ref*w[ix])))<1e-7)
 checks<-checks+6L
}
x<-matrix(rnorm(100),50,2);x<-scale(x,TRUE,FALSE);V<-crossprod(x)/50
for(l in 1:49){a<-crossprod(x[(l+1):50,,drop=FALSE],x[1:(50-l),,drop=FALSE])/50;V<-V+hac_kernel_weight(l/5,'qs')*(a+t(a))}
stopifnot(max(abs(hac_lrv(x,L=4,kernel='qs')-V))<1e-10)
nested<-function(){m<-40;par_lapply_safe(1:4,function(i)hac_lrv(matrix(rnorm(m)))[1,1],n_cores=2,seed=37)}
stopifnot(all(unlist(nested())>0))
stopifnot(inherits(try(rsms_path(matrix(0,60,1),20,2,0),silent=TRUE),'try-error'))
x<-matrix(rnorm(180),90,2);ob<-persistent_training_objects(x,30);z<-matrix(rnorm(20),10,2)
stopifnot(max(abs(persistent_base_statistics(ob,z)[,1]-aersn::aersn_gauge(aersn_hull_from_G(ob$G),z)^2))<1e-7)
cat('Passed',checks,'path/normalizer checks, affine transformations, exact extrema, first crossings, QS, PSOCK and failure handling.\n')
