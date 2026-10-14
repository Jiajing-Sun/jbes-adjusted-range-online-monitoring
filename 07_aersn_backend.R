# aersn 0.2.3 normalizers for the existing sequential statistics.
# No offline aersn p-value or critical value is used for an online maximum.
stopifnot(requireNamespace('aersn',quietly=TRUE), packageVersion('aersn')=='0.2.3')
# The private environment travels with this closure to PSOCK workers.
aersn_training_fit <- local({
 cache <- new.env(parent=emptyenv())
 function(train) {
  train<-as.matrix(train)
  if(!is.null(cache$train) && identical(train,cache$train))return(cache$fit)
  fit<-aersn::aersn(rep(0,ncol(train)),train)
  cache$train<-train;cache$fit<-fit;fit
 }
})
# Convert an already scaled closed path to public aersn path/hull objects.
# Adding zero endpoints leaves its increment hull unchanged. This also handles
# unordered vertex-only inputs used by certified interval screening.
aersn_hull_from_G <- function(G) {
 G<-as.matrix(G);G<-rbind(rep(0,ncol(G)),G,rep(0,ncol(G)))
 n<-nrow(G)-1L;psi<-(G[-1,,drop=FALSE]-G[-nrow(G),,drop=FALSE])*sqrt(n)
 h<-aersn::aersn_hull(aersn::aersn_path(psi),n_directions=0L)
 if(!h$full_rank)stop('aersn: rank-deficient training hull; path is retained as a failure, not resampled')
 h
}
# Optional reusable LP for the same public aersn hull in q > 2.
# Every new hull is checked against aersn_gauge at its first objective.
aersn_highs_lp_engine <- ae_lp_engine
# q=2: exact facets of the polygon returned by the public package interface.
# Batch evaluation avoids solving the same LP at every monitoring time.
aersn_gauge_evaluator <- function(h) {
 if(h$q==1L) {
  unit<-as.numeric(aersn::aersn_gauge(h,1))
  return(function(Z)as.numeric(abs(as.matrix(Z))*unit))
 }
 if(h$q==2L) {
  V<-aersn::aersn_vertices(h);W<-V[c(2:nrow(V),1),,drop=FALSE];E<-W-V
  normals<-cbind(E[,2],-E[,1]);c0<-rowSums(normals*V)
  if(any(c0<=0)) {normals<- -normals;c0<- -c0}
  stopifnot(all(c0>0));A<-normals/c0
  return(function(Z) {
   Z<-matrix(Z,ncol=2);ans<-rep(0,nrow(Z))
   for(j in seq_len(nrow(A)))ans<-pmax(ans,as.numeric(Z%*%A[j,]))
   ans
  })
 }
 if(identical(Sys.getenv('JBES_AERSN_GAUGE','highs'),'highs') && requireNamespace('highs',quietly=TRUE)) {
  engine<-aersn_highs_lp_engine(h$G);checked<-FALSE
  return(function(Z) {
   Z<-matrix(Z,ncol=h$q)
   out<-vapply(seq_len(nrow(Z)),function(i)engine(Z[i,]),numeric(1))
   if(!checked && nrow(Z)) {
    ref<-as.numeric(aersn::aersn_gauge(h,Z[1,,drop=FALSE]))
    if(abs(out[1]-ref)>1e-7*(1+abs(ref)))stop('HiGHS and public aersn gauge disagree')
    checked<<-TRUE
   }
   out
  })
 }
 function(Z)as.numeric(aersn::aersn_gauge(h,matrix(Z,ncol=h$q)))
}
ae_lp_engine <- function(G) {
 eval<-aersn_gauge_evaluator(aersn_hull_from_G(G))
 function(z)eval(matrix(z,nrow=1))
}
ae_gauge_values <- function(G,Z) aersn_gauge_evaluator(aersn_hull_from_G(G))(as.matrix(Z))
rsms_path <- function(psi_mat,m,T,gamma,ridge=1e-8) {
 psi_mat<-as.matrix(psi_mat);N<-as.integer(round(m*T));q<-ncol(psi_mat)
 stopifnot(m>q,N>0,nrow(psi_mat)>=m+N,gamma>=0,gamma<.5)
 h<-aersn::aersn_hull(aersn::aersn_path(psi_mat[seq_len(m),,drop=FALSE]),n_directions=0)
 if(!h$full_rank)stop('aersn: rank deficient training hull')
 Z<-matrix(apply(psi_mat[m+seq_len(N),,drop=FALSE],2,cumsum),ncol=q)/sqrt(m)
 k<-seq_len(N);aersn_gauge_evaluator(h)(Z)^2/((1+k/m)^2*(k/(m+k))^(2*gamma))
}
ssms_path <- function(psi_mat,m,T,gamma,ridge=1e-8) {
 psi_mat<-as.matrix(psi_mat);N<-as.integer(round(m*T));q<-ncol(psi_mat)
 fit<-aersn_training_fit(psi_mat[seq_len(m),,drop=FALSE])
 V<-aersn::aersn_shao_normalizer(fit,integration='calendar')$matrix
 Z<-matrix(apply(psi_mat[m+seq_len(N),,drop=FALSE],2,cumsum),ncol=q)/sqrt(m)
 k<-seq_len(N);rowSums((Z%*%inv_pd(V,ridge=ridge))*Z)/((1+k/m)^2*(k/(m+k))^(2*gamma))
}
hac_lrv <- function(psi_train,L=NULL,kernel='bartlett') {
 m<-nrow(psi_train);if(is.null(L))L<-ceiling(1.06*m^(1/5))
 L<-max(1L,min(as.integer(L),m-1L));kernel<-normalize_hac_kernel(kernel)
 kernels<-c(bartlett='Bartlett',parzen='Parzen',qs='Quadratic Spectral')
 stopifnot(kernel%in%names(kernels))
 # Package applies all-lag QS weights. The old finite-truncation QS sensitivity
 # is intentionally corrected and must be recomputed, not silently relabeled.
 aersn::aersn_hac_lrv(aersn_training_fit(psi_train),kernel=kernels[[kernel]],bandwidth=L+1,
                    prewhite=FALSE,adjust=FALSE,center=TRUE)$matrix
}
