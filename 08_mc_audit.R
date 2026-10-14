# Event counts and paired uncertainty conditional on the calibrated boundaries.
# Separate detection events are retained when comparing conditional mean delays.
mc_event_audit <- function(stops,k_star,horizon) {
 stops<-as.matrix(stops);stopifnot(all(is.finite(stops)),all(stops>=1))
 n<-nrow(stops);names<-colnames(stops);post<-stops>k_star&stops<=horizon;pre<-stops<=k_star
 per<-do.call(rbind,lapply(seq_len(ncol(stops)),function(j){
  d<-stops[post[,j],j]-k_star;p<-mean(post[,j])
  data.frame(method=names[j],n=n,pre_count=sum(pre[,j]),post_count=sum(post[,j]),no_crossing_count=sum(stops[,j]>horizon),post=p,post_se=sqrt(p*(1-p)/n),survival_count=sum(!pre[,j]),conditional_post=if(any(!pre[,j]))mean(post[!pre[,j],j]) else NA_real_,cadd=if(length(d))mean(d) else NA_real_,cadd_se=if(length(d)>1)sd(d)/sqrt(length(d)) else NA_real_)
 }))
 pairs<-combn(seq_len(ncol(stops)),2,simplify=FALSE)
 paired<-do.call(rbind,lapply(pairs,function(ix){
  a<-ix[1];b<-ix[2];gap<-as.integer(post[,a])-as.integer(post[,b]);pa<-mean(post[,a]);pb<-mean(post[,b])
  da<-ifelse(post[,a],stops[,a]-k_star,0);db<-ifelse(post[,b],stops[,b]-k_star,0)
  ma<-if(pa>0)sum(da)/sum(post[,a]) else NA_real_;mb<-if(pb>0)sum(db)/sum(post[,b]) else NA_real_
  influence<-if(pa>0&&pb>0)post[,a]*(da-ma)/pa-post[,b]*(db-mb)/pb else rep(NA_real_,n)
  data.frame(first=names[a],second=names[b],post_difference=mean(gap),paired_se=sd(gap)/sqrt(n),cadd_difference=ma-mb,cadd_difference_se=if(all(is.finite(influence)))sd(influence)/sqrt(n) else NA_real_)
 }))
 list(events=per,paired=paired)
}
save_mc_event_audit <- function(reps,methods,k_star,horizon,tag,fields=c('stop_asym','stop_adj')) {
 directory<-file.path(jbes_root(),'results','event_audit');dir.create(directory,recursive=TRUE,showWarnings=FALSE)
 saveRDS(list(tag=tag,k_star=k_star,horizon=horizon,replications=reps),file.path(directory,paste0(tag,'.rds')))
 for(field in fields){
  st<-sapply(methods,function(method)vapply(reps,function(z)unname(z[[method]][field]),numeric(1)))
  colnames(st)<-methods;a<-mc_event_audit(st,k_star,horizon)
  write.csv(a$events,file.path(directory,paste0(tag,'_',field,'_events.csv')),row.names=FALSE)
  write.csv(a$paired,file.path(directory,paste0(tag,'_',field,'_paired.csv')),row.names=FALSE)
 }
}
