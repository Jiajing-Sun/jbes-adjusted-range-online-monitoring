# Black-and-white figures of the complete B grid at the existing m=500, gamma=0 setting.
library(ggplot2)
src <- 'results/persistent_training'
d <- read.csv(file.path(src,'persistent_summary.csv'))
stopifnot(nrow(d)==1260,all(d$N==1000))
selected <- subset(d,m==500 & gamma==0 & calibration=='oracle' & scenario!='none')
selected$method <- factor(selected$method,levels=c('RSMS','SSMS','CSMS'))
selected$direction <- factor(sub('_(early|late)$','',selected$scenario),
 levels=c('parallel','reverse','transverse'),labels=c('Same direction','Opposite direction','Transverse direction'))
selected$timing <- ifelse(grepl('early',selected$scenario),'Early change','Late change')
selected$panel <- factor(paste(selected$timing,selected$direction,sep=': '),
 levels=c(outer(c('Same direction','Opposite direction','Transverse direction'),
 c('Early change','Late change'),function(a,b)paste(b,a,sep=': '))))
makeplot <- function(y,label,name) {
 p <- ggplot(selected,aes(x=B,y=.data[[y]],linetype=method,shape=method,group=method))+
  geom_line(colour='black',linewidth=.45,na.rm=TRUE)+geom_point(colour='black',size=1.7,na.rm=TRUE)+
  facet_wrap(~panel,ncol=3,scales='free_y')+
  scale_linetype_manual(values=c('solid','longdash','dotdash'))+
  scale_shape_manual(values=c(16,2,0))+
  scale_x_continuous(breaks=c(0,.5,1))+
  scale_y_continuous(breaks=scales::breaks_pretty(n=3),expand=expansion(mult=c(.08,.10)))+
  labs(x='Persistent training-shift magnitude B',y=label,linetype=NULL,shape=NULL)+
  theme_bw(base_size=10)+theme(legend.position='bottom',legend.key.width=grid::unit(1.4,'cm'),
   strip.text=element_text(size=8),panel.grid.minor=element_blank())
 ggsave(file.path('figures',name),p,width=8,height=5.4,device='pdf',useDingbats=FALSE)
}
makeplot('post','Post-break detection probability','persistent_power.pdf')
makeplot('CADD','CADD (monitoring observations)','persistent_cadd.pdf')
for(pair in list(c('persistent_summary.csv','persistent_training_summary.csv'),
 c('persistent_critical_values.csv','persistent_training_critical_values.csv'),
 c('persistent_geometry.csv','persistent_training_geometry.csv'),
 c('persistent_paired_comparison.csv','persistent_training_paired_comparison.csv')))
 stopifnot(file.copy(file.path(src,pair[1]),file.path('results',pair[2]),overwrite=TRUE))
