# Same curves and 121-point display mesh as script 43. The 2016 display
# window is extended through August to include all illustrated alarm dates.
source('00_utils.R');source('50_empirical_utils.R');jbes_set_root(getwd())
raw<-load_gbpusd_data(Sys.getenv('JBES_GBPUSD_PATH', unset='data/GBPUSD.csv'));manifest<-list()
alarms<-read.csv('results/empirical_monitoring_results.csv',stringsAsFactors=FALSE)
for(year in c(2016,2020)) {
 event<-as.Date(if(year==2016)'2016-06-23' else '2020-03-17')
 fd<-build_daily_fd(prepare_gbpusd_year(raw,year),nbasis_cap=21L,min_points=10L)
 bundle<-evaluate_fd_matrix(fd,n_grid=121L)
 j<-match(event,fd$dates)
 first<-fd$dates[max(1L,j-60L)]
 last<-if(year==2016)as.Date('2016-08-31') else fd$dates[min(length(fd$dates),j+20L)]
 keep<-which(fd$dates>=first & fd$dates<=last)
 horizon<-if(year==2016)5 else 1
 selected<-alarms[alarms$year==year & alarms$m==50 & alarms$T==horizon & alarms$gamma==0,]
 stopifnot(nrow(selected)==3L,all(selected$reject),all(as.Date(selected$first_detection_date)%in%fd$dates[keep]))
 out<-data.frame(Date=as.character(fd$dates[keep]),bundle$values[keep,])
 write.csv(out,sprintf('results/fx_%d_multimedia_mesh.csv',year),row.names=FALSE)
 manifest[[as.character(year)]]<-list(year=year,event=as.character(event),first=as.character(min(fd$dates[keep])),last=as.character(max(fd$dates[keep])),dates=length(keep),grid=121,nbasis=fd$nbasis,unit='100/Close; GBP per 100 USD',intraday_coordinate='Original recorded order mapped to [0,1]',m=50,T=horizon,gamma=0,alpha=0.05,training_first=as.character(fd$dates[1]),training_last=as.character(fd$dates[50]),monitoring_first=as.character(fd$dates[51]),monitoring_last=as.character(fd$dates[50+50*horizon]),alarms=selected[,c('method','first_detection_date')])
}
jsonlite::write_json(manifest,'results/fx_multimedia_mesh_manifest.json',pretty=TRUE,auto_unbox=TRUE)
print(manifest)
