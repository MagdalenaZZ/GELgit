require(opencgaR)
require(getPass)
require(dplyr)
require(ggplot2)
require(lubridate)
options(stringsAsFactors=F,scipen=12)

username="mzarowiecki" ##change here

ifnull=function(x) {ifelse(length(x)!=8,NA,x)}
ifnull2=function(x) {ifelse(is.null(x),NA,x)}
ifnull3=function(x) {if(inherits(x,"try-error")) {return(NA)} else {return(x)}}
ifGEL_METRICS=function(x,y) {
  if(length(x)!=8) {
    return(NA)
  } else {
    return(x[[y]])
  }
}

filename=today()

#connect to OpenCGA
con=initOpencgaR(host = "https://opencgainternal.gel.zone/opencga/",version="v1",user=username)


studies=c(study_cs=1000000038, #ddf writes here
            study_cg=1000000034, #cancer germline 38
	    ) 

limit=NULL #NULL for all otherwise change

stats_list=vector("list",length(studies))
for(i in 1:length(studies)) {
  
  con=opencgaLogin(opencga = con,userid=username,passwd = getPass())
    study=names(studies)[i]
  stats_list[[study]]=fileClient(OpencgaR = con,action="search",params=list(study=as.character(studies[i]),format="BAM",include="creationDate,status.date,size,sampleIds,stats.SAMTOOLS_STATS_FILTERED,stats.SAMTOOLS_STATS_ALL,stats.WHOLE_GENOME_COVERAGE,stats.CDS_COVERAGE,stats.ILLUMINA_SUMMARY_REPORT,stats.ILLUMINA_SUMMARY_REPORT_CANCER,stats.GEL_METRICS,name,stats.BAM_HEADER_MACHINE,stats.AT_GC_DROP,stats.VERIFY_BAM_ID.FREEMIX,attributes.runId",status="READY",limit=limit)) # retrieve all records
  save(stats_list,file=paste0("/Users/magz/",filename,"-stats_list.RData")) #change to your path here
  gc()
  }







  stats=stats_list[["study_cs"]]
  out=lapply(1:nrow(stats),function(i) {
  x=stats[i,]
  coverageSummary=x$stats$WHOLE_GENOME_COVERAGE$coverageSummary[[1]]
  data.frame(date=parse_date_time(x$creationDate,"YmdHMS"),
             sample=ifnull2(ifnull3(try(x$stats$ILLUMINA_SUMMARY_REPORT_CANCER[[1]]$tumour_stats$SAMPLE_ID,T))),
             plate=ifnull2(ifnull3(try(sub(pattern="(LP.*)-.*","\\1",x$stats$ILLUMINA_SUMMARY_REPORT_CANCER[[1]]$tumour_stats$SAMPLE_ID),T))),
             chimeric=x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_PAIRS_ON_DIFFERENT_CHROMOSOMES/x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_RAW_TOTAL_SEQUENCES,
             contamination=x$stats$VERIFY_BAM_ID$FREEMIX,
             mapped=x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_READS_MAPPED/x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_RAW_TOTAL_SEQUENCES,
             insert=x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_INSERT_SIZE_AVERAGE,
             dups=x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_READS_DUPLICATED/x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_RAW_TOTAL_SEQUENCES,
             gbq30=x$stats$GEL_METRICS$GbQ30NoDupsNoClip,
             gt15x=x$stats$GEL_METRICS$perc_bases_ge_15x_mapQ_ge11,
             gc=x$stats$AT_GC_DROP$gc_drop,
             at=x$stats$AT_GC_DROP$at_drop,
             error=x$stats$SAMTOOLS_STATS_FILTERED$SAMTOOLS_ERROR_RATE,
             coverage=ifnull2(coverageSummary$med[coverageSummary$scope=="autosomes"]),
             rmsd=ifnull2(coverageSummary$localRMSD[coverageSummary$scope=="autosomes"]),
             snvs=ifnull2(ifnull3(try(x$stats$ILLUMINA_SUMMARY_REPORT_CANCER[[1]]$pair_stats$SOMATIC_SNVS,T))),
             indels=ifnull2(ifnull3(try(x$stats$ILLUMINA_SUMMARY_REPORT_CANCER[[1]]$pair_stats$SOMATIC_INDELS,T))),
             cnvs=ifnull2(ifnull3(try(x$stats$ILLUMINA_SUMMARY_REPORT_CANCER[[1]]$pair_stats$SOMATIC_CNVS,T))),
             svs=ifnull2(ifnull3(try(x$stats$ILLUMINA_SUMMARY_REPORT_CANCER[[1]]$pair_stats$SOMATIC_SV_BREAKENDS,T))),
             diversity=ifnull2(ifnull3(try(x$stats$ILLUMINA_SUMMARY_REPORT_CANCER[[1]]$tumour_stats$DIVERSITY,T))),
             purity=ifnull2(ifnull3(try(x$stats$ILLUMINA_SUMMARY_REPORT_CANCER[[1]]$pair_stats$ESTIMATED_PURITY,T)))
             )
  })
  out.cs.df=bind_rows(out)
  out.cs.df=mutate(out.cs.df,SNR=coverage/error)  
  
  
  stats=stats_list[["study_cg"]]
  out=lapply(1:nrow(stats),function(i) {
    x=stats[i,]
    coverageSummary=x$stats$WHOLE_GENOME_COVERAGE$coverageSummary[[1]]
    data.frame(date=parse_date_time(x$creationDate,"YmdHMS"),
               sample=x$stats$ILLUMINA_SUMMARY_REPORT$SAMPLE_ID,
               plate=sub(pattern="(LP.*)-.*","\\1",x$stats$ILLUMINA_SUMMARY_REPORT$SAMPLE_ID),
               chimeric=x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_PAIRS_ON_DIFFERENT_CHROMOSOMES/x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_RAW_TOTAL_SEQUENCES,
               contamination=x$stats$VERIFY_BAM_ID$FREEMIX,
               mapped=x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_READS_MAPPED/x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_RAW_TOTAL_SEQUENCES,
               insert=x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_INSERT_SIZE_AVERAGE,
               dups=x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_READS_DUPLICATED/x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_RAW_TOTAL_SEQUENCES,
               gbq30=ifGEL_METRICS(x$stats$GEL_METRICS[[1]],"GbQ30NoDupsNoClip"),
               gt15x=ifGEL_METRICS(x$stats$GEL_METRICS[[1]],"perc_bases_ge_15x_mapQ_ge11"),
               gc=x$stats$AT_GC_DROP$gc_drop,
               at=x$stats$AT_GC_DROP$at_drop,
               error=x$stats$SAMTOOLS_STATS_FILTERED$SAMTOOLS_ERROR_RATE,
               coverage=ifnull2(coverageSummary$med[coverageSummary$scope=="autosomes"]),
               rmsd=ifnull2(coverageSummary$localRMSD[coverageSummary$scope=="autosomes"]),
               snvs=x$stats$ILLUMINA_SUMMARY_REPORT$SNVS,
               indels=x$stats$ILLUMINA_SUMMARY_REPORT$INDELS,
               cnvs=x$stats$ILLUMINA_SUMMARY_REPORT$CNVS,
               svs=x$stats$ILLUMINA_SUMMARY_REPORT$SV_BREAKENDS,
               diversity=x$stats$ILLUMINA_SUMMARY_REPORT$DIVERSITY,
               purity=NA
    )
  })
  out.cg.df=bind_rows(out)
  out.cg.df=mutate(out.cg.df,SNR=coverage/error)
  
  
  stats=stats_list[["study_g8"]]
  out=lapply(1:nrow(stats),function(i) {
    x=stats[i,]
    coverageSummary=x$stats$WHOLE_GENOME_COVERAGE$coverageSummary[[1]]
    data.frame(date=parse_date_time(x$creationDate,"YmdHMS"),
               sample=x$stats$ILLUMINA_SUMMARY_REPORT$SAMPLE_ID,
               plate=sub(pattern="(LP.*)-.*","\\1",x$stats$ILLUMINA_SUMMARY_REPORT$SAMPLE_ID),
               chimeric=x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_PAIRS_ON_DIFFERENT_CHROMOSOMES/x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_RAW_TOTAL_SEQUENCES,
               contamination=x$stats$VERIFY_BAM_ID$FREEMIX,
               mapped=x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_READS_MAPPED/x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_RAW_TOTAL_SEQUENCES,
               insert=x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_INSERT_SIZE_AVERAGE,
               dups=x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_READS_DUPLICATED/x$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_RAW_TOTAL_SEQUENCES,
               gbq30=ifGEL_METRICS(x$stats$GEL_METRICS[[1]],"GbQ30NoDupsNoClip"),
               gt15x=ifGEL_METRICS(x$stats$GEL_METRICS[[1]],"perc_bases_ge_15x_mapQ_ge11"),
               gc=x$stats$AT_GC_DROP$gc_drop,
               at=x$stats$AT_GC_DROP$at_drop,
               error=x$stats$SAMTOOLS_STATS_FILTERED$SAMTOOLS_ERROR_RATE,
               coverage=ifnull2(coverageSummary$med[coverageSummary$scope=="autosomes"]),
               rmsd=ifnull2(coverageSummary$localRMSD[coverageSummary$scope=="autosomes"]),
               snvs=x$stats$ILLUMINA_SUMMARY_REPORT$SNVS,
               indels=x$stats$ILLUMINA_SUMMARY_REPORT$INDELS,
               cnvs=x$stats$ILLUMINA_SUMMARY_REPORT$CNVS,
               svs=x$stats$ILLUMINA_SUMMARY_REPORT$SV_BREAKENDS,
               diversity=x$stats$ILLUMINA_SUMMARY_REPORT$DIVERSITY,
               purity=NA
    )  
  })
  out.g8.df=bind_rows(out)
  out.g8.df=mutate(out.g8.df,SNR=coverage/error)

  out.df=rbind(cbind(out.cs.df,study="Cancer Somatic"),cbind(out.cg.df,study="Cancer Germline"),cbind(out.g8.df,study="RD Germline 38"))  


pdf( file=paste0("c:/Users/arendon/OneDrive/Desktop/Metrics-",filename,".pdf"),paper = "a4r")

plot.new()
text(x=.5,y=.5,paste0("Sequencing Metrics Report\n",filename))

plot.new()
text(x=.5,y=.5,"Trends")

ggplot(out.df,aes(date,gbq30,colour=study))+geom_point()+geom_smooth(level=0.99)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_hline(yintercept = c(85,212.5),colour=c(4,3))

ggplot(out.df,aes(date,gt15x,colour=study))+geom_point()+geom_smooth(level=0.99)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_hline(yintercept = 95)

ggplot(out.df,aes(date,contamination,colour=study))+geom_point(level=0.99)+geom_smooth()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_hline(yintercept = 0.03)
  
ggplot(out.df,aes(date,mapped,colour=study))+geom_point()+geom_smooth(level=0.99)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(out.df,aes(date,chimeric,colour=study))+geom_point()+geom_smooth(level=0.99)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylim(0,0.025)
  
ggplot(out.df,aes(date,dups,colour=study))+geom_point()+geom_smooth(level=0.99)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(out.df,aes(date,insert,colour=study))+geom_point()+geom_smooth(level=0.99)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(out.df,aes(date,gc,colour=study))+geom_point()+geom_smooth(level=0.99)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(out.df,aes(date,at,colour=study))+geom_point()+geom_smooth(level=0.99)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(out.df,aes(date,SNR,colour=study))+geom_point()+geom_smooth(level=0.99)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(out.df,aes(date,coverage,colour=study))+geom_point()+geom_smooth(level=0.99)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(out.df,aes(date,log10(snvs),colour=study))+geom_point()+geom_smooth(level=0.99)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(out.df,aes(date,log10(indels),colour=study))+geom_point()+geom_smooth(level=0.99)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(out.df,aes(date,log10(svs),colour=study))+geom_point()+geom_smooth(level=0.99)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(out.df,aes(date,log10(cnvs),colour=study))+geom_point()+geom_smooth(level=0.99)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(out.df,aes(date,diversity,colour=study))+geom_point()+geom_smooth(level=0.99)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

out.df.plate=out.df %>%
  group_by(plate,study) %>%
  summarise(mean(chimeric),mean(contamination),mean(gc),mean(insert),mean(snvs),mean(cnvs),mean(svs),mean(date))

out.df.plate=out.df.plate[order(out.df.plate$`mean(date)`),]
out.df.plate$plate=factor(out.df.plate$plate,ordered = T)


plot.new()
text(x=.5,y=.5,"Plate means")
ggplot(out.df.plate,aes(plate,`mean(chimeric)`))+geom_col()+facet_grid(study~.)
ggplot(out.df.plate,aes(plate,`mean(contamination)`))+geom_col()+facet_grid(study~.)
ggplot(out.df.plate,aes(plate,`mean(gc)`))+geom_col()+facet_grid(study~.)
ggplot(out.df.plate,aes(plate,`mean(insert)`))+geom_col()+facet_grid(study~.)
ggplot(out.df.plate,aes(plate,`mean(svs)`))+geom_col()+facet_grid(study~.)
ggplot(out.df.plate,aes(plate,`mean(cnvs)`))+geom_col()+facet_grid(study~.)
ggplot(out.df.plate,aes(plate,`mean(svs)`))+geom_col()+facet_grid(study~.)

dev.off()


write.csv(out.df,file=paste0("c:/Users/arendon/OneDrive/Desktop/Metrics-",filename,".csv"))  #change path here


# #read bertha
# 
# 
# con=initOpencgaR(host = "https://opencgainternal.gel.zone/opencga/",version="v1",user=username)
# con=opencgaLogin(opencga = con,userid=username,passwd = getPass())
# 
# studies=c(study_cs=1000000038, #ddf writes here
#           study_cg=1000000034, #cancer germline 38
#           study_g8=1000000032) #germline 38
# 
# limit=NULL #NULL for all otherwise change
# 
# token=con@sessionId
# 
# url_base="https://opencgainternal.gel.zone/opencga/webservices/rest/v1/files/search?name=params.json&include=creationDate,attributes.runId&status=READY&limit="
# i=0
# page=5000
# dates=vector("list",40)
# keep_going=T
# study="study_g8"
# while(keep_going){
#   refresh=content(httr::POST(paste0("https://opencgainternal.gel.zone/opencga/webservices/rest/v1/users/arendon/login?sid=",token),httr::add_headers("Content-Type"="application/json","Accept"="application/json","Authorisation"="Bearer"),body="{}"))
#   token=refresh$response[[1]]$result[[1]]$token
#   out=fromJSON(paste0(url_base,page,"&sid=",token,"&skip=",as.character(i*page),"&study=",as.character(studies[study])))$response
#   if(out$numResults==0) {
#     keep_going=F
#   } else {
#     dates[[i+1]]=out$result[[1]]
#     i=i+1
#     cat(paste0(i," "))
#   }
# }
# 
# dates_df=bind_rows(lapply(dates,function(x) {data.frame(runId=as.character(x$attributes[[1]]),creationDate=parse_date_time(x$creationDate,"YmdHMS"))}))
# 
# dates_df_cs=dates_df
# 
# 
# dates_df_dedup=dates_df %>%
#   group_by(runId) %>%
#   slice(which.min(creationDate))
# 
# load("c:/Users/arendon/OneDrive/Desktop/Bertha-20180605.RData")
# 
# bertha_sum_df=merge(bertha_sum_df,dates_df_dedup,by.x="x.workflow_run_id",by.y="runId",all.x=T,all.y=F)
# 
# 
# 
# baseurl="http://bio-prod-bertha-orch-01.gel.zone"
# next_page="/api/1/workflowrunsummary/?page=1&page_size=100"
# bertha=vector("list",1200)
# page=1
# while(!is.null(next_page)){
#   input=fromJSON(paste0(baseurl,next_page))
#   next_page=input$pagination[['next']]
#   bertha[[page]]=input$data
#   page=page+1
#   message(paste("page:",page))
# }
# 
# save(bertha,file="~/bertha-2018-06-12.RData")
# save.image(file="c:/Users/arendon/test.RData")
# 
# save(stats,file=paste0("c:/Users/Augusto Rendon/OneDrive/Desktop/Metrics-",filname,".rda"))
# 
# for(i in 660:1200) {bertha[[i]]=bertha2[[i]]}
# rm(bertha2)
# 
# 
# bertha_sum=vector("list",1200)
# 
# for(i in 1:length(bertha)){
#   x=bertha[[i]]
#   bertha_sum[[i]]=data.frame(x$workflow_run_id,x$status,x$workflow_name,x$total_walltime,x$workflow_run$bertha_version$version_name)
# }
# 
# bertha_sum_df=bind_rows(bertha_sum)
# bertha_sum_df$x.total_walltime_dur=duration(sub(pattern="(.*):(.*):(.*)",replacement = "\\1H\\2M\\3S",x=bertha_sum_df$x.total_walltime))
# 
# load("~/../Desktop/Bertha-20180605.RData")
# 
# table(bertha_sum_df$x.status)
# ggplot(subset(bertha_sum_df,x.status=="completed" & x.workflow_name %in% c("cancer_1-0-0","cancer-no-dispatch_1-0-0")),aes(x.workflow_run_id,as.numeric(x.total_walltime_dur)/3600/24,colour=x.workflow_name))+geom_point()+ylim(0,4)
# 
# ggplot(subset(bertha_sum_df, x.workflow_name %in% c("cancer_1-0-0","cancer-no-dispatch_1-0-0")),aes(x.workflow_name,as.numeric(x.total_walltime_dur)/3600/24))+geom_violin()+ylim(0,4)
#   
# ggplot(subset(bertha_sum_df,x.status%in%c("completed","failed") & x.workflow_name %in% c("cancer_1-0-0","cancer-no-dispatch_1-0-0")),aes(x.workflow_run_id,colour=x.status))+geom_ma(n=50)
# 
