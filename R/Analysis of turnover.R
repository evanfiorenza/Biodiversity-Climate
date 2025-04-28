library(tidyverse)
library(glmmTMB)
library(targets)


#Load Pieces of data from intermediate step
tar_load(combined.beta)
orig=read.csv("Data/BioTIME_Renamed_Gridded_12_Filtered.csv")
nond.distance1000<-read.csv("Data/BioTime_Renamed_Gridded_12_Non_Directional_Distance_All_1000.csv")
d.distance10000<-read.csv("Data/BioTime_Renamed_Gridded_12_Directional_Distance_All_10000.csv")
elev.bath<-read.csv("Data/BioTime_Renamed_Gridded_12_Elevation_Bathymetry_All_1855.csv")
d.distance5000<-read.csv("Data/BioTime_Renamed_Gridded_12_Directional_Distance_All_5000.csv")
mar.sst.had<-read.csv("Data/BioTime_Renamed_Gridded_12_Filtered_Marine_Temperature_Trends.csv")
mar.sst.er<-read.csv("Data/BioTime_Renamed_Gridded_12_Filtered_Marine_Temperature_ERSST_Trends.csv")
land.water<-read.csv("Data/BioTime_Renamed_Gridded_12_Land_Water_All_500.csv")
land.cover<-read.csv("Data/BioTime_Renamed_Gridded_12_Land_Cover_All_500.csv")
terr.temp<-read.csv("Data/BioTime_Renamed_Gridded_12_Filtered_Terrestrial_Temperature_Trends.csv")
precip<-read.csv("Data/BioTime_Renamed_Gridded_12_Filtered_Precipitation_Trends.csv")
meta<-read.csv("Data/BioTIMEMetadata_24_06_2021.csv")


#Clean and rename data
mar.sst.had.clean<-mar.sst.had%>%select(assemblageID,modelType,aicModel,modelInterceptEstimate,modelYearEstimate,modelYearSE,modelYearZ,ModelYearP)%>%rename(Type.had=modelType,aic.had=aicModel,intEst.had=modelInterceptEstimate,yearEst.had=modelYearEstimate,yearSE.had=modelYearSE,yearZ.had=modelYearZ,yearP.had=ModelYearP)

mar.sst.er.clean<-mar.sst.er%>%select(assemblageID,modelType,aicModel,modelInterceptEstimate,modelYearEstimate,modelYearSE,modelYearZ,ModelYearP)%>%rename(Type.er=modelType,aic.er=aicModel,intEst.er=modelInterceptEstimate,yearEst.er=modelYearEstimate,yearSE.er=modelYearSE,yearZ.er=modelYearZ,yearP.er=ModelYearP)

terr.temp.clean<-terr.temp%>%select(assemblageID,modelType,aicModel,modelInterceptEstimate,modelYearEstimate,modelYearSE,modelYearZ,ModelYearP)%>%rename(Type.terr=modelType,aic.terr=aicModel,intEst.terr=modelInterceptEstimate,yearEst.terr=modelYearEstimate,yearSE.terr=modelYearSE,yearZ.terr=modelYearZ,yearP.terr=ModelYearP)

precip.clean<-precip%>%select(assemblageID,modelType,aicModel,modelInterceptEstimate,modelYearEstimate,modelYearSE,modelYearZ,ModelYearP)%>%rename(Type.precip=modelType,aic.precip=aicModel,intEst.precip=modelInterceptEstimate,yearEst.precip=modelYearEstimate,yearSE.precip=modelYearSE,yearZ.precip=modelYearZ,yearP.precip=ModelYearP)

meta.clean<-meta%>%select(STUDY_ID,HABITAT,PROTECTED_AREA,BIOME_MAP,ORGANISMS,START_YEAR,END_YEAR,ABUNDANCE_TYPE)

orig.clean<-orig%>%select(assemblageID,LATITUDE,LONGITUDE)%>%group_by(assemblageID)%>%summarise(lat=mean(LATITUDE,na.rm=T),lon=mean(LONGITUDE,na.rm=T))

#Combine the various datasets
comb.dd<-left_join(combined.beta,d.distance10000)
comb.dd.nd<-left_join(comb.dd,nond.distance1000)
comb.dd.nd.eb<-left_join(comb.dd.nd,elev.bath)
comb.dd.nd.eb.had<-left_join(comb.dd.nd.eb,mar.sst.had.clean)
comb.dd.nd.eb.had.er<-left_join(comb.dd.nd.eb.had,mar.sst.er.clean)
comb.dd.nd.eb.had.er.terr<-left_join(comb.dd.nd.eb.had.er,terr.temp.clean)
comb.dd.nd.eb.had.er.terr.lc<-left_join(comb.dd.nd.eb.had.er.terr,land.cover)
comb.dd.nd.eb.had.er.terr.lc.lw<-left_join(comb.dd.nd.eb.had.er.terr.lc,land.water)
comb.dd.nd.eb.had.er.terr.lc.lw.precip<-left_join(comb.dd.nd.eb.had.er.terr.lc.lw,precip.clean)
comb.dd.nd.eb.had.er.terr.lc.lw.precip.orig<-left_join(comb.dd.nd.eb.had.er.terr.lc.lw.precip,orig.clean)
Master.Data<-left_join(comb.dd.nd.eb.had.er.terr.lc.lw.precip.orig,meta.clean)


#write.csv(Master.Data,"Outputs/Turnover_Master_23082024.csv")
Master.Data<-read.csv("Outputs/Turnover_Master_23082024.csv")
Master.Data$time<-Master.Data$end.year-Master.Data$start.year

Terrestrial.Data<-Master.Data%>%filter(REALM=="Terrestrial" | REALM=="Freshwater")
Marine.Data<-Master.Data%>%filter(REALM=="Marine")

check<-Marine.Data[which(Marine.Data$Elevation_Bathymetry>0),]

ggplot()+geom_histogram(aes(x=end.year-start.year,fill=REALM),data=Master.Data)+theme_minimal()

min(Master.Data$end.year-Master.Data$start.year)
max(Master.Data$end.year-Master.Data$start.year)



#Is turnover related to nearest ecoregion in the terrestrial realm
Terrestrial.Mod.Nearest<-glmmTMB(turnover~I(Terrestrial_Ecoregions_Distance/1000)*time+lat+end.year+(1|STUDY_ID)+(1|cell),data=Terrestrial.Data,family=ordbeta())

Terrestrial.Mod.Nearest.add<-glmmTMB(turnover~I(Terrestrial_Ecoregions_Distance/1000)+time+lat+end.year+(1|STUDY_ID)+(1|cell),data=Terrestrial.Data,family=ordbeta())

Terrestrial.Mod.Nearest.time<-glmmTMB(turnover~time+lat+end.year+(1|STUDY_ID)+(1|cell),data=Terrestrial.Data,family=ordbeta())



AIC(Terrestrial.Mod.Nearest.null,Terrestrial.Mod.Nearest.time,Terrestrial.Mod.Nearest.dist,Terrestrial.Mod.Nearest.add,Terrestrial.Mod.Nearest)


library(MuMIn)

Terr.Nearest.Avg=model.avg(Terrestrial.Mod.Nearest.time,Terrestrial.Mod.Nearest.add,Terrestrial.Mod.Nearest)

summary(Terr.Nearest.Avg)

summary(Terrestrial.Mod.Nearest)

#Is turnover related to elevation in the terrestrial realm

Terrestrial.Mod.Elevation<-glmmTMB(turnover~I(Elevation_Bathymetry/1000)*time+lat+end.year+(1|STUDY_ID)+(1|cell),data=Terrestrial.Data,family=ordbeta())


Terrestrial.Mod.Elevation.add<-glmmTMB(turnover~I(Elevation_Bathymetry/1000)+time+lat+end.year+(1|STUDY_ID)+(1|cell),data=Terrestrial.Data,family=ordbeta())



Terrestrial.Mod.Elevation.null<-glmmTMB(turnover~time+lat+end.year+(1|STUDY_ID)+(1|cell),data=Terrestrial.Data,family=ordbeta())
summary(Terrestrial.Mod.Elevation)

AIC(Terrestrial.Mod.Elevation,Terrestrial.Mod.Elevation.add,Terrestrial.Mod.Elevation.time,Terrestrial.Mod.Elevation.elev,Terrestrial.Mod.Elevation.null)


Terr.Elev.Avg=model.avg(Terrestrial.Mod.Elevation,Terrestrial.Mod.Elevation.add,Terrestrial.Mod.Elevation.null)

summary(Terr.Elev.Avg)



newdat.elev<-rbind(data.frame(Elevation_Bathymetry=seq(0,2750,10),time=10,end.year=2020,STUDY_ID=NA,cell=NA,lat=40),data.frame(Elevation_Bathymetry=seq(0,2750,10),time=20,end.year=2020,STUDY_ID=NA,cell=NA,lat=40),data.frame(Elevation_Bathymetry=seq(0,2750,10),time=30,end.year=2020,STUDY_ID=NA,cell=NA,lat=40),data.frame(Elevation_Bathymetry=seq(0,2750,10),time=40,end.year=2020,STUDY_ID=NA,cell=NA,lat=40))


pred.elev<-cbind(newdat.elev,predict(Terrestrial.Mod.Elevation,newdata = newdat.elev,se.fit = T))
pred.elev$latit=as.factor(pred.elev$lat)
ggplot(aes(x=Elevation_Bathymetry,y=plogis(fit),color=time),data=pred.elev)+geom_point(data=Terrestrial.Data,aes(y=turnover))+geom_ribbon(aes(ymin=plogis(fit-se.fit),ymax=plogis(fit+se.fit),fill=time,group=as.factor(time)),alpha=.25)+geom_line(aes(group=as.factor(time)))+scale_color_gradientn(colors = c('#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58'))+scale_fill_gradientn(colors = c('#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58'))+theme_classic()+scale_y_continuous(expand = c(0.01,0))+scale_x_continuous(expand=c(0,0))+xlab('Elevation (m)')+ylab('Species Turnover (Jaccard)')



#Is turnover related to diseatnce to nearest boundary or coastline in the marine realm
Marine.Mod.Coast.Nearest<-glmmTMB(turnover~I(Coastline_Distance/1000)*time+lat+end.year+(1|STUDY_ID)+(1|cell),data=Marine.Data[which(Marine.Data$Elevation_Bathymetry<=0),],family=ordbeta())

Marine.Mod.Coast.Nearest.add<-glmmTMB(turnover~I(Coastline_Distance/1000)+time+lat+end.year+(1|STUDY_ID)+(1|cell),data=Marine.Data[which(Marine.Data$Elevation_Bathymetry<=0),],family=ordbeta())

Marine.Mod.Coast.Nearest.null<-glmmTMB(turnover~time+lat+end.year+(1|STUDY_ID)+(1|cell),data=Marine.Data[which(Marine.Data$Elevation_Bathymetry<=0),],family=ordbeta())

Marine.Mod.Eco.Nearest<-glmmTMB(turnover~I(Marine_Ecoregions_Distance/1000)*time+lat+end.year+(1|STUDY_ID)+(1|cell),data=Marine.Data[which(Marine.Data$Elevation_Bathymetry<=0),],family=ordbeta())

Marine.Mod.Eco.Nearest.add<-glmmTMB(turnover~I(Marine_Ecoregions_Distance/1000)+time+lat+end.year+(1|STUDY_ID)+(1|cell),data=Marine.Data[which(Marine.Data$Elevation_Bathymetry<=0),],family=ordbeta())

Marine.Mod.Eco.Nearest.null<-glmmTMB(turnover~time+lat+end.year+(1|STUDY_ID)+(1|cell),data=Marine.Data[which(Marine.Data$Elevation_Bathymetry<=0),],family=ordbeta())

Marine.Nearest.Avg=model.avg(Marine.Mod.Eco.Nearest.null,Marine.Mod.Eco.Nearest.add,Marine.Mod.Eco.Nearest,Marine.Mod.Coast.Nearest.add,Marine.Mod.Coast.Nearest)

summary(Marine.Nearest.Avg)


newdat.coast.n<-rbind(data.frame(Coastline_Distance=seq(0,1300,10),time=10,end.year=2020,STUDY_ID=NA,cell=NA,lat=40),data.frame(Coastline_Distance=seq(0,1300,10),time=20,end.year=2020,STUDY_ID=NA,cell=NA,lat=40),data.frame(Coastline_Distance=seq(0,1300,10),time=30,end.year=2020,STUDY_ID=NA,cell=NA,lat=40),data.frame(Coastline_Distance=seq(0,1300,10),time=40,end.year=2020,STUDY_ID=NA,cell=NA,lat=40))

pred.coast.n<-cbind(newdat.coast.n,predict(Marine.Mod.Coast.Nearest,newdata = newdat.coast.n,se.fit = T))

ggplot(aes(x=Coastline_Distance,y=plogis(fit),color=time),data=pred.coast.n)+geom_point(data=Marine.Data,aes(y=turnover),alpha=.25)+geom_ribbon(aes(ymin=plogis(fit-se.fit),ymax=plogis(fit+se.fit),fill=time,group=as.factor(time)),alpha=.25)+geom_line(aes(group=as.factor(time)))+scale_color_gradientn(colors = c('#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58'))+scale_fill_gradientn(colors = c('#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58'))+theme_classic()+scale_y_continuous(expand = c(0.01,0))+scale_x_continuous(expand=c(0,0))+xlab("Distance to Coastline (km)")+ylab("Species Turnover (Jaccard)")

library(emmeans)








