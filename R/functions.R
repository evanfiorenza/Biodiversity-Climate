
#Functions to accompany an analysis of shifts in biodiversity related to climate change in North America



get_data_csv<-function(file){
 #Read data into R using  remove accidental row number columns
  #Packages: vroom, tidyverse
  dat<-read.csv(file)
  return(dat)
}

group_data_csv<-function(data){
 data2<-data%>%select(-DAY,-MONTH,-BIOMASS,-resolution,-taxon,-PLOT,-SAMPLE_DESC,-X,-`...1`,-ABUNDANCE_TYPE,-BIOMASS_TYPE)%>%group_by(CLIMATE,REALM,TAXA,StudyMethod,assemblageID,STUDY_ID_REPEAT,cell,STUDY_ID,YEAR,Species)%>%summarise(ABUNDANCE=sum(ABUNDANCE,na.rm = T))%>%as.data.frame()
  
  return(data2)
}


clean_and_pivot<-function(data,names.from,values.from,values.fill=0,assemblage.ID){
  #For complex/multi data sources, use branching to create multiple targets. Function takes data in long format and returns it in wide format with the removal of columns that have zero values
  #Packages: tidyverse
  data2<-data%>%dplyr::filter(assemblageID==assemblage.ID)
 
  Spp.Mat<-pivot_wider(data2,names_from =!! sym(names.from),values_from = !! sym(values.from),values_fill =values.fill)
  Spp.Meta<-Spp.Mat[,1:9]
  Spp.Count<-Spp.Mat[,10:length(colnames(Spp.Mat))]
  Spp.Count<-Spp.Count[,which(colSums(Spp.Count)>0)]
  Spp.Mat<-cbind(Spp.Meta,Spp.Count)
  return(Spp.Mat)
}


get_decomposed_beta_pairwise<-function(spp.matrix,index.family,spp.col.start=10){
  #Take a species abundance/presence-absence matrix and get the decomposed dissimilarity matrices. Can be supplied as a list of matrices or as a single matrix. Abundance matrix can be supplied even if looking for presence-absence metrics (sorenson, jaccard)
  #packages: vegan, betapart, tidyverse
  
      xmat<-spp.matrix
      xmat<-xmat[,spp.col.start:length(xmat[1,])]
      if(index.family %in% c("sorenson",'jaccard')){
        #convert to presence-absence matrix before passing to betapart
        xmat<-decostand(xmat,method='pa')
       decomposed<- beta.pair(xmat,index.family = index.family)
      }#pa
    
      if(index.family %in% c("bray",'ruzicka')){
        decomposed<- beta.pair.abund(xmat,index.family = index.family)
        
      }#pa
      decomposed$index.family<-index.family
      decomposed$meta.dat<-as.data.frame(spp.matrix[,1:spp.col.start-1])
return(decomposed)
  
}


extract_decomposed_beta_start_to_end<-function(decomposed.beta){
 #Extract the start to end pairwise beta values.
  #packages: tidyverse
  mat.1<-as.matrix(decomposed.beta[[1]])
  mat.2<-as.matrix(decomposed.beta[[2]])
  mat.3<-as.matrix(decomposed.beta[[3]])
  start.year<-min(decomposed.beta[[5]]$YEAR)
  end.year<-max(decomposed.beta[[5]]$YEAR)
  meta<-decomposed.beta[[5]]%>%select(-YEAR)
  if(decomposed.beta[[4]] %in% c("sorenson",'jaccard')){
   extracted<-data.frame(meta[1,],start.year=start.year,end.year=end.year,index.family=decomposed.beta[[4]],turnover=mat.1[nrow(mat.1),1],nestedness=mat.2[nrow(mat.2),1],total.disimilarity=mat.3[nrow(mat.3),1]) 
  }#pa
  
  if(decomposed.beta[[4]] %in% c("bray",'ruzicka')){
    extracted<-data.frame(meta[1,],start.year=start.year,end.year=end.year,index.family=decomposed.beta[[4]],balenced=mat.1[nrow(mat.1),1],unidirectional=mat.2[nrow(mat.2),1],total.disimilarity=mat.3[nrow(mat.3),1]) 
    
  }#pa
  return(extracted)
  
}


get_alpha_div_metrics<-function(spp.matrix,spp.col.start=10){
  #calculate and return alpha diversity metrics
  #packages: tidyverse, vegan, breakaway
  xmat<-spp.matrix
  xmat<-xmat[,spp.col.start:length(xmat[1,])]
  breakaway.est<-summary(breakaway(t(xmat)))
  
  raw.richness<-t(estimateR(xmat))
  
  shannon.div<-diversity(xmat,index="shannon")
  
  alpha.table=cbind(spp.matrix[,1:spp.col.start-1],breakaway.est,raw.richness,shannon.div)
  
  return(alpha.table)
}


get_alpha_div_trends<-function(alpha.table){
  
  model.breakaway<-betta(formula = estimate ~I(YEAR-min(YEAR)),ses=error,data=alpha.table)
  
  model.raw<-glmmTMB(S.obs~I(YEAR-min(YEAR)),data=alpha.table,family = gaussian())
raw.summary<-summary(model.raw)


model.shannon<-glmmTMB(shannon.div~I(YEAR-min(YEAR)),data=alpha.table,family = gaussian())
shannon.summary<-summary(model.shannon)
  


model.chao<-glmmTMB(S.chao1~I(YEAR-min(YEAR)),data=alpha.table,family = gaussian(),weights = 1/(se.chao1+1))
chao.summary<-summary(model.chao)


model.ace<-glmmTMB(S.ACE~I(YEAR-min(YEAR)),data=alpha.table,family = gaussian(),weights = 1/(se.ACE+1))
ace.summary<-summary(model.ace)

start.year<-min(alpha.table$YEAR,na.rm=T)
end.year<-max(alpha.table$YEAR,na.rm=T)
n.samples<-length(alpha.table$YEAR)

alpha.trends=data.frame(resamp=alpha.table[1,1],assemblageID=alpha.table[1,2],STUDY_ID=alpha.table[1,3],cell=alpha.table[1,4],start.year=start.year,end.year=end.year,n.samples=n.samples,breakaway.int.est=model.breakaway$table[1,1],breakaway.int.se=model.breakaway$table[1,2],breakaway.int.p=model.breakaway$table[1,3],breakaway.year.est=model.breakaway$table[2,1],breakaway.year.se=model.breakaway$table[2,2],breakaway.year.p=model.breakaway$table[2,3],raw.int.est=raw.summary$coefficients$cond[1,1],raw.int.se=raw.summary$coefficients$cond[1,2],raw.int.p=raw.summary$coefficients$cond[1,4],raw.year.est=raw.summary$coefficients$cond[2,1],raw.year.se=raw.summary$coefficients$cond[2,2],raw.year.p=raw.summary$coefficients$cond[2,4],shannon.int.est=shannon.summary$coefficients$cond[1,1],shannon.int.se=shannon.summary$coefficients$cond[1,2],shannon.int.p=shannon.summary$coefficients$cond[1,4],shannon.year.est=shannon.summary$coefficients$cond[2,1],shannon.year.se=shannon.summary$coefficients$cond[2,2],shannon.year.p=shannon.summary$coefficients$cond[2,4],chao.int.est=chao.summary$coefficients$cond[1,1],chao.int.se=chao.summary$coefficients$cond[1,2],chao.int.p=chao.summary$coefficients$cond[1,4],chao.year.est=chao.summary$coefficients$cond[2,1],chao.year.se=chao.summary$coefficients$cond[2,2],chao.year.p=chao.summary$coefficients$cond[2,4],ace.int.est=ace.summary$coefficients$cond[1,1],ace.int.se=ace.summary$coefficients$cond[1,2],ace.int.p=ace.summary$coefficients$cond[1,4],ace.year.est=ace.summary$coefficients$cond[2,1],ace.year.se=ace.summary$coefficients$cond[2,2],ace.year.p=ace.summary$coefficients$cond[2,4])


return(alpha.trends)
}

  
  





















