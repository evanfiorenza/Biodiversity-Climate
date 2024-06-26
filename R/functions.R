
#Functions to accompany an analysis of shifts in biodiversity related to climate change in North America



get_data_csv<-function(file){
 #Read data into R using vroom, remove accidental row number columns
  #Packages: vroom, tidyverse
  dat<-vroom(file,show_col_types = F)
  if(colnames(dat)[1]=='...1'){
    dat<-dat%>%select(-`...1`)}
  return(dat)
}


clean_and_pivot<-function(data,names.from,values.from,values.fill=0,resamp,assemblageID){
  #For complex/multi data sources, use branching to create multiple targets. Function takes data in long format and returns it in wide format with the removeal of columns that have zero values
  #Packages: tidyverse
  data2<-data%>%filter(resamp==resamp & assemblageID==assemblageID)
  start.col=length(colnames(data))-1
  Spp.Mat<-pivot_wider(data2,names_from =!! sym(names.from),values_from = !! sym(values.from),values_fill =values.fill)
  Spp.Meta<-Spp.Mat[,1:start.col]
  Spp.Count<-Spp.Mat[,start.col:length(colnames(Spp.Mat))]
  Spp.Count<-Spp.Count[,which(colSums(Spp.Count)>0)]
  Spp.Mat<-cbind(Spp.Meta,Spp.Count)
  return(Spp.Mat)
}


get_decomposed_beta_pairwise<-function(spp.matrix,index.family,spp.col.start=6){
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



























