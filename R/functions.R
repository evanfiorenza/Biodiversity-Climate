
#Functions to accompany an analysis of shifts in biodiversity related to climate change in North America



get_data_csv<-function(file){
 #Read data into R using vroom, remove accidental row number columns
  #Packages: vroom, tidyverse
  dat<-vroom(file,show_col_types = F)
  if(colnames(dat)[1]=='...1'){
    dat<-dat%>%select(-`...1`)}
  return(dat)
}











