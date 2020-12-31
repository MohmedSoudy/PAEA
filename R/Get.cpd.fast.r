Get.cpd.fast <- function(pathway.df ,  save_db = ""){
  
  if ( class(pathway.df)[1] != "character" ){
    
    
  len.df <- nrow(pathway.df)
  
  if (len.df %% 2 == 1){
    len.df <- len.df / 2
    len.df1 <- round(len.df,0)
    len.df2 <- round(len.df,0) + 1 
  }else{
    len.df1 <- len.df / 2
    len.df2 <- len.df1 + 1
  }
  
  
  pathway.df.1 <- pathway.df[1:len.df1,]
  pathway.df.2 <- pathway.df[len.df2:nrow(pathway.df),]
  
  pathway.df.list <- list(pathway.df.1, pathway.df.2)
  temp.list.final <- mclapply(pathway.df.list , Get.cpd) 

  
  temp.list.final <- rbind(data.frame(temp.list.final[[1]]) 
                           ,data.frame(temp.list.final[[2]]))
  
  if (save_db != "" & class(temp.list.final) != "character") {
    write.csv(temp.list.final , paste0("get.cpd_" , save_db , ".csv"), 
              row.names=FALSE)
  }
  
  }else{
    temp.list.final <- read.csv(pathway.df , header = T)
  }
  
  return(temp.list.final)
  
  
  
}