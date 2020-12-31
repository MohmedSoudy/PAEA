Get.cpd <- function(pathway.df){
  

    pathway.cpd.df <- data.frame()
    

    for (i_pathway in pathway.df$IDs_map){ 

      cpd.url <- paste0("https://www.genome.jp/dbget-bin/get_linkdb?-t+compound+path:",i_pathway)
      
      cpd.read <- fread(cpd.url, fill=T , quote = "" , sep=",", showProgress = T, header = F )
      
      colnames(cpd.read) <- "V1"
      
      #cpd.read <- read.csv(cpd.url)
      
      #get cPD for each pathway
      
      cpd.get <- str_extract_all(cpd.read$V1, "(cpd:C[0-9]+)")
      
      cpd.join <- paste0(cpd.get[[1]], collapse = ";")
      cpd.join <- str_remove_all(cpd.join , fixed("NA"))
      cpd.join <- str_remove_all(cpd.join , fixed("NA;"))
      
      cpd.join <- str_remove_all(cpd.join, "cpd:")
      
      
      
      
      ### merge 
      df.temp <- data.frame(IDs_map = i_pathway, cpd = cpd.join)
      pathway.cpd.df <- rbind(df.temp, pathway.cpd.df)
    }
    
    merged.df <- merge(pathway.cpd.df , pathway.df , by = "IDs_map")
    
    id.org <- max(merged.df$org_acc_ids)
    
    
    merged.df$Name <- ifelse(merged.df$org_acc_ids != id.org, paste0(merged.df$Name, merged.df$org_acc_ids) , merged.df$Name )
    
    merged.df$org_acc_ids <- id.org
    
    
    
  return(merged.df)
  
}