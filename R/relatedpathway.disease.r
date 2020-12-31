get.info.relatedpathway.disease <- function(sig.map, org.df){
  
  
  org.df.pathways.names <- org.df %>% select(IDs_map, Name)
  
  df.all.info <- data.frame()
  
  for (i_pathway in sig.map){
    
    info.url <- paste0("https://www.genome.jp/dbget-bin/www_bget?pathway+",i_pathway)
    info.url <- fread(info.url, fill=T , quote = "" , sep=",", showProgress = T, header = F )

    colnames(info.url) <- "V1"
  ## get related pathways
  
  if (any(grepl("Related<br>pathway",info.url[,1], fixed = T))){
    
    
    related.pathway <- info.url[which(grepl("Related<br>pathway",info.url$V1, fixed = T)):nrow(info.url),1]
    
    extract.related.pathway <- unlist(str_extract_all(related.pathway$V1, "(map[0-9]+)"))
    extract.related.pathway <- unique(extract.related.pathway)
    extract.related.pathway <- extract.related.pathway[! extract.related.pathway %in% c(i_pathway)]
    
    extract.related.pathway <- paste0(extract.related.pathway, collapse = ";")
  }else{
    extract.related.pathway <- ""
  }
  
  ## get Disease  
  
  disease.get <- unlist(str_extract_all(info.url$V1, "(H[0-9]+)"))
  disease.get <- unique(disease.get)
  disease.get <- paste0(disease.get, collapse = ";")
    
  
  ## collect
  
  df.temp <- data.frame(IDs_map = i_pathway , related.pathways = extract.related.pathway, disease = disease.get)
  
  df.all.info <- rbind(df.all.info, df.temp)
  }
  
  df.all.info <- merge(org.df.pathways.names, df.all.info , by = "IDs_map")
  return(df.all.info)
}