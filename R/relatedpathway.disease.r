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




## plot
plot.pathway.related.bar <- function(related.pathway.info){
  
  related.pathway.info <- related.pathway.info %>% select(Name, related.pathways)
  related.pathway.info$related.count <- str_count(related.pathway.info$related.pathways , fixed(";")) +1
  
  pathway.related.plot <- ggplot(related.pathway.info, aes(y = as.factor(related.count), 
                                                           x = reorder(Name,related.count) ))+
    
    geom_bar(stat="identity", position=position_dodge(), width  = 0.5) + 
    
    theme_classic() + coord_flip() + ylab("# Of Related Pathways") + xlab("") +
    
    theme(text = element_text(size=15, color  = "black", face = "bold"),
          axis.text = element_text(size=10, color  = "black", face = "bold"))
  
  return(pathway.related.plot)

} 


plot.disease.bar <- function(related.disease.info){
  
  related.disease.info <- related.disease.info %>% select(Name, disease )
  related.disease.info$related.count <- str_count(related.disease.info$disease , fixed(";")) +1
  
  disease.related.plot <- ggplot(related.disease.info, aes(y = as.factor(related.count),
                                                           x = reorder(Name,related.count) ))+

    geom_bar(stat="identity", position=position_dodge(), width  = 0.5) + 
  
    theme_classic() + coord_flip() + ylab("# Of Related Diseases") + xlab("") +
    
    theme(text = element_text(size=15, color  = "black", face = "bold"),
          axis.text = element_text(size=10, color  = "black", face = "bold")) 
  
  return(disease.related.plot)
  
} 
