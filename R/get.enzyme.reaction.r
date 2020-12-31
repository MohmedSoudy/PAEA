get.enzymes <- function(sig.map, org.df){
  

  df.all <- data.frame()
  org.df.pathways.names <- org.df %>% select(IDs_map, Name)
  
  
  for (i_pathway in sig.map){ 
    
    url.info <- paste0("https://www.genome.jp/dbget-bin/get_linkdb?-t+enzyme+path:",i_pathway)
    
    url.info <- fread(url.info, fill=T , quote = "" , sep=",", showProgress = T, header = F )
    
    colnames(url.info) <- "V1"
    
    
    url.info <- url.info[which(grepl("Definition",url.info$V1, fixed = T)):nrow(url.info),1]
    
    
    get.info <- url.info[grepl("ec:", url.info$V1 , fixed = T), ]
    
    
    parse.info <- data.frame(str_split_fixed(get.info$V1 , fixed(" "), n = 2))
    
    parse.info <- data.frame(str_split_fixed(parse.info$X2  , fixed(" "), n = 2))
    
    colnames(parse.info)[2] <- "Definition"
    
    parse.info$Definition <- trimws(parse.info$Definition)
    
    parse.info$enzyme.id <- str_remove_all(unlist(str_extract_all(parse.info$X1, "ec.*\">")), fixed('\">'))
    
    parse.info$X1 <- NULL
    
    parse.info$pathway <- org.df.pathways.names[org.df.pathways.names$IDs_map == i_pathway]$Name
    
    parse.info$pathway.IDs <- i_pathway
    
    parse.info <- parse.info %>% select(pathway.IDs, pathway, enzyme.id , Definition)
    
    df.all <- rbind(df.all, parse.info)
  }
   return(df.all)
}

#unlist(str_extract_all(url.info$V1, "(ec:[0-9]+)"))                      


get.reaction <- function(sig.map, org.df){
  
  
  df.all <- data.frame()
  org.df.pathways.names <- org.df %>% select(IDs_map, Name)
  
  
  for (i_pathway in sig.map){ 
    
    url.info <- paste0("https://www.genome.jp/dbget-bin/get_linkdb?-t+reaction+path:",i_pathway)
    
    url.info <- fread(url.info, fill=T , quote = "" , sep=",", showProgress = T, header = F )
    
    colnames(url.info) <- "V1"
    
    
    url.info <- url.info[which(grepl("Definition",url.info$V1, fixed = T)):nrow(url.info),1]
    
    
    get.info <- url.info[grepl("rn:", url.info$V1 , fixed = T), ]
    
    
    parse.info <- data.frame(str_split_fixed(get.info$V1 , fixed(" "), n = 2))
    
    parse.info <- data.frame(str_split_fixed(parse.info$X2  , fixed(" "), n = 2))
    
    colnames(parse.info)[2] <- "Definition"
    
    parse.info$Definition <- trimws(parse.info$Definition)
    
    
    parse.info$reaction.id <- str_remove_all(unlist(str_extract_all(parse.info$X1, "rn.*\">")), fixed('\">'))
    
    parse.info$X1 <- NULL
    
    parse.info$pathway <- org.df.pathways.names[org.df.pathways.names$IDs_map == i_pathway]$Name
    
    parse.info$pathway.IDs <- i_pathway
    
    parse.info <- parse.info %>% select(pathway.IDs, pathway, reaction.id , Definition)
    
    df.all <- rbind(df.all, parse.info)
  }
  return(df.all)
}