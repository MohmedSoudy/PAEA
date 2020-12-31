Do.associate <- function(pathway.df.merged, input.ids , alpha = 0.01, uni= F){
  
  split.cpd <- function(cpds){
    cpds.count <- str_count(cpds, ";") + 1
    cpds.split <- as.character(str_split_fixed(cpds, ";", cpds.count))
    return(data.frame(cpds.split))
  }
  
  if (uni == F){
  unique.ids.database <- lapply(pathway.df.merged$cpd,
                                function(x) split.cpd(x))
  }else if(uni == T){
    unique.ids.database <- lapply(pathway.df.merged$uniprot_id,
                                  function(x) split.cpd(x))
  }
  unique.ids.database <- do.call(rbind.data.frame, unique.ids.database)
  
  unique.ids.database <- unique(unique.ids.database$cpds.split)
  
  unique.ids.database <- unique.ids.database[!unique.ids.database %in% input.ids]
  
  new.input.ids <- c()
  
  for (x.ids in input.ids){
    unique.ids.database <- unique.ids.database[!unique.ids.database %in% c(x.ids)]
    
    for(z.ids in unique.ids.database){
      
    
    if (uni == F){
    p_value <- Get.associate1.1(pathway.df.merged$cpd, x.input = x.ids ,
                             z.database =  z.ids)
    }else if(uni == T){
      p_value <- Get.associate1.1(pathway.df.merged$uniprot_id, x.input = x.ids ,
                               z.database =  z.ids)
    }
    if (p_value < alpha){
      new.input.ids <- c(new.input.ids , c(z.ids))
      unique.ids.database <- unique.ids.database[!unique.ids.database %in% c(z.ids)]
    }
    
    }
    
  }
  new.input.ids <- unique(new.input.ids)
  new.input.ids <- c(new.input.ids , input.ids)
  return(new.input.ids)
}