

Run.PEA <- function(cpd.search.input, merged.df, N = 2, method = "BH"){
  
  pathway.df.merged <- merged.df
  
  pathwaykegg.to.cpd <- merged.df[,1:2]
  #pathwaykegg.to.cpd_ <- merged.df %>% dplyr:select(1:2)
  
  pathwaykegg.to.name <- merged.df[,c(1,4)]
  
  
  pathwaykegg.to.cpd.melted <- data.frame()
  

  for (i_cpd_index in 1:nrow(pathwaykegg.to.cpd)){
    
    i_cpd <- pathwaykegg.to.cpd[i_cpd_index,]
    
    cpd.split <- str_split(i_cpd$cpd, ";")[[1]]
    
    replicate.Ids <- replicate(length(cpd.split), i_cpd$IDs_map)
    
    df.temp <- data.frame(IDs_map = replicate.Ids , cpd_ids = cpd.split)
    
    pathwaykegg.to.cpd.melted <- rbind(df.temp, pathwaykegg.to.cpd.melted)
  }
  

  
  all.unique.pathway <- unique(pathwaykegg.to.cpd.melted$cpd_ids)
  all.unique.cpd <- unique(cpd.search.input)
  
  
  
  all_KEGG_cnt <- pathwaykegg.to.cpd.melted %>% dplyr::filter(.data$cpd_ids %in% 
                                                                all.unique.pathway) %>% dplyr::group_by(.data$IDs_map) %>% dplyr::summarize(KEGG_cnt = length(.data$cpd_ids))
  

  sig_KEGG_cnt <- pathwaykegg.to.cpd.melted %>% dplyr::filter(.data$cpd_ids %in% 
                                                                all.unique.cpd) %>% dplyr::group_by(.data$IDs_map) %>% dplyr::summarize(KEGG_in_list = length(.data$cpd_ids))
  
  
  enrich_table <- merge(all_KEGG_cnt, sig_KEGG_cnt, by = "IDs_map", 
                        all = TRUE)
  
  
  enrich_table$KEGG_in_list[is.na(enrich_table$KEGG_in_list)] = 0
  
  
  enrich_table <- enrich_table %>% dplyr::filter(.data$KEGG_in_list >= 
                                                   N) %>% dplyr::mutate(numTested = length(all.unique.pathway), numSig = length(all.unique.pathway[which(all.unique.pathway %in% 
                                                                                                                                                           all.unique.cpd)]), expected = (.data$numSig/.data$numTested) * 
                                                                          .data$KEGG_cnt)
  
  enrich_table <- merge(pathwaykegg.to.name, enrich_table, 
                        by = "IDs_map")
  
  
  enrich_table$enrich_p <- NA
  
  
  for (i in 1:nrow(enrich_table)) {
    a = enrich_table[i, "KEGG_in_list"]
    b = enrich_table[i, "KEGG_cnt"] - enrich_table[i, "KEGG_in_list"]
    c = enrich_table[i, "numSig"] - enrich_table[i, "KEGG_in_list"]
    d = enrich_table[i, "numTested"] - a - b -c 
    enrich_table$enrich_p[i] = fisher.test(matrix(c(a, b, c, d), nrow = 2), 
                                           alternative = "greater")$p.value
    
    
  }

  
  enrich_table <- enrich_table[order(enrich_table$enrich_p), 
                               ]
  
  #methods  c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
  #   "fdr", "none")
  
  enrich_table$fdr <- p.adjust(enrich_table$enrich_p, 
                               method = method)

  colnames(enrich_table) <- c("KEGG_PATHWAY_ID", "KEGG_PATHWAY_description", 
                              "KEGG_PATHWAY_cnt", "KEGG_PATHWAY_in_list", "KEGG_DATABASE_cnt", 
                              "KEGG_DATABASE_in_list", "expected", "enrich_p", "p_adj")
  
  species <- unlist(strsplit(as.character(enrich_table$KEGG_PATHWAY_description[1]), 
                             split = " - ", fixed = TRUE))[2]
  org.name.only <- trimws(species, which = "both")
  
  enrich_table$org_name <- org.name.only
  
  return (enrich_table)
  
}