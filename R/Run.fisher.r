Run.fisher <- function(list.enr1 , list.enr2 , method = "BH"){
  
  if (dim(list.enr1)[2] != 11 & dim(list.enr2)[2]){
    return(print("list.enr1 and list.enr is out of dimension"))
  }
  both.list1.list2 <- merge(list.enr1, list.enr2 ,by = c("KEGG_PATHWAY_ID", "KEGG_PATHWAY_description", 
                                                                "KEGG_PATHWAY_cnt", "KEGG_DATABASE_cnt"))
  
  both.list1.list2 <- both.list1.list2 %>% dplyr::filter(!(.data$KEGG_PATHWAY_in_list.x == 
                                                               0 & .data$KEGG_PATHWAY_in_list.y == 0))
  
  paths.removed <- dim(both.list1.list2)[1] - dim(both.list1.list2)[1]
  
  if (min(list.enr1$KEGG_PATHWAY_in_list) == 0 & min(list.enr2$KEGG_PATHWAY_in_list) == 0) {
    warning(paste0("KEGG pathways that do not contain any genes from either list provided will be removed as these cannot be tested.", 
                   paths.removed, " pathways were removed."))
  }
  
  
  colnames(both.list1.list2) <- gsub(".x", "1", colnames(both.list1.list2), 
                                    fixed = TRUE)
  colnames(both.list1.list2) <- gsub(".y", "2", colnames(both.list1.list2), 
                                    fixed = TRUE)
  
  list.12.fisher <- cbind(both.list1.list2, do.call("rbind", apply(both.list1.list2[, c("KEGG_PATHWAY_in_list2", 
                                                 "KEGG_PATHWAY_in_list1", "KEGG_DATABASE_in_list2", "KEGG_DATABASE_in_list1")], 
                                          1, function(a) {
                                            fisher.test.enr(a[1], a[2], a[3], a[4])
                                          })))
  
  colnames(list.12.fisher)[19:21] <- c("diff.odds.ratio","diff.p.value", "diff.p.adj")
  
  row.names(list.12.fisher) <- list.12.fisher$KEGG_PATHWAY_ID
  
  return(list.12.fisher)
  
}