Fold.Enrich <- function(enrich_table){
  enrich_table$fold_enrichment <- enrich_table$KEGG_PATHWAY_in_list/enrich_table$expected
  
  return(enrich_table)
}