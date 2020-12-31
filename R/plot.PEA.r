Plot.PEA <- function(enrich_table, type = "cpd"){
  
  
  
  string_remove <- paste0(" - ",unique(enrich_table$org_name))
  enrich_table$new_pathway_name <- str_remove(enrich_table$KEGG_PATHWAY_description ,
                                              string_remove)
  
  if (type == "cpd"){
    leg.label <- "Metabolites Count"
  }else if(type == "uni"){
    leg.label <- "Proteins Count"
  }
  plot<- ggplot(data=enrich_table, aes(y=reorder(new_pathway_name,fold_enrichment), x= fold_enrichment,
                                       size = factor(KEGG_PATHWAY_in_list), 
                                       color = p_adj)) +
    geom_point() +
    theme_bw() + xlab("Fold Enrichment") + ylab("KEGG Pathway") +
    theme(legend.position="right", text = element_text(face="bold"),
          axis.text = element_text(color = "black", face = "bold")) + 
    labs(size = leg.label, color = "p_adj.")+
    scale_color_gradient(low="darkblue", high="red")
  
  return(plot)
}