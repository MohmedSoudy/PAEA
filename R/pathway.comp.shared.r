pathway.comp.shared.pvalue <- function(enr.result , ass.result, org.data){
  
  enr.result <- enr.result %>% select(KEGG_PATHWAY_description , p_adj)
  colnames(enr.result)[2] <- c("Enrichment")
  enr.result$Enrichment <- -log10(enr.result$Enrichment)
  
  
  ass.result <- ass.result %>% select(KEGG_PATHWAY_description,  p_adj)
  colnames(ass.result)[2] <- c("Association")
  ass.result$Association <- -log10(ass.result$Association)
  
  
  
  both.result <- merge(enr.result, ass.result , by = "KEGG_PATHWAY_description")
  both.result <- reshape2::melt( data.frame(both.result) )
  
  
  
  p.value.plot <- ggplot(both.result, aes(y = value, x = KEGG_PATHWAY_description , fill = variable))+
    
    geom_bar(stat="identity", position=position_dodge(), width  = 0.5) + 
    
    theme_classic() + coord_flip() + ylab("p-adj. Value") + xlab("") +
    
    geom_hline(yintercept= -log10(0.05), linetype="dashed", color = "red", size=1 , alpha= 0.7)

  p.value.plot <- p.value.plot + scale_fill_manual(values=c('#99A3A4','#17202A')) + 
    guides(fill=guide_legend(title="")) + 
    
    theme(text = element_text(size=15, color  = "black", face = "bold"),
         axis.text = element_text(size=10, color  = "black", face = "bold"))

  return( p.value.plot )

}

pathway.comp.shared.fold <- function(enr.result , ass.result, org.data){

  enr.result <- enr.result %>% select(KEGG_PATHWAY_description , fold_enrichment)
  colnames(enr.result)[2] <- c("Enrichment")

  
  ass.result <- ass.result %>% select(KEGG_PATHWAY_description,  fold_enrichment)
  colnames(ass.result)[2] <- c("Association")

  
  
  both.result <- merge(enr.result, ass.result , by = "KEGG_PATHWAY_description")
  both.result <- reshape2::melt( data.frame(both.result) )
  
  
  
  fold.value.plot <- ggplot(both.result, aes(y = value, x = KEGG_PATHWAY_description , fill = variable))+
    
    geom_bar(stat="identity", position=position_dodge(), width  = 0.5) + 
    
    theme_classic() + coord_flip() + ylab("Fold Value") + xlab("") 
    

    fold.value.plot <- fold.value.plot + scale_fill_manual(values=c('#99A3A4','#17202A')) + 
    guides(fill=guide_legend(title="")) + 
    
    theme(text = element_text(size=15, color  = "black", face = "bold"),
          axis.text = element_text(size=10, color  = "black", face = "bold"))
  
  return( fold.value.plot )
  
}