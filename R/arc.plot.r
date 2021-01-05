arc.pathways.plot <- function(related.info, pathway.sig, top.n){
  
  related.info <- related.info[!related.info$related.pathways == "*",]
  
  pathway.sig <- pathway.sig[pathway.sig$KEGG_PATHWAY_ID %in% related.info$IDs_map,]

  if (top.n > nrow(pathway.sig)){
    
    pathway.sig.top <- pathway.sig
  }else{
    pathway.sig.top <- top_n(x = pathway.sig, n = -top.n,
                             wt =p_adj )
  }
  
  related.info <- related.info[related.info$IDs_map %in% pathway.sig.top$KEGG_PATHWAY_ID,]
  
  network.df <- NULL
  
  network.df <- data.frame()
  
  related.info <- related.info[!is.na(related.info$related.pathways),]
  
  related.info <- related.info[!related.info$disease == "",]
  
  related.info$pathway_name <- trimws(str_split_fixed(related.info$Name, fixed(" - "),2)[,1])
  
  
  for(i_pathway in related.info$pathway_name){
    
    df.temp <- related.info[related.info$pathway_name == i_pathway,]
    
    net.to <- unlist( str_split(df.temp$related.pathways, fixed(";")) )
    
    new.from <- rep(i_pathway, length(net.to))
    
    df.temp.temp <- data.frame(source = new.from, target = net.to)
    
    network.df <- rbind(network.df , df.temp.temp)
  }
  
  
  mygraph <- graph_from_data_frame(network.df)
  
  pathway.arc.plot <- ggraph(mygraph, layout="linear") + 
    geom_edge_arc(edge_colour="black", edge_alpha=0.3, edge_width=0.25) +
    geom_node_point( color="#A93226", size=7, alpha=0.5) +
    geom_node_text( aes(label=name), repel = FALSE,  size=4, color="black", nudge_y=-0.4, angle=60, hjust=1) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0.8,0), "null"),
      panel.spacing=unit(c(0,0,3.4,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-6.5, 1.2)) 
  
  return(pathway.arc.plot)
}


arc.disease.plot <- function(related.info, pathway.sig, top.n){
  
  related.info <- related.info[!related.info$disease == "*",]
  
  pathway.sig <- pathway.sig[pathway.sig$KEGG_PATHWAY_ID %in% related.info$IDs_map,]
  
  if (top.n > nrow(pathway.sig)){
    
    pathway.sig.top <- pathway.sig
  }else{
    pathway.sig.top <- top_n(x = pathway.sig, n = -top.n,
                             wt =p_adj )
  }
  
  related.info <- related.info[related.info$IDs_map %in% pathway.sig.top$KEGG_PATHWAY_ID,]
  
  network.df <- NULL
  
  network.df <- data.frame()
  
  related.info <- related.info[!is.na(related.info$disease),]
  
  related.info <- related.info[!related.info$disease == "",]
  
  related.info$pathway_name <- trimws(str_split_fixed(related.info$Name, fixed(" - "),2)[,1])
  
  
  for(i_pathway in related.info$pathway_name){
    
    df.temp <- related.info[related.info$pathway_name == i_pathway,]
    
    net.to <- unlist( str_split(df.temp$disease, fixed(";")) )
    
    new.from <- rep(i_pathway, length(net.to))
    
    df.temp.temp <- data.frame(source = new.from, target = net.to)
    
    network.df <- rbind(network.df , df.temp.temp)
  }
  
  
  mygraph <- graph_from_data_frame(network.df)
  
  disease.arc.plot <- ggraph(mygraph, layout="linear") + 
    geom_edge_arc(edge_colour="black", edge_alpha=0.3, edge_width=0.25) +
    geom_node_point( color="#A93226", size=7, alpha=0.5) +
    geom_node_text( aes(label=name), repel = FALSE, size=4, color="black", nudge_y=-0.4, angle=60, hjust=1) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0.8,0), "null"),
      panel.spacing=unit(c(0,0,3.4,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-6.5, 1.2))
  
  return(disease.arc.plot)
  
}