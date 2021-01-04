sanky.cpd.plot <- function(pathway.sig, org.df, top.n , cpd.atleast ){
  
  if (top.n > nrow(pathway.sig)){
    
    pathway.sig.top <- pathway.sig
  }else{
    pathway.sig.top <- top_n(x = pathway.sig, n = -top.n,
                                wt =p_adj )
  }
  
  org.df.top <- org.df[org.df$IDs_map %in% pathway.sig.top$KEGG_PATHWAY_ID,]
  org.df.top$pathway_name <- trimws(str_split_fixed(org.df.top$Name, fixed(" - "),2)[,1])
  
  network.df <- data.frame()
  for(i_pathway in org.df.top$pathway_name){
    
    df.temp <- org.df.top[org.df.top$pathway_name == i_pathway,]
    
    net.to <- unlist( str_split(df.temp$cpd, fixed(";")) )
    
    new.from <- rep(i_pathway, length(net.to))
    
    df.temp.temp <- data.frame(source = new.from, target = net.to, value = log10(length(new.from)+1))
    
    network.df <- rbind(network.df , df.temp.temp)
  }
  
  atlest.df <- network.df %>% group_by(target) %>% summarise(n = n())
  atlest.df <- atlest.df[atlest.df$n >= cpd.atleast,]
  network.df <- network.df[network.df$target %in% atlest.df$target,]
  
  
  temp <- network.df[!duplicated(network.df$source),]
  temp <- temp %>% select(source, value)
  temp$source <- seq(0,nrow(temp)-1)
  # make a nodes data frame out of all unique nodes in networkData
  nodes <- data.frame(name = unique(c(network.df$source, network.df$target)))
  
  # make a group variable where nodes in networkData$src are identified
  nodes$group <- nodes$name %in% network.df$source
  
  # make a links data frame using the indexes (0-based) of nodes in 'nodes'
  links <- data.frame(source = match(network.df$source, nodes$name) - 1,
                      target = match(network.df$target, nodes$name) - 1)
  links <- merge(links, temp, by= "source")
  
  links <- links %>% filter(!is.na(source) & !is.na(target))
  
  snakey.plot <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
                Target = "target", Value = "value", NodeID = "name",
                units = "TWh", fontSize = 12, nodeWidth = 20)
  
  return(snakey.plot)
}

