plot.network.pathway <- function(related.info, charge, linkDistance){
  
  network.df <- NULL
  
  network.df <- data.frame()
  
  related.info <- related.info[!is.na(related.info$disease),]
  
  related.info <- related.info[!related.info$disease == "",]
  
  related.info$pathway_name <- trimws(str_split_fixed(related.info$Name, fixed(" - "),2)[,1])
    
    
  for(i_pathway in related.info$pathway_name){
    
    df.temp <- related.info[related.info$pathway_name == i_pathway,]
    
    net.to <- unlist( str_split(df.temp$related.pathways, fixed(";")) )
    
    new.from <- rep(i_pathway, length(net.to))
    
    df.temp.temp <- data.frame(src = new.from, target = net.to)
    
    network.df <- rbind(network.df , df.temp.temp)
  }
  
  
  
  ColourScale <- 'd3.scaleOrdinal()
            .domain(["Source", "Related"])
           .range(["#C0392B", "#707B7C"]);'
  
  # make a nodes data frame out of all unique nodes in networkData
  nodes <- data.frame(name = unique(c(network.df$src, network.df$target)))
  
  # make a group variable where nodes in networkData$src are identified
  nodes$group <- nodes$name %in% network.df$src
  
  # make a links data frame using the indexes (0-based) of nodes in 'nodes'
  links <- data.frame(source = match(network.df$src, nodes$name) - 1,
                      target = match(network.df$target, nodes$name) - 1)
  
  
  links <- links %>% filter(!is.na(source) & !is.na(target))
  
  network.pathway.plot <- forceNetwork(Links = links, Nodes = nodes, Source = "source",
               Target = "target", NodeID ="name", Group = "group",
               opacity = 1, opacityNoHover = 1, zoom = T ,
               charge = charge,
               linkDistance = linkDistance,
               linkWidth = 2,
               fontSize = 14,
               fontFamily = "impact",  
               linkColour = "#5D6D7E",
               colourScale = JS(ColourScale),
               legend = F,
               arrows = F)
  
  #saveWidget(network.pathway.plot, file= "y.html")
 
  return(network.pathway.plot) 
}


plot.network.disease <- function(related.info, charge, linkDistance){
  
  network.df <- NULL
  
  network.df <- data.frame()
  
  related.info <- related.info[!is.na(related.info$disease),]
  related.info <- related.info[!related.info$disease == "",]
  
  
  
  related.info$pathway_name <- trimws(str_split_fixed(related.info$Name, fixed(" - "),2)[,1])
  
  
  for(i_pathway in related.info$pathway_name){
    
    df.temp <- related.info[related.info$pathway_name == i_pathway,]
    net.to <- unlist( str_split(df.temp$disease, fixed(";")) )
    
    new.from <- rep(i_pathway, length(net.to))
    
    df.temp.temp <- data.frame(src = new.from, target = net.to)
    
    network.df <- rbind(network.df , df.temp.temp)
  }
  
  
  
  ColourScale <- 'd3.scaleOrdinal()
            .domain(["Source", "Related"])
           .range(["#C0392B", "#707B7C"]);'
  
  # make a nodes data frame out of all unique nodes in networkData
  nodes <- data.frame(name = unique(c(network.df$src, network.df$target)))
  
  # make a group variable where nodes in networkData$src are identified
  nodes$group <- nodes$name %in% network.df$src
  
  # make a links data frame using the indexes (0-based) of nodes in 'nodes'
  links <- data.frame(source = match(network.df$src, nodes$name) - 1,
                      target = match(network.df$target, nodes$name) - 1)
  
  links <- links %>% filter(!is.na(source) & !is.na(target))
  
  
  network.disease.plot <- forceNetwork(Links = links, Nodes = nodes, Source = "source",
                                       Target = "target", NodeID ="name", Group = "group",
                                       opacity = 1, opacityNoHover = 1, zoom = T ,
                                       charge = charge,
                                       linkDistance = linkDistance,
                                       linkWidth = 2,
                                       fontSize = 14,
                                       fontFamily = "impact",  
                                       linkColour = "#5D6D7E",
                                       colourScale = JS(ColourScale),
                                       legend = F,
                                       arrows = F)
  
  
  #saveWidget(network.disease.plot, file= "y.html")
  
  return(network.disease.plot) 
}

