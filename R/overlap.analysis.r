overlap.analysis.cpd <- function(enr.data, ass.data, org.data){
  
  org.data <- org.data[org.data$cpd != "",]
  #1st overlap analysis
  org.enr <- org.data[org.data$IDs_map %in% enr.data$KEGG_PATHWAY_ID,]
  org.enr.path <- org.enr$IDs_map
  org.ass <- org.data[org.data$IDs_map %in% ass.data$KEGG_PATHWAY_ID,]
  org.ass.path <- org.ass$IDs_map
  
  
  
  cpd.enr <- unlist( str_split(org.enr$cpd , ";") )
  cpd.ass <- unlist(str_split(org.ass$cpd , ";"))
  
  # venn dig. 
  
  venn.plot.cpd <- draw.pairwise.venn(area1= length(cpd.enr[!cpd.enr %in% cpd.ass]) + length(Reduce(intersect,list(cpd.enr, cpd.ass))),
                                  area2= length(cpd.ass[!cpd.ass %in% cpd.enr]) + length(Reduce(intersect,list(cpd.enr, cpd.ass))),
                                  cross.area= length(Reduce(intersect,list(cpd.enr, cpd.ass))),
                                  category=c('Enrichment Analysis', 'Association Analysis'),
                                  scaled=F, 
                                  alpha = c(0.25,0.25),
                                  fill = c("#007bff", "#28a745"),

                                lwd = c(0, 0),
                                  lty = 'blank',
                                  
                                  cex = 3,
                                  fontface = "bold",
                                  fontfamily = "sans",

                                  # Set names
                                  cat.cex = 2,
                                  cat.fontface = "bold",
                                  cat.default.pos = "outer",
                                  cat.pos = c(-27, 27),
                                  cat.dist = c(0.055, 0.055),

                                  cat.fontfamily = "sans",
                                  
                                  filename = NULL,
                                  output=F,

                                  )
  #venn.plot.cpd <- grid.arrange(gTree(children=venn.plot.cpd), top="Overlap analysis for the CPDs\n(Enrichment Vs. Association)\n")
  
  return(venn.plot.cpd)
}

overlap.analysis.path <- function(enr.data, ass.data, org.data){
  
  org.data <- org.data[org.data$cpd != "",]
  #1st overlap analysis
  org.enr <- org.data[org.data$IDs_map %in% enr.data$KEGG_PATHWAY_ID,]
  org.enr.path <- org.enr$IDs_map
  org.ass <- org.data[org.data$IDs_map %in% ass.data$KEGG_PATHWAY_ID,]
  org.ass.path <- org.ass$IDs_map
  
  
  
  cpd.enr <- unlist( str_split(org.enr$cpd , ";") )
  cpd.ass <- unlist(str_split(org.ass$cpd , ";"))
  
  # venn dig. 

  venn.plot.path <- draw.pairwise.venn(area1= length(org.enr.path[!org.enr.path %in% org.ass.path]) + length(Reduce(intersect,list(org.enr.path, org.ass.path))),
                                      area2= length(org.ass.path[!org.ass.path %in% org.enr.path]) + length(Reduce(intersect,list(org.enr.path, org.ass.path))),
                                      cross.area= length(Reduce(intersect,list(org.enr.path, org.ass.path))),
                                      category=c('Enrichment Analysis', 'Association Analysis'),
                                      scaled=F,
                                      alpha = c(0.25,0.25),
                                      fill = c("#007bff", "#28a745"),

                                      lwd = c(0, 0),
                                      lty = 'blank',

                                      cex = 3,
                                      fontface = "bold",
                                      fontfamily = "sans",

                                      # Set names
                                      cat.cex = 2,
                                      cat.fontface = "bold",
                                      cat.default.pos = "outer",
                                      cat.pos = c(-27, 27),
                                      cat.dist = c(0.055, 0.055),

                                      cat.fontfamily = "sans",

                                      filename = NULL,
                                      output=F,

  )
  
  #venn.plot.path <- grid.arrange(gTree(children=venn.plot.path), top="Overlap analysis for the pathways\n(Enrichment Vs. Association)\n" )
  
  return(venn.plot.path)
  
  
}
