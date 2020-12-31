Get.associate1.1 <- function(pathway.metabolites, x.input, z.database){
  
  pathway.metabolites <- data.frame(pathway.metabolites)
  
  both.there <- pathway.metabolites[str_detect(pathway.metabolites[,1], fixed(x.input)) 
                                    & str_detect( pathway.metabolites[,1], fixed(z.database)),]
  
  
  x.only <- pathway.metabolites[str_detect(pathway.metabolites[,1],fixed(x.input)),]
  
  z.only <- pathway.metabolites[str_detect(pathway.metabolites[,1],fixed(z.database)),]
  
  #both.there <- data.frame(both.there)
  #x.only <- data.frame(x.only)
  #z.only <- data.frame(z.only)
  
  
  test.matrix <- matrix(c(length(both.there)
                          ,length(x.only),
                          length(z.only),
                          abs((length(both.there)+length(x.only)+length(z.only))-nrow(pathway.metabolites))
  ),
  nrow = 2)
  
  test.result <- stats::fisher.test(test.matrix)
  
  p.value <- test.result$p.value
  
  return(p.value)
  
}
