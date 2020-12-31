fisher.test.enr <- function(found.list2, found.list1, list2, list1, method = "BH"){
  
  
  not.found.list2 <- list2 - found.list2
  not.found.list1 <- list1 - found.list1
  
  test.matrix <- matrix(c(found.list2, found.list1, not.found.list2, not.found.list1), 
                        nrow = 2)
  
  test.result <- stats::fisher.test(test.matrix)

  odds.ratio <- test.result$estimate
  p.value <- test.result$p.value
  
  p.adj <- stats::p.adjust(p.value, method = method)

  test.out <- data.frame(odds.ratio = odds.ratio, p.value = p.value, p.adj = p.adj)
  
  return(test.out)
  
}