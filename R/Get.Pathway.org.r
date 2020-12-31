
# get org. pathways infromation from kegg using org. ID
# IDs: https://www.pnas.org/content/suppl/2008/09/11/0806162105.DCSupplemental/ST1_PDF.pdf


Get.Pathway.org <- function(org.id){
  

  pathway.url <- paste0("http://rest.kegg.jp/list/pathway/" , org.id)
  

  pathway.df <-  fread(pathway.url, fill=F , quote = "" , sep="\t", showProgress = T, header = F )
  
  colnames(pathway.df) <- c("IDs_org" , "Name")
  

  pathway.df$IDs_org <- str_remove(pathway.df$IDs_org , "path:")
  
  pathway.df$IDs_map <- str_replace(pathway.df$IDs_org ,org.id , "map")
  
  org.name <- read.csv( paste0("http://rest.kegg.jp/find/genome/",org.id), sep = ";", header = F)
  
  #org.name$V1
  
  pathway.df$org_acc_ids <- org.id
    
  return(pathway.df)
  
}