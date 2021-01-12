# PAEA (Pathway Association Enrichment Analysis)

<p align="center">
  <img width="250" height="250" src="https://raw.githubusercontent.com/AliYoussef96/PAEA/main/logos/paealogo.png">
</p>


**If you find PAEA useful, please cite our paper:**


# Description

PAEA (Pathway Association Enrichment Analysis), is a new method for metabolomics pathway enrichment analysis using a novel association algorithm to increase the initial input list with significant associates metabolites (compounds), to give a bigger picture on possible activated pathways. In this version KEGG (Kyoto Encyclopedia of Genes and Genomes) is used as background database, for enrishment and association analysis. All the species available on KEGG can be used by the user.

For more information about the PAEA implementation and how the algorithm work, please take a look on the original [manuscript]() publish at [journal name]().

# Availability and Installation 


**PAEA is available in two forms on R;**

## 1- Offline version as a Shiny app.

### Step 1. Install package dependencies 

To use PAEA 1.0.0, first install all package dependencies. Ensure that you have necessary system environment configured. R base with version > 3.6.1 is recommended.


```R
dependencies.pkg.cran <- c("shiny", "shinycustomloader", "shinyFiles", "rhandsontable", "shinythemes", "shinyjs", "ggplot2" , 
                           "dplyr", "ggnewscale", "data.table", "qdapRegex", "VennDiagram" , "reshape2", "igraph",
                           "networkD3", "htmlwidgets", "stringr", "ggraph")

install.packages(dependencies.pkg.cran)

```
### Step 2. 

2.1. Download the zip file of this repository 

![](https://raw.githubusercontent.com/AliYoussef96/PAEA/main/for_installation/1.jpg)

2.2. Extract the zip file in your computer at any directory

![](https://raw.githubusercontent.com/AliYoussef96/PAEA/main/for_installation/2.jpg)

2.3 Locate the folder named R 

![](https://raw.githubusercontent.com/AliYoussef96/PAEA/main/for_installation/3.jpg)

2.4. Open the R file named PAEA.r

![](https://raw.githubusercontent.com/AliYoussef96/PAEA/main/for_installation/4.jpg)

2.5. Run the whole script by pressing Ctrl+A to select all the script lines, and then press run (or press Ctrl+Enter)

![](https://raw.githubusercontent.com/AliYoussef96/PAEA/main/for_installation/5.jpg)

**The PAEA Shiny app should run locally, if everything is fine**

**See the Full Documentation for all the [details]() ; how to use the PAEA Shiny app. and how to save the results.**

## 2-	Online version as a Shiny app.

**PAEA can also run online at**

# Documentation Availability

For step by step documentation on how to use PAEA shiny app, follow this [link](https://github.com/AliYoussef96/PAEA/blob/main/Documentation.pdf). Also, a [tutorial video](https://drive.google.com/file/d/1OKhkrU4O7niL-UA25bD4dy4t1VU9UiBq/view) is available.



# Package information

- Version: **v1.0.0**
- License: GPL-3
- Encoding: UTF-8


# Contribution Guidelines

For bugs and suggestions, the most effective way is by raising an issue on the github issue tracker. Github allows you to classify your issues so that we know if it is a bug report, feature request or feedback to the authors.

**Email: Proteomicslab2017@gmail.com**
**Email: ali.anwar@57357.org**
**Email: aliali.mostafa99@gmail.com**
**Email: MohmedSoudy2009@gmail.com**

# Future plans

# Citation

![](https://www.57357.org/app/uploads/2019/12/logo-2.png)
