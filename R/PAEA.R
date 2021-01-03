#Cran
library(shiny)
library(shinycustomloader)
library(shinyFiles)
library(rhandsontable)
library(shinythemes)
library(shinyjs)
library(stringr)
library(ggplot2)
library(dplyr)
library(ggnewscale)
library(data.table)
library(qdapRegex)
library(VennDiagram)
library(reshape2)
library(igraph)
library(networkD3)
library(htmlwidgets)
#github
library(shinyDirectoryInput)
#base
library(stats)
library(parallel)


source("Get.Pathway.org.r")
source("Get.cpd.fast.r")
source("Get.cpd.r")
source("Run.PEA.r")
source("Fold.Enrich.r")
source("plot.PEA.r")
source("fisher.test.r")
source("Run.fisher.r")
source("Get.associate1.1.r")
source("Do.associate.r")
source("overlap.analysis.r")
source("pathway.comp.shared.r")
source("relatedpathway.disease.r")
source("plot.network.r")
source("get.enzyme.reaction.r")
source("sankey.cpd.plot.r")


run.shiny.paea <- function(){
  app = shinyApp(

    ui = fluidPage(

      theme = shinytheme("sandstone"),

      sidebarLayout(
        sidebarPanel(width = 3,

                     shinyjs::useShinyjs(),

                     conditionalPanel(condition = "input.inTabset == 'Input_Summary_Table' | input.inTabset == 'Retrieved_cpd_Table' | input.inTabset == 'Enrichment_Analysis_Result' | input.inTabset == 'Association_Analysis_Result'  ",

                                      textInput("org_id","Enter the Species KEGG ID:",placeholder = "e,g. 'hsa' for human"),
                                      actionButton("IDs_Guide", "Open Species IDs Guide", style='padding:10px; font-size:100%', width = 200, icon = icon("file-pdf")),

                                      h3("OR"),
                                      fileInput("species_cpd", "Enter the Species KEGG CPD", multiple=FALSE, accept=c(".txt",".csv"),),


                                      fileInput('upload_cpd', 'Choose The CPD File', multiple=FALSE,
                                                accept=c(".txt",".csv"),),

                                      selectInput("p_adj_methods", "Adjusted p-value method:",
                                                  selected = "none",
                                                  c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                                                    "fdr", "none"),
                                                  multiple = FALSE),

                                      sliderInput("adj_p_en", "Adjusted p-value (enrichment analysis)", 0, 1, 0.01, value = 0.05),
                                      sliderInput("p_as", "p-value (association analysis)", 0, 1, 0.01, value = 0.01),
                                      sliderInput("N_cutoff", "Only use metabolite sets containing at least", 0, 10, 1, value = 2),



                                      actionButton("Run_en", "enrichment analysis", style='padding:10px; font-size:100%; margin-top:0.5em', width = 200, icon = icon("spinner")),
                                      actionButton("Run_as", "association analysis", style='padding:10px; font-size:100%; margin-top:0.5em', width = 200, icon = icon("spinner")),
                                      actionButton("download_session", "Save The Result", style='padding:10px; font-size:100%; margin-top:0.5em', width = 200, icon = icon("save")),
                                      actionButton("restart_session", "Restart The Session", class = "btn-success", style='padding:10px; font-size:100%; margin-top:0.5em', width = 200, icon = icon("refresh")),

                     ),
                     conditionalPanel(condition = "input.inTabset == 'save_sesults'",

                                      h2("Saving Results Tab"),


                                      sliderInput("width", "Plot Width (in): ", 1, 100, 1, value = 12),
                                      sliderInput("height", "Plot Height (in)", 1, 100, 1, value = 10),
                                      sliderInput("dpi", "Plot dpi", 50, 1000, 50, value = 300),


                                      textInput("save_result_by","Save the result by:",placeholder = "e,g. 'Sample1"),

                                      # shinyDirButton('folder_action', 'Save Result In',
                                      #                'Please select a folder',
                                      #                FALSE,
                                      #                style='padding:10px; font-size:100%; margin-top:0.5em',
                                      #                width = 200, icon = icon("download")),

                                      #directoryInput('folder_action', label = 'Please select a directory',value = '~' ),

                                      actionButton("save_now", "Download",
                                                   style='padding:10px; font-size:100%; margin-top:0.5em',
                                                   width = 200, icon = icon("download"))),


                     conditionalPanel(condition = "input.inTabset == 'welcome'",
                                      img(src = "https://raw.githubusercontent.com/AliYoussef96/PAEA/main/logos/paealogo.png" ,
                                          align = "center", style="width:200px;height:200px; margin-top:0.5em; block; margin-left: auto; margin-right: auto"),
                                      actionButton("welcome_start", "Start PAEA",
                                                   style='padding:10px; font-size:150%; margin-top:0.5em',
                                                   width = 200, icon = icon("door-open"))
                     ),


                     conditionalPanel(condition = "input.inTabset == 'additionalanalysis1' ",

                                      h1("Additional Analysis 1 (This may take longer time)"),

                                      # checkboxInput("enzymes", label = "Get Related Enzymes", value = FALSE),
                                      # checkboxInput("reaction", label = "Get Related Eeactions", value = FALSE),

                                      actionButton("additionalanalysis_start1", "Start",
                                                   style='padding:10px; font-size:100%; margin-top:0.5em', width = 200, icon = icon("tasks")),


                                      radioButtons("add_on_plots1", "Type of Plot:",
                                                   c(
                                                     "Comparative analysis for the P-values" = "pvalues",
                                                     "Comparative analysis for the Fold-values" = "foldvalues",
                                                     "Related Diseases" = "diseasescheckbox",
                                                     "Related Pathways" = "pathwayscheckbox")
                                      ),


                                      actionButton("nextnext", "Next", style='padding:10px; font-size:100%; margin-top:0.5em', width = 200, icon = icon("forward"))

                     ),

                     conditionalPanel(condition = "input.inTabset == 'additionalanalysis2' ",

                                      h1("Additional Analysis 2 (This may take longer time)"),

                                      h2("Parameters For The Network Plots"),

                                      sliderInput("charger", "The strength of the node repulsion", -1000, 0, -50, value = -100),

                                      sliderInput("linkDistance", "The distance between the links in pixels", 0, 100, 10, value = 20),

                                      h2("Parameters For The Sankey Plot"),

                                      sliderInput("toppathway", "Top N Pathways (P-adj.)", 0, 50, 1, value = 5),

                                      sliderInput("atleastcpd", "The number of pathways the CPD should have a role in;", 0, 10, 1, value = 1),

                                      actionButton("additionalanalysis_start2", "Start",
                                                   style='padding:10px; font-size:100%; margin-top:0.5em', width = 200, icon = icon("tasks")),


                                      radioButtons("add_on_plots2", "Type of Plot:",
                                                   c("Related Diseases Network" = "networkdiseasescheckbox",
                                                     "Related Pathways Network" = "networkpathwayscheckbox",
                                                     "Related CPDs Network" = "sankeycpdplotcheckbox")
                                      ),


                                      actionButton("backback", "Back", style='padding:10px; font-size:100%; margin-top:0.5em', width = 200, icon = icon("backward"))

                     ),



        ) ,

        mainPanel(
          img(src = "https://www.57357.org/app/uploads/2019/12/logo-2.png" , align = "right", style="width:150px"),
          img(src = "https://www.57357.org/app/uploads/2019/12/afnci.jpg" , align = "right", style="width:150px"),

          img(src = "https://www.57357.org/app/uploads/2019/12/ECN-Canada.jpg" , align = "right", style="width:150px"),
          img(src = "https://www.57357.org/app/uploads/2019/12/2-1.jpg" , align = "right", style="width:150px"),


          tabsetPanel(id = "inTabset",



                      tabPanel("Welcome To PAEA",
                               value = "welcome",
                               h1("Welcome to PAEA (Pathway Association Enrichment Analysis)"),
                               h2("Please Cite: ")),

                      tabPanel("Input Summary", withLoader(dataTableOutput("Input_Summary_Table"),
                                                           type="html", loader="loader2"), value = "Input_Summary_Table"),

                      tabPanel("Retrieved KEGG CPD", withLoader(dataTableOutput("Retrieved_cpd_Table"),
                                                                type="html", loader="loader2"), value = "Retrieved_cpd_Table"),


                      tabPanel("Enrichment Analysis Result",
                               plotOutput("Enrichment_Analysis_Result1"),
                               value = "Enrichment_Analysis_Result") ,


                      tabPanel("Association Analysis Result",
                               plotOutput("Association_Analysis_Result"),
                               value = "Association_Analysis_Result"),

                      tabPanel("Additional Analysis 1", h1(""),

                               div(plotOutput("additionalanalysisplots1",  height = "100%",
                                              width = "100%" ), align = "center"),


                               value = "additionalanalysis1") ,

                      tabPanel("Additional Analysis 2", h1(""),

                               forceNetworkOutput("additionalanalysisplots2",  height = "1000px",
                                                  width = "100%" ),
                               h1(""),
                               sankeyNetworkOutput("additionalanalysisplots2.2",  height = "1000px",
                                                   width = "100%" ),

                               value = "additionalanalysis2"),


                      tabPanel("Save Results",
                               value = "save_sesults")

          )
        )

      )
    ),


    server = function(input, output , session ) {
      options(shiny.maxRequestSize = 102400*1024^2)

      shinyjs::disable("Run_as")
      shinyjs::disable("download_session")
      shinyjs::disable("save_now")
      #shinyjs::disable("folder_action")
      shinyjs::disable("Run_en")
      shinyjs::disable("additionalanalysis_start1")
      shinyjs::disable("additionalanalysis_start2")



      #welcome tab
      observe({
        if(input$welcome_start[1] > 0){
          observeEvent(input$welcome_start,{

            updateTabsetPanel(session, "inTabset",
                              selected = 'Input_Summary_Table'
            )})

        }
      })

      ## slider
      sliderValues <- reactive({

        data.frame(
          Name = c("adj_p_en",
                   "p_as",
                   "N_cutoff",
                   "width",
                   "height",
                   "dpi",
                   "charger",
                   "linkDistance",
                   "toppathway",
                   "atleastcpd"),

          Value = as.numeric(c(input$adj_p_en,
                               input$p_as,
                               input$N_cutoff,
                               input$width,
                               input$height,
                               input$dpi,
                               input$charger,
                               input$linkDistance,
                               input$toppathway,
                               input$atleastcpd)),

          stringsAsFactors = FALSE)


      })

      ## restart
      observe({
        if(input$restart_session > 0){
          session$reload()
        }

      })


      ## org IDs input

      pathways.df <- NULL
      pathway.information.found <- NULL

      observe({
        if(input$org_id != ""){
          if(nchar(input$org_id) == 3){

            showModal(modalDialog("Retrieving The Species Information", footer=NULL))

            pathways.df <<- Get.Pathway.org(input$org_id)

            removeModal()


            output$Input_Summary_Table <- renderDataTable({

              pathways.df })

            shinyjs::disable("org_id")
            shinyjs::disable("species_cpd")

            pathway.information.found <<- "NO"
          }}

      })


      ## read CPD infor for the org. from csv

      pathway.cpd.kegg <- NULL

      observe({
        if(!is.null(input$species_cpd)){

          observeEvent(input$species_cpd,{

            updateTabsetPanel(session, "inTabset",
                              selected = 'Retrieved_cpd_Table'
            )})

          pathway.cpd.kegg <<- read.csv(input$species_cpd$datapath)

          output$Retrieved_cpd_Table <- renderDataTable({


            pathway.cpd.kegg})

          shinyjs::disable("org_id")
          shinyjs::disable("species_cpd")

          pathway.information.found <<- "Yes"
        }

      })

      ## upload the cpd for analysis and retrive information about pathways' cpd of the org. if not found
      input.cpd.user <- NULL
      observe({
        if(!is.null(input$upload_cpd)){

          shinyjs::enable("Run_en")


          if(is.null(pathway.cpd.kegg)){
            showModal(modalDialog("Retrieving The CPD Information", footer=NULL))


            pathway.cpd.kegg <<- Get.cpd.fast(pathways.df)
            pathway.cpd.kegg <<- pathway.cpd.kegg[!pathway.cpd.kegg$cpd == "",]


            input.cpd.user <<- read.csv(input$upload_cpd$datapath, header= F)
            input.cpd.user <<- as.character(input.cpd.user[,c(1)])

            observeEvent(input$upload_cpd,{

              updateTabsetPanel(session, "inTabset",
                                selected = 'Retrieved_cpd_Table'
              )})

            output$Retrieved_cpd_Table <- renderDataTable({

              removeModal()

              pathway.cpd.kegg})


          }else{
            input.cpd.user <<- read.csv(input$upload_cpd$datapath, header= F)
            input.cpd.user <<- as.character(input.cpd.user[,c(1)])
          }
        }
      })

      ## RUN normal Enrichment Analysis
      PEA.normal.result <- NULL
      PEA.normal.result.plot <- NULL
      observe({

        if( input$Run_en[1] > 0){

          showModal(modalDialog("Runing Enrichment Analysis...", footer=NULL))

          observeEvent(input$Run_en,{

            updateTabsetPanel(session, "inTabset",
                              selected = 'Enrichment_Analysis_Result'
            )})

          PEA.normal.result <<- Run.PEA(cpd.search.input = input.cpd.user,
                                        merged.df = pathway.cpd.kegg ,
                                        method = input$p_adj_methods,
                                        N = as.numeric(input$N_cutoff))

          PEA.normal.result <<- PEA.normal.result[PEA.normal.result$p_adj < as.numeric(input$adj_p_en),]

          PEA.normal.result <<- Fold.Enrich(PEA.normal.result)

          PEA.normal.result.plot <<- Plot.PEA(PEA.normal.result)


          output$Enrichment_Analysis_Result1 <- renderPlot({
            PEA.normal.result.plot
          }, height = 800, width = 1000 ,  )

          removeModal()
          shinyjs::enable("Run_as")

        }
      })

      ## association analysis

      pathway.cpd.kegg.do <- NULL
      association.CPD.IDs <- NULL
      PEA.association.result <- NULL
      PEA.association.result.plot <- NULL



      observe({
        if(input$Run_as[1] > 0){
          showModal(modalDialog("Runing Association Analysis...", footer=NULL))

          observeEvent(input$Run_en,{

            updateTabsetPanel(session, "inTabset",
                              selected = 'Association_Analysis_Result'
            )})

          pathway.cpd.kegg.do <<- pathway.cpd.kegg[pathway.cpd.kegg$IDs_map %in% PEA.normal.result$KEGG_PATHWAY_ID,]

          association.CPD.IDs <<- Do.associate(pathway.df.merged = pathway.cpd.kegg.do,
                                               input.ids = input.cpd.user,
                                               alpha = as.numeric(input$p_as))

          PEA.association.result <<- Run.PEA(cpd.search.input = association.CPD.IDs,
                                             merged.df = pathway.cpd.kegg ,
                                             method = input$p_adj_methods,
                                             N = as.numeric(input$N_cutoff))

          PEA.association.result <<- Fold.Enrich(PEA.association.result)


          PEA.association.result.plot <<- Plot.PEA(PEA.association.result)


          output$Association_Analysis_Result <- renderPlot({
            PEA.association.result.plot
          }, height = 800, width = 1000 )






          removeModal()
          shinyjs::enable("download_session")
          shinyjs::enable("additionalanalysis_start1")
          shinyjs::enable("additionalanalysis_start2")


        }
      })

      ## additional analysis1

      CPDs.plot <- NULL
      pathways.plot <- NULL
      foldvalues.plot <- NULL
      foldvalues.plot <- NULL
      pathway.related.plot <- NULL
      disease.related.plot <- NULL


      association.map.sig <- NULL
      related.pathways.disease <- NULL
      # enzymes.related <- NULL
      # reactions.related <- NULL

      observe({

        if(input$additionalanalysis_start1[1] >0 ){
          showModal(modalDialog("Loading...", footer=NULL))

          CPDs.plot <<- overlap.analysis.cpd(enr.data = PEA.normal.result, ass.data = PEA.association.result, org.data = pathway.cpd.kegg)
          pathways.plot <<- overlap.analysis.path(enr.data = PEA.normal.result, ass.data = PEA.association.result, org.data = pathway.cpd.kegg)

          pvalues.plot <<- pathway.comp.shared.pvalue( enr.result =  PEA.normal.result, ass.result =  PEA.association.result, org.data =  pathway.cpd.kegg)
          foldvalues.plot <<- pathway.comp.shared.fold( enr.result =  PEA.normal.result, ass.result =  PEA.association.result, org.data =  pathway.cpd.kegg)


          association.map.sig <<- PEA.association.result[PEA.association.result$p_adj < input$adj_p_en,]$KEGG_PATHWAY_ID

          related.pathways.disease <<- get.info.relatedpathway.disease(association.map.sig, pathway.cpd.kegg)

          pathway.related.plot <<- plot.pathway.related.bar(related.pathways.disease)
          disease.related.plot <<- plot.disease.bar(related.pathways.disease)

          # if(input$enzymes == TRUE){
          #
          #   pathway.sig <<- PEA.association.result[PEA.association.result$p_adj < input$adj_p_en,]
          #
          #   enzymes.related <<- get.enzymes(pathway.sig, pathway.cpd.kegg)
          #
          # }
          #
          # if(input$reaction == TRUE){
          #
          #   pathway.sig <<- PEA.association.result[PEA.association.result$p_adj < input$adj_p_en,]
          #
          #   reactions.related <<- get.reaction(pathway.sig, pathway.cpd.kegg)
          #
          # }

          output$additionalanalysisplots1 <- renderPlot({ pvalues.plot },height = 800, width = 1200)


          removeModal()
        }
      })


      observe({

        if (input$additionalanalysis_start1[1] >0 ){

          if(input$add_on_plots1 == "pvalues"){

            output$additionalanalysisplots1 <- renderPlot({ pvalues.plot },height = 800, width = 1200)
          }

          if(input$add_on_plots1 == "foldvalues"){

            output$additionalanalysisplots1 <- renderPlot({ foldvalues.plot },height = 800, width = 1200)
          }

          if(input$add_on_plots1 == "diseasescheckbox"){

            output$additionalanalysisplots1 <- renderPlot({ disease.related.plot },height = 800, width = 1200)
          }

          if(input$add_on_plots1 == "pathwayscheckbox"){

            output$additionalanalysisplots1 <- renderPlot({ pathway.related.plot },height = 800, width = 1200)
          }


        }
      })

      ## additional analysis2 only on the sig. pathways


      pathway.network.plot <- NULL
      disease.network.plot <- NULL
      snakey.plot <- NULL


      observe({

        if(input$additionalanalysis_start2[1] > 0){



          if(input$add_on_plots2 == "networkdiseasescheckbox"){
            if(!is.null(related.pathways.disease)){

              disease.network.plot <<- plot.network.disease(related.pathways.disease, charge = input$charger, linkDistance = input$linkDistance)

              output$additionalanalysisplots2 <- renderForceNetwork( disease.network.plot )
            }}


          if(input$add_on_plots2 == "networkpathwayscheckbox"){
            if(!is.null(related.pathways.disease) ){

              pathway.network.plot <<- plot.network.pathway(related.pathways.disease, charge = input$charger, linkDistance = input$linkDistance)

              output$additionalanalysisplots2 <- renderForceNetwork( pathway.network.plot )
            }}

          if(input$add_on_plots2 == "sankeycpdplotcheckbox"){
            if(!is.null(related.pathways.disease) ){
              output$additionalanalysisplots2 <- renderForceNetwork( NULL )

              pathway.sig <<- PEA.association.result[PEA.association.result$p_adj < input$adj_p_en,]


              tryCatch(
                expr = {

                  snakey.plot <<- sanky.cpd.plot(pathway.sig = pathway.sig, org.df = pathway.cpd.kegg,
                                                 top.n = input$toppathway,
                                                 cpd.atleast = input$atleastcpd)
                  output$additionalanalysisplots2.2 <- renderSankeyNetwork( snakey.plot )


                },
                warning = function(w){
                  output$additionalanalysisplots2.2 <- renderSankeyNetwork( NULL )

                },
                error = function(e){
                  output$additionalanalysisplots2.2 <- renderSankeyNetwork( NULL )

                }
              )



            }
          }
        }

      })
      #next

      observe({
        if(input$nextnext[1] > 0){
          observeEvent(input$nextnext,{

            updateTabsetPanel(session, "inTabset",
                              selected = "additionalanalysis2"
            )})
        }
      })

      #back

      observe({
        if(input$backback[1] > 0){
          observeEvent(input$backback,{

            updateTabsetPanel(session, "inTabset",
                              selected = "additionalanalysis1"
            )})
        }
      })


      observe({
        if(input$IDs_Guide){
          browseURL("https://www.pnas.org/content/suppl/2008/09/11/0806162105.DCSupplemental/ST1_PDF.pdf")
        }
      })
      ## save result


      project_folder <- NULL


      observe({
        if(input$download_session[1] > 0){
          shinyjs::enable("save_now")
          #shinyjs::enable("folder_action")

          project_folder <<- choose.dir()

          observeEvent(input$download_session,{

            updateTabsetPanel(session, "inTabset",
                              selected = "save_sesults"
            )})



        }
      })


      ## Folder Selection


      # observe({
      #   #volumes <- getVolumes()
      #   #root <- getVolumes()()
      #   #shinyDirChoose(input, 'folder_action', roots=volumes() , session = session )
      #
      #   if (input$folder_action[1] > 0) {
      #     # condition prevents handler execution on initial app launch
      #
      #     # launch the directory selection dialog with initial path read from the widget
      #     project_folder <<- readDirectoryInput(session, 'folder_action')
      #
      #     # update the widget value
      #     updateDirectoryInput(session, 'folder_action', value = project_folder)
      #   }
      #
      #
      # })


      observe({
        if(input$save_now[1] > 0){

          showModal(modalDialog("Saving...", footer=NULL))

          save.result.path <- project_folder


          ## plots

          ggsave(filename = paste0(save.result.path, "/" , input$save_result_by, "_Enrichmentplot.tiff"), plot = PEA.normal.result.plot ,
                 dpi = input$dpi, width = input$width,
                 height = input$height, limitsize = FALSE)

          ggsave(filename = paste0(save.result.path, "/" , input$save_result_by, "_Associationplot.tiff"), plot = PEA.association.result.plot ,
                 dpi = input$dpi, width = input$width,
                 height = input$height, limitsize = FALSE)


          ggsave(filename = paste0(save.result.path, "/" , input$save_result_by, "_Overlap analysis for the CPDs.tiff"), plot = CPDs.plot ,
                 dpi = input$dpi, width = input$width,
                 height = input$height, limitsize = FALSE)



          ggsave(filename = paste0(save.result.path, "/" , input$save_result_by, "_Overlap analysis for the pathways.tiff"), plot = pathways.plot ,
                 dpi = input$dpi, width = input$width,
                 height = input$height, limitsize = FALSE)


          ggsave(filename = paste0(save.result.path, "/" , input$save_result_by, "_Comparative analysis for the P-values.tiff"), plot = pvalues.plot ,
                 dpi = input$dpi, width = input$width,
                 height = input$height, limitsize = FALSE)


          ggsave(filename = paste0(save.result.path, "/" , input$save_result_by, "_Comparative analysis for the Fold-values.tiff"), plot = foldvalues.plot ,
                 dpi = input$dpi, width = input$width,
                 height = input$height, limitsize = FALSE)


          if(!is.null(pathway.network.plot)){
            saveWidget(pathway.network.plot, paste0(save.result.path, "/" , input$save_result_by, "_PathWayNetwork.html"))
          }
          if(!is.null(disease.network.plot)){
            saveWidget(disease.network.plot, paste0(save.result.path, "/" , input$save_result_by, "_DiseaseNetwork.html"))
          }

          if(!is.null(snakey.plot)){
            saveWidget(snakey.plot, paste0(save.result.path, "/" , input$save_result_by, "_CDPsSnaky.html"))
          }

          if(!is.null(pathway.related.plot)){
            ggsave(filename = paste0(save.result.path, "/" ,
                                     input$save_result_by, "_Related Pathways.tiff"), plot = pathway.related.plot ,
                   dpi = input$dpi, width = input$width,
                   height = input$height, limitsize = FALSE)
          }

          if(!is.null(disease.related.plot)){
            ggsave(filename = paste0(save.result.path, "/" ,
                                     input$save_result_by, "_Related Diseases.tiff"), plot = disease.related.plot ,
                   dpi = input$dpi, width = input$width,
                   height = input$height, limitsize = FALSE)
          }


          #### csv
          write.csv(PEA.normal.result, paste0(save.result.path, "/" , input$save_result_by, "_Enrichmentplot.csv"))

          write.csv(PEA.association.result, paste0(save.result.path, "/" , input$save_result_by, "_Associationplot.csv"))

          if (pathway.information.found == "No"){
            write.csv(pathway.cpd.kegg , paste0(save.result.path, "/" , input$save_result_by, "_Host_info.csv") )
          }



          if(!is.null(related.pathways.disease)){


            write.csv( related.pathways.disease, paste0(save.result.path, "/" , input$save_result_by, "_Related Pathways and Diseases.csv") )

          }

          # if(!is.null(enzymes.related)){
          #   write.csv(enzymes.related, paste0(save.result.path, "/" , input$save_result_by, "Related Enzymes.csv") )
          #
          # }
          #
          #
          # if(!is.null(reactions.related)){
          #   write.csv(reactions.related, paste0(save.result.path, "/" , input$save_result_by, "Related Reactions.csv") )
          #
          # }

          removeModal()

        }
      })

    }
  )
  runApp(app)
}

run.shiny.paea()
