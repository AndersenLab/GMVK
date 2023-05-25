library(shiny)
library(readr)
library(dplyr)
library(tidyr)

orthos <- read_tsv("data/masterOrthoDB_wAlias.tsv")
tair <-read_tsv("data/master_DB.tsv")
source("helper.R")

ui <- fluidPage(
  
  titlePanel(
    img(src = "gmvk_logo.png", height = 250, width = 400)
  ),
  
  sidebarLayout(
    
    sidebarPanel(
      
      checkboxGroupInput(input="spp", 
                         label = p(em("Caenorhabditis")," Ortholog and Paralog finder:"), 
                         choices = list("C. elegans"="N2",
                                        "C. briggsae"="QX1410",
                                        "C. tropicalis"="NIC58")),
      
      textInput(input="geneid",
                label=p("Enter ",em("C. elegans"), "gene name or ID:"),
                value=""
      )
    ),
    
    mainPanel(
      h3("Orthology:"),
      tableOutput("view"),
      
      h3("Gene Models:"),
      #tableOutput("genes"),
      plotOutput("models",height=700)
      
      
    )
  )
)

server <- function(input,output) {
  
  datasetOrtho <- reactive({
    req(input$geneid)
    req(input$spp)
    
    orthos %>% 
      dplyr::select(c(input$spp,WB_id,WB_alias,seqname)) %>%
      dplyr::filter(if_any(all_of(c("seqname", "WB_id", "WB_alias")), ~grepl(input$geneid, .))) %>%
      dplyr::select(-seqname,-WB_id)
    
  })
  
  datasetTair <- reactive({

    req(input$geneid)
    req(input$spp)
    # ODB_filter <- orthos %>%
    #   dplyr::select(c(input$spp,WB_id,WB_alias,seqname)) %>%
    #   dplyr::select(-seqname,-WB_id)
    # 
    # cbrig <- ODB_filter$QX1410
    # ctrop <- ODB_filter$NIC58
    # cele <- ODB_filter$N2


    tair %>% dplyr::filter(mAlias==input$geneid)
  })

  
  output$view <- renderTable({
    datasetOrtho()
 })
  
  output$genes <- renderTable({
    datasetTair()
  })

  output$models <- renderPlot({
    req(input$geneid)
    req(input$spp)
    
    subset_genes <- orthos %>% 
      dplyr::select(c(input$spp,c("seqname", "WB_id", "WB_alias"))) %>%
      dplyr::filter(if_any(all_of(c("seqname", "WB_id", "WB_alias")), ~grepl(input$geneid, .)))
    
    filter <- c(unlist(strsplit(subset_genes$WB_id,", ")),unlist(strsplit(subset_genes$QX1410,", ")),unlist(strsplit(subset_genes$NIC58,", ")))
    
    plot_df <- tair %>% dplyr::filter(L1 %in% filter)
    
    plotNormGeneFeatures(plot_df)
  })
}

shinyApp(ui = ui,server=server)

