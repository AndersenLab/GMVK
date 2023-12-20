library(shiny)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

orthos <- read_tsv("data/masterOrthoDB_wAlias.tsv")
#orthos <- read_tsv("/projects/b1059/projects/Nicolas/shiny_gmvk/data/masterOrthoDB_wAlias.tsv")
tair <-read_tsv("data/master_DB_V2.tsv")
#tair <- read_tsv("/projects/b1059/projects/Nicolas/shiny_gmvk/data/master_DB.tsv")
source("helper.R")

ui <- fluidPage(
  
  titlePanel(
    img(src = "gmvk_logo.png", height = 250, width = 400)
  ),
  
  sidebarLayout(
    
    sidebarPanel(
      
      
      textInput(input="geneid",
                label=p("Enter ",em("C. elegans"), "gene name(s) or ID(s):"),
                value=""
      ),
      
      checkboxGroupInput(input="spp", 
                         label = "Select species:", 
                         choices = list("C. elegans"="N2",
                                        "C. briggsae"="QX1410",
                                        "C. tropicalis"="NIC58")
                         ),
      
      h4("Additional Options"),
      
      radioButtons(input="ort",
                   label = "Orthogroup:",
                   choices = list("All orthologs"="orth",
                                  "Target gene only"="targ")
      ),
      
      
      radioButtons(input="mode",
                   label = "Alternative splice isoforms:",
                   choices = list("All isoforms"="iso",
                                  "Longest isoform only"="lon")
                   ),
      
      radioButtons(input="cut",
                   label = "Display options:",
                   choices = list("Joint"="joint",
                                  "Separate by species"="spp",
                                  "Sepearate by orthology"="gsep",
                                  "Separate by species and orthology"="gspp")
                   ),
                   
      radioButtons(input="lab",
                   label = "Axis freedom:",
                   choices = list("Free Y"="free_y",
                                  "Free X & Y"="free")
                  ),
      
      sliderInput(input="hslider", 
                  label = "Plot height:", 
                  min = 300, 
                  max = 1000,
                  value = 700,
                  step=100
                  ),
      
      sliderInput(input="cexslider", 
                  label = "Label size:", 
                  min = 1, 
                  max = 10, 
                  value = 7,
                  step=0.5
                  ),
      sliderInput(input="legendy", 
                  label = "Legend horizontal position:", 
                  min = 0, 
                  max = 1, 
                  value = 0.95,
                  step=0.05
                  ),
      sliderInput(input="legendx", 
                  label = "Legend vertical position:", 
                  min = 0, 
                  max = 1, 
                  value = 0.15,
                  step=0.05
                  ),
      checkboxGroupInput(input="anno", 
                         label = "Variant Annotations:", 
                         choices = list("Enabled"="anno")
      )
    ),
    
    mainPanel(
      h3("Orthology:"),
      tableOutput("view"),
      
      h3("Gene Models:"),
      #tableOutput("genes"),
      uiOutput("plot.ui")
      
      
    )
  )
)

server <- function(input,output) {
  
  datasetOrtho <- reactive({
    req(input$geneid)
    req(input$spp)
    
    gene_id <- gsub("[[:space:]]","",input$geneid)
    
    orthos %>% 
      dplyr::filter(if_any(all_of(c("seqname", "WB_id", "WB_alias")), ~grepl(paste0("\\<",paste(unlist(str_split(gene_id,",")),collapse = "\\>|\\<"),"\\>"), .))) %>%
      dplyr::select(c(input$spp,WB_id,WB_alias,seqname)) %>%
      dplyr::select(-seqname,-WB_id)
    
  })
  
  output$view <- renderTable({
    datasetOrtho()
  })
  
  plotht <- reactive({
    req(input$hslider)
    as.numeric(input$hslider)
  })
  
  plotcex <- reactive({
    req(input$cexslider)
    as.numeric(input$cexslider)
  })
  
  plotlx <- reactive({
    req(input$legendx)
    as.numeric(input$legendx)
  })
  
  plotly <- reactive({
    req(input$legendy)
    as.numeric(input$legendy)
  })

  output$models <- renderPlot({
    req(input$geneid)
    req(input$spp)
    
    gene_id <- gsub("[[:space:]]","",input$geneid)
    #gene_id <- "ben-1,isw-1"
    spp<- "QX1410"
    subset_genes <- orthos %>% 
      dplyr::filter(if_any(all_of(c("seqname", "WB_id", "WB_alias")), ~grepl(paste0("\\<",paste(unlist(str_split(gene_id,",")),collapse = "\\>|\\<"),"\\>"), .))) %>%
      #dplyr::filter(if_any(all_of(c("seqname", "WB_id", "WB_alias")), ~grepl(paste0("\\<",paste(unlist(str_split(gene_id,",")),collapse = "\\>|\\<"),"\\>"), .))) %>%
      dplyr::select(c("Orthogroup", "WB_id","QX1410","NIC58")) %>%
      dplyr::group_by(Orthogroup) %>%
      dplyr::mutate(all=paste(c(unlist(strsplit(QX1410,", ")),unlist(strsplit(NIC58,", ")),unlist(strsplit(WB_id,", "))),collapse = ", ")) %>%
      tidyr::separate_rows(all,sep=", ") %>%
      dplyr::select(Orthogroup,all)
    
    ref <- data.frame(tag=c("N2","QX1410","NIC58"),colID=c("WB_id","QX1410","NIC58"),prefix=c("WBGene","QX1410","NIC58")) %>%
      dplyr::filter(tag %in% input$spp)
      #dplyr::filter(tag %in% spp)
    
    subset_genes2 <- subset_genes %>% dplyr::filter(grepl(paste(ref$prefix,collapse = "|"),all))
    
    plot_df <- tair %>%
      dplyr::left_join(subset_genes2,by=c("L1"="all")) %>%
      dplyr::filter(!is.na(Orthogroup))
    
    if (input$ort == "targ") {
      plot_df <- plot_df %>% dplyr::filter(grepl(gsub(",","|",gene_id),mAlias))
    }
    
    plotNormGeneFeatures(plot_df,input$mode,input$cut,input$lab,plotcex(),plotlx(),plotly())
  })
  
  output$plot.ui <- renderUI ({
    plotOutput("models", height=plotht())
  })
}

shinyApp(ui = ui,server=server)

