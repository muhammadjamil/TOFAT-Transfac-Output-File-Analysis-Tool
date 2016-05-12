options(rgl.useNULL=TRUE)
library(shinyRGL)
library(shinyjs)
shinyUI(
  fluidPage(
    titlePanel(
      strong((h1("TOFAT: TRANSFAC Output File Analysis Tool",style="color:Blue",align="center"))),windowTitle = "TOFAT")

    ##### left first panel###
    ,column(2,wellPanel(height=20,
                      fluidRow(
                        fileInput("CompleteReport",label = h6("CompleteReport")) 
                      ),
                      fluidRow(
                        fileInput("TargetReport",label = h6("Target Genes"))
                      ),
                      fluidRow(fileInput("TFReport",label = h6("TranscriptionFactor Report"))),
                      fluidRow(fileInput("FuncAnalysis",label=h6("Functional Analysis Report"))),
                      fluidRow(uiOutput("byAnno2"))
                      

      
    ))

    ### Plot Panel
    ,column(10,tabsetPanel(
      tabPanel("TranscriptionFactor Report",fluidRow(
        #column(1,
        #       downloadButton("SaveTable", "save")),
        tableOutput("limmaTable")
      )),
      tabPanel("Target Genes",fluidRow(
        #column(1,
        #       downloadButton("SaveTable", "save")),
        column(4,tableOutput("CReport")),
        column(6,
               fluidRow(tableOutput("CDetails")),
               column(5,fluidRow(plotOutput("line")))
      ))),
      tabPanel("Gene Ontology",
               tabsetPanel(
               tabPanel("Molecular Function",tableOutput("GOMF")),
               tabPanel("Biological Process",tableOutput("GOBP")),
               tabPanel("Cellular Component",tableOutput("GOCC")))),
      tabPanel("Disease or Tumors",
               tabsetPanel(
                 tabPanel("Disease",tableOutput("Dis")),
                 tabPanel("Tumors",tableOutput("Tum"))
               )),
      tabPanel("Location",
               tabsetPanel(
                 tabPanel("Organ or Tissue",tableOutput("OGT")),
                 tabPanel("Cells",tableOutput("CL"))
               )),
      tabPanel("Pathways",tableOutput("pways"))

     
    )
    
    )
    ))
