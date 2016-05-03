options(rgl.useNULL=TRUE)
library(shinyRGL)
library(shinyjs)
shinyUI(
  fluidPage(
    titlePanel(
      strong((h1("TRANSFAC Analysis",style="color:Blue",align="center"))),windowTitle = "UKB TRANSFAC")

    ##### left first panel###
    ,column(2,wellPanel(height=20,
                      fluidRow(
                        fileInput("CompleteReport",label = h6("CompleteReport")) 
                      ),
                      fluidRow(
                        fileInput("TargetReport",label = h6("Target Genes"))
                      ),
                      fluidRow(fileInput("TFReport",label = h6("TranscriptionFactor Report"))),
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
      )))

     
    )
    
    )
    ))