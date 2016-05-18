library(ggplot2)
options(shiny.maxRequestSize=Inf)
x = read.csv("BKL Accession and IlluminaIds.txt",sep="\t",header=TRUE)
shinyServer(function(input, output,session) {
  dataSet = reactive({
    inFile = input$TFReport
    inFile2 = input$CompleteReport
    inFile3 = input$TargetReport
    inFile4 = input$FuncAnalysis
    list(inFile=inFile,inFile2=inFile2,inFile3=inFile3,inFile4=inFile4)
  })
  TFVariable <- reactiveValues()
  TFVariable$TF = 0
  TFVariable$Genes = 0
  ########Select Transcription Factor Matrix for Target Genes###########
  output$byAnno2 <- renderUI({
    inFile = dataSet()$inFile
    if (is.null(inFile))
      return(NULL)
    else{
      annotation = read.csv(inFile$datapath,sep="\t",header=TRUE)
      selectInput(inputId="byAnno",label=h6("Select TF Matrix"),choices=as.character(annotation[,1]))
    }
  })
  output$limmaTable <- renderTable(
    if (is.null(dataSet()$inFile))
      return(NULL)
    else{

      annotation = read.csv(dataSet()$inFile$datapath,header=TRUE,sep="\t")
      TFVariable$TF = annotation
      annotation
    }
  )
  output$CReport <- renderTable({
    if (is.null(dataSet()$inFile2) || is.null(dataSet()$inFile3) || is.null(dataSet()$inFile))
      return(NULL)
    else{
      x = readLines(dataSet()$inFile2$datapath)
      proteinNamesPos = grep("Inspecting",x)
      proteinNames = lapply(1:length(proteinNamesPos),function(i) unlist(strsplit(x[proteinNamesPos[i]]," "))[6])
      res = rep(list(list()),length(proteinNamesPos))
      for (i in 1:(length(proteinNamesPos)-1)){
        proteinMap = x[(proteinNamesPos[i]+2) : (proteinNamesPos[i+1]-2)]
        res[[i]] = unlist(lapply(1:length(proteinMap),function(j) unlist(strsplit(proteinMap[j]," "))[2]))
      }
      proteinMap = x[(proteinNamesPos[length(res)]+2) : (length(x)-6)]
      res[[length(res)]] = unlist(lapply(1:length(proteinMap),function(j) unlist(strsplit(proteinMap[j]," "))[2]))
      names(res) = proteinNames
      ProteinDetails = read.csv(dataSet()$inFile3$datapath,sep="\t")
      TFDetails = TFVariable$TF
      #output$byAnno2 <- renderUI({
      #  inFile = dataSet()$inFile
      #    annotation = TFVariable$TF
      #    selectInput(inputId="byAnno",label=h6("Select TranscriptionFactor"),choices=annotation[,4])
      #})
        TFList = lapply(res, function(j) grep(unlist(strsplit(as.character(input$byAnno),"\\$"))[2], j))
        TargetedGenes = names(unlist(lapply(TFList,function(j) if(length(j) > 0) return(TRUE))))
        resultGenes = ProteinDetails[unlist(lapply(TargetedGenes,function(j) grep(j,ProteinDetails[,1]))),4]
        names(resultGenes) = TargetedGenes
        TFVariable$Genes = resultGenes
        res2 = cbind(resultGenes)
        res2 = data.frame(res2[order(as.numeric(as.character(res2[,1]))),])
        colnames(res2) = ("logFC")
        res2
        
        
        
    }
  })
  output$CDetails <- renderTable({
    if (is.null(dataSet()$inFile2) || is.null(dataSet()$inFile3) || is.null(dataSet()$inFile))
      return(NULL)
    else{
      TFDetails = TFVariable$TF
      transfac2illumina2 <- function(x){
        anno = read.csv("C:/Users/Osman/Desktop/Ahmer/Shiny Tutorial/transfac_app/enterzId and illumina.txt",sep="\t",header=TRUE)
        entrezIDs = unlist(strsplit(toString(x[1,4])," "))
        names(entrezIDs) = x[,1]
        res = matrix(nrow=0,ncol=2)
        for(i in 1:length(entrezIDs)){
          res = rbind(res,
                      cbind(rep(names(entrezIDs)[[i]],length(entrezIDs[[i]])),entrezIDs[[i]])
          )
        }
        res = data.frame(res)
        colnames(res) = c("Matrices","EntrezIDs")
        res = merge(res,anno,by.x="EntrezIDs",by.y="ENTREZ_GENE_ID")
        res = res[,c(2,1,3,4)]
        return (res)
      }
      transfac2illumina2(TFDetails[which(TFDetails == input$byAnno),])

    }
  })
  output$line = renderPlot(
    if (TFVariable$Genes == 0)
      return (NULL)
    else{
         boxplot(TFVariable$Genes,
         xlab=paste(length(TFVariable$Genes)," Targeted",sep=""),
         main = input$byAnno,pch=20,cex=0.2,
         ylab="LogRatio",xaxt="n")
      #qplot(x=c(0,1),y=TFVariable$Genes,geom=c("boxplot","jitter"),xlab="",ylab="logFC",main=input$byAnno)
      }
  )
  ###### Functional Analysis Tables Output########
  output$GOMF <- renderTable(
    if (is.null(input$FuncAnalysis))
      return(NULL)
    else{
      annotation = read.csv(input$FuncAnalysis$datapath,header=TRUE,sep="\t")
      res = subset(annotation,annotation[,1] == "GOMF")
      GeneList = c()
      for (i in 1:nrow(res)){
        GeneList[i] = paste(x[which(x$BIOBASE.accession %in% unlist(strsplit(toString(res$Genes[i])," "))),"HGNC"],collapse = " ")
      }
      res$Genes = GeneList
      res
    }
  )
  output$GOBP <- renderTable(
    if (is.null(input$FuncAnalysis))
      return(NULL)
    else{
      annotation = read.csv(input$FuncAnalysis$datapath,header=TRUE,sep="\t")
      res = subset(annotation,annotation[,1] == "GOBP")
      GeneList = c()
      for (i in 1:nrow(res)){
        GeneList[i] = paste(x[which(x$BIOBASE.accession %in% unlist(strsplit(toString(res$Genes[i])," "))),"HGNC"],collapse = " ")
      }
      res$Genes = GeneList
      res
    }
  )
  output$GOCC <- renderTable(
    if (is.null(input$FuncAnalysis))
      return(NULL)
    else{
      annotation = read.csv(input$FuncAnalysis$datapath,header=TRUE,sep="\t")
      res = subset(annotation,annotation[,1] == "GOCC")
      GeneList = c()
      for (i in 1:nrow(res)){
        GeneList[i] = paste(x[which(x$BIOBASE.accession %in% unlist(strsplit(toString(res$Genes[i])," "))),"HGNC"],collapse = " ")
      }
      res$Genes = GeneList
      res
    }
  )
  output$Dis <- renderTable(
    if (is.null(input$FuncAnalysis))
      return(NULL)
    else{
      annotation = read.csv(input$FuncAnalysis$datapath,header=TRUE,sep="\t")
      res = subset(annotation,annotation[,1] == "disease")
      GeneList = c()
      for (i in 1:nrow(res)){
        GeneList[i] = paste(x[which(x$BIOBASE.accession %in% unlist(strsplit(toString(res$Genes[i])," "))),"HGNC"],collapse = " ")
      }
      res$Genes = GeneList
      res
    }
  )
  output$Tum <- renderTable(
    if (is.null(input$FuncAnalysis))
      return(NULL)
    else{
      annotation = read.csv(input$FuncAnalysis$datapath,header=TRUE,sep="\t")
      res = subset(annotation,annotation[,1] == "tumor")
      GeneList = c()
      for (i in 1:nrow(res)){
        GeneList[i] = paste(x[which(x$BIOBASE.accession %in% unlist(strsplit(toString(res$Genes[i])," "))),"HGNC"],collapse = " ")
      }
      res$Genes = GeneList
      res
    }
  )
  output$OGT <- renderTable(
    if (is.null(input$FuncAnalysis))
      return(NULL)
    else{
      annotation = read.csv(input$FuncAnalysis$datapath,header=TRUE,sep="\t")
      res = subset(annotation,annotation[,1] == "tissue")
      GeneList = c()
      for (i in 1:nrow(res)){
        GeneList[i] = paste(x[which(x$BIOBASE.accession %in% unlist(strsplit(toString(res$Genes[i])," "))),"HGNC"],collapse = " ")
      }
      res$Genes = GeneList
      res
    }
  )
  output$CL <- renderTable(
    if (is.null(input$FuncAnalysis))
      return(NULL)
    else{
      annotation = read.csv(input$FuncAnalysis$datapath,header=TRUE,sep="\t")
      res = subset(annotation,annotation[,1] == "cell")
      GeneList = c()
      for (i in 1:nrow(res)){
        GeneList[i] = paste(x[which(x$BIOBASE.accession %in% unlist(strsplit(toString(res$Genes[i])," "))),"HGNC"],collapse = " ")
      }
      res$Genes = GeneList
      res
    }
  )
  output$pways <- renderTable(
    if (is.null(input$FuncAnalysis))
      return(NULL)
    else{
      annotation = read.csv(input$FuncAnalysis$datapath,header=TRUE,sep="\t")
      res = subset(annotation,annotation[,1] == "pathway")
      GeneList = c()
      for (i in 1:nrow(res)){
        GeneList[i] = paste(x[which(x$BIOBASE.accession %in% unlist(strsplit(toString(res$Genes[i])," "))),"HGNC"],collapse = " ")
      }
      res$Genes = GeneList
      res
    }
  )
  
})
