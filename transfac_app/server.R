library(ggplot2)
options(shiny.maxRequestSize=Inf)

shinyServer(function(input, output,session) {
  dataSet = reactive({
    inFile = input$TFReport
    inFile2 = input$CompleteReport
    inFile3 = input$TargetReport
    list(inFile=inFile,inFile2=inFile2,inFile3=inFile3)
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

})