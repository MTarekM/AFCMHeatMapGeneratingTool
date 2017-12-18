# Copyright (C) 2016-2017 Mohammad Tarek Mansour - AFCM BioInformatics Lab

# This file is the server file of the Application

# ------------------------------------------------------------------------------------

library(shiny)
library(d3heatmap)
library(htmlwidgets)
library(edgeR)
library(dplyr)
library(png)


server <- shinyServer(function(input, output) {	
  
  # sample Matrix download
  output$downloadSample <- downloadHandler(
    filename <- function() {
      paste('SampleFile', '.csv', sep='')
    },
    content <- function(file) {
      file.copy("SampleFile.csv", file)
    },
    contentType = "text/csv"
  )
 
  #  user file upload
  datasetInput <- reactive({
    validate(
      need(input$filename != 0, "To generate a heatmap using AFCM Generation tool, please upload your *.CSV file") 
    )
    inFile <- input$filename
    if (is.null(inFile)) return(NULL)
    read.table(inFile$datapath, header= TRUE, sep=",", quote='"', row.names=1)
  })
 
  # filter AFCMstats table based on user cutoffs
  AFCMstats <- reactive({
    AFCMres_df <- AFCMstatse()
    datasetInputTwo <- datasetInput()
    
    if (input$bothPvalueFDR == "Pvalue"){
      pNAFCMres_df <- AFCMres_df[which(AFCMres_df$PValue<input$userPvalue),]
      pafcmno <- row.names(pNAFCMres_df)
      pFinalResult <- datasetInputTwo[pafcmno,]
    }
    else if(input$bothPvalueFDR == "FDR"){
      fNAFCMres_df <- AFCMres_df[which(AFCMres_df$FDR<input$userFDR),]
      fafcmno <- row.names(fNAFCMres_df)
      fFinalResult <- datasetInputTwo[fafcmno,]
    }
    else if(input$bothPvalueFDR == "both"){
      bNAFCMres_df <- AFCMres_df[which(AFCMres_df$FDR<input$userFDR & AFCMres_df$PValue<input$userPvalue),]
      bafcmno <- row.names(bNAFCMres_df)
      bFinalResult <- datasetInputTwo[bafcmno,]
    }
  })
  
  
  # heatmap height
  output$pixelation <- renderUI({
     
    AFCMstats()
    inputLines <- NROW(AFCMstats())
    if(inputLines >= 0 && inputLines <= 1000)
      d3heatmapOutput("heatmap", width = "100%", height = "1000px")
    else
      d3heatmapOutput("heatmap", width = "100%", height = "10000px")
    
  })	
  
  
  # data Log transformation 
  log2_datasetInput <- reactive({
    AFCMstats()
    cpm(AFCMstats(), prior.count=2, log=TRUE)
  })
  
  
  # d3heatmap functions					
  plot <- reactive({
    df = AFCMstats()
    if (input$goButtonHeat == 0) {return(validate(
      need(input$filename != 0, "To generate a heatmap using AFCM Generation tool, please  select a *.CSV"),
      need(input$goButtonHeat !=0 , "for viewing statistically significant genes, please select desired parameters then  press 'Draw Heatmap'"),
      need(!is.null(df), "Unfortunately I can't generate a heatmap for you with these statistical cutoff parameters")
    ))}   
    else {
      d3heatmap( 
        if (input$log2_transform) log2_datasetInput() else AFCMstats(),
      #  cexRow = as.numeric(as.character(input$xfontsize)),
       # cexCol = as.numeric(as.character(input$yfontsize)),
      cexRow = 0.5,
      cexCol = 0.5,
      k_row = input$color_row_branches,
      k_col = input$color_column_branches,
      userdendrogram = input$userdendrogram,
      colors = input$choosecolor
       
      )  
    }
  })
  
  # edgeR data preparation for Normalization
  output$noncontrol <- renderUI({
    df <- datasetInput()
    if (is.null(df)) return(NULL)
    noncontrol <- names(df)
    selectInput('noncontrol', 'Please select the noncontrol samples:', choices = noncontrol, multiple = TRUE)
  })			
  
  
  #  edgeR function for statistics
  AFCMstatse <- reactive({
    if(!is.null(datasetInput()) & !is.null(input$noncontrol)){
      
      group <- as.numeric(names(datasetInput()) %in% input$noncontrol)
      y <- DGEList(counts = datasetInput(), group = group)
      y <- calcNormFactors(y)
      y <- estimateCommonDisp(y)
      y <- estimateTagwiseDisp(y)
      et <- exactTest(y)
      results <- topTags(et, n=1000000)
      results_df <- as.data.frame(results)
      
    }
  })
  # d3heatmap output
  output$heatmap <- renderD3heatmap({
    if(!is.null(datasetInput()))
      plot()
  })
 
  # d3heatmap download								
  output$downloadHeatmap <- downloadHandler(
    filename = function() {
      paste0(basename(file_path_sans_ext(input$filename)), '.html')
    },
    content = function(file) {
      saveWidget(plot(), file)
    }
  )
 
  
  
  
  
})
