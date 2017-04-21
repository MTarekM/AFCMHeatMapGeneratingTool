# Copyright (C) 2016-2017 Mohammad Tarek Mansour - AFCM BioInformatics Lab

# This file is the server file of the Application

# ------------------------------------------------------------------------------------

library(shiny)
library(d3heatmap)
library(htmlwidgets)
library(tools)
library(edgeR)
library(dplyr)
library(png)

# backend 
server <- shinyServer(function(input, output) {	
  
  
  
  # sample file download
  output$downloadData <- downloadHandler(
    filename <- function() {
      paste('SampleFile', '.csv', sep='')
    },
    content <- function(file) {
      file.copy("SampleFile.csv", file)
    },
    contentType = "text/csv"
  )
  
  
  # file upload
  datasetInput <- reactive({
    validate(
      need(input$filename != 0, "To generate a heatmap using AFCM Generation tool, please upload your *.CSV file") 
    )
    inFile <- input$filename
    if (is.null(inFile)) return(NULL)
    read.table(inFile$datapath, header= TRUE, sep=",", quote='"', row.names=1)
  })
  
  
  # filter stats table based on cutoffs
  slimStats <- reactive({
    results_df <- stats()
    datasetInputTwo <- datasetInput()
    
    if (input$pvFDRchoose == "Pvalue"){
      pNewResults_df <- results_df[which(results_df$PValue<input$statPV),]
      pNames <- row.names(pNewResults_df)
      pFinalResult <- datasetInputTwo[pNames,]
    }
    else if(input$pvFDRchoose == "FDR"){
      fNewResults_df <- results_df[which(results_df$FDR<input$statFDR),]
      fNames <- row.names(fNewResults_df)
      fFinalResult <- datasetInputTwo[fNames,]
    }
    else if(input$pvFDRchoose == "both"){
      bNewResults_df <- results_df[which(results_df$FDR<input$statFDR & results_df$PValue<input$statPV),]
      bNames <- row.names(bNewResults_df)
      bFinalResult <- datasetInputTwo[bNames,]
    }
  })
  
  
  # heatmap height
  output$pixelation <- renderUI({
     
    slimStats()
    inputLines <- NROW(slimStats())
    if(inputLines >= 0 && inputLines <= 2001){
      d3heatmapOutput("heatmap", width = "100%", height = "900px")
    }
    else if(inputLines > 2001 && inputLines <= 5001){
      d3heatmapOutput("heatmap", width = "100%", height = "1700px")
    }
    else if(inputLines > 5001 && inputLines <= 10001){
      d3heatmapOutput("heatmap", width = "100%", height = "3200px")
    }
    else if(inputLines > 10001 && inputLines <= 20001){
      d3heatmapOutput("heatmap", width = "100%", height = "9000px")
    }
    else if(inputLines > 20001 && inputLines <= 40001){
      d3heatmapOutput("heatmap", width = "100%", height = "13000px")
    }
    else if(inputLines > 40001 && inputLines <= 60001){
      d3heatmapOutput("heatmap", width = "100%", height = "18000px")
    }
    else
      d3heatmapOutput("heatmap", width = "100%", height = "30000px")
    
  })	
  
  
  # log2 data transformation
  log2_datasetInput <- reactive({
    slimStats()
    cpm(slimStats(), prior.count=2, log=TRUE)
  })
  
  
  # d3heatmap prep					
  plot <- reactive({
    df = slimStats()
    if (input$goButtonHeat == 0) {return(validate(
      need(input$filename != 0, "To generate a heatmap using AFCM Generation tool, please  select a *.CSV"),
      need(input$goButtonHeat !=0 , "for viewing statistically significant genes, please select desired parameters then  press 'Draw Heatmap'"),
      need(!is.null(df), "Unfortunately I can't generate a heatmap for you with these statistical cutoff parameters")
    ))}   
    else {
      d3heatmap( 
        if (input$log2_transformed_data) log2_datasetInput() else slimStats(),
        cexRow = as.numeric(as.character(input$xfontsize)),
        cexCol = as.numeric(as.character(input$yfontsize)),
        colors = input$choose,
        k_row = input$color_row_branches,
        k_col = input$color_column_branches,
        dendrogram = input$dendrogram
      )  
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
  
  
 
  
  
  # edgeR prep
  output$expcolumns <- renderUI({
    df <- datasetInput()
    if (is.null(df)) return(NULL)
    expcolumns <- names(df)
    selectInput('expcolumns', 'Please select the non-control samples:', choices = expcolumns, multiple = TRUE)
  })			
  
  
  # edgeR statistical engine
  stats <- reactive({
    if(!is.null(datasetInput()) & !is.null(input$expcolumns)){
      
      group <- as.numeric(names(datasetInput()) %in% input$expcolumns)
      y <- DGEList(counts = datasetInput(), group = group)
      y <- calcNormFactors(y)
      y <- estimateCommonDisp(y)
      y <- estimateTagwiseDisp(y)
      et <- exactTest(y)
      results <- topTags(et, n=50000)
      results_df <- as.data.frame(results)
      
    }
  })
  
  
 
  
  
  
  
  
})