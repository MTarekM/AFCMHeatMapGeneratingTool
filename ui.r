# Copyright (C) 2016-2017 Mohammad Tarek Mansour - AFCM BioInformatics Lab

# This file is the UI part of the platform.

# ------------------------------------------------------------------------------------

library(shiny)
library(d3heatmap)
library(RColorBrewer)
library(edgeR)
library(shinythemes)


ui <- shinyUI(fluidPage(theme = shinytheme("paper"),
    list(tags$head(HTML('<link rel="icon", href="afcm.png", 
                        type="image/png" />'))),
    div(style="padding: 1px 0px; width: '100%'",
        titlePanel(
          title="", windowTitle="AFCM HeatMap Tool"
        )
    ),
  pageWithSidebar(
    
  headerPanel(h1(a(img(src = "header.png", align = "center", width = "100%",height="250px"), href = "http://www.afcm.ac.eg/EN/"))),
  sidebarPanel( 
    helpText(a("AFCM HeatMap Generation tool Source Code on GitHub", href = "https://github.com/MTarekM/AFCMHeatMapGeneratingTool", target = "_blank")),
    downloadButton("downloadSample", label = "Download a Sample  File"),
    fileInput("filename", "please upload your input file:", accept = c('.csv')),
    uiOutput("noncontrol"),
    selectInput("bothPvalueFDR", "Modify Heatmap generation paramters", selected = "Pvalue", c("Pvalue", "FDR", "FDR and PValue" = "both")),  
   numericInput("userFDR",  "Please select Cutoff for FDR in Heatmap", 0.1, min = 0.001, max = 0.1, step = 0.001),
  numericInput("userPvalue", "Please select Cutoff for Pvalues in Heatmap", 0.05, min = 0.001, max = 0.5, step = 0.001),
    selectInput("choosecolor", "Select your Favourite Color Scheme:", c( "Oranges", "Greys", "Greens","YlOrRd", "Purples", "PuRd", "PuBu", "OrRd", "GnBu", "BuPu", "BuGn", "Blues")),
    checkboxInput("log2_transform", "log2 Transform Heatmap"),  
    selectInput("userdendrogram", "Clustering:", c("none","row", "column", "both")),
    actionButton("goButtonHeat", "Generate My Heatmap"),
    downloadButton("downloadHeatmap", "Download Heatmap")
  ),
    mainPanel( "My Heatmap", uiOutput("pixelation") )
    ),
  tags$footer("Developed at Bioinformatics & Computational Biology Lab - AFCM 2017 under the
GNU General Public License v3.0", align = "center", style = "
              position:bottom;
              bottom:0;
              width:100%;
              height:25px;   /* Height of the footer */
              color: black;
              padding: 10px;
              background-color: white;
              z-index: 1000;")
)

)
