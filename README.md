# AFCMHeatMapGeneratingTool

Generating heatmaps of genetic datasets is a 2D graphical visualization of data where the individual expression values contained in a matrix are represented as colors. Herein, we describe AFCMHeatMap a shiny web App that integrates quantitative interaction of genomics data and results from microarrays or RNA-Seq to highlight expression levels of various genetic datasets with a *.CSV input file. The application also facilitates downloading heatmaps as a supplementary material for user's publications. Written in R using Shiny framework, it is a user-friendly framework for interactive expression data  visualization that can be easily deployed without any restrictions to any operating system used by any online user.

Availability: It is an open-source and available at https://github.com/MTarekM/AFCMHeatMapGeneratingTool
Corresponding Author: Mohammad M. Tarek1
Contact: mohammadtareq459@gmail.com

The Authors are very thankful to the authors of the paper entitled “MicroScope: ChIP-seq and RNA-seq software analysis suite for gene expression heatmaps” [1], for providing their source code freely available on github[2], in which their developed functions have helped so much to use the required package for this ShinyApp.

In server file the AFCMstats Function was inspired by the authors SlimStats function submitted to their referenced repository with customization to limit visualization spaces and give more space to user specific zooming.
The referenced repositiory was also of great help considering Usage of edgeR package functions for data analysis.

References:

[1] Khomtchouk BB, Hennessy JR, Wahlestedt C. Microscope: ChIP-seq and RNA-seq software analysis suite for gene expression heatmaps. BMC Bioinformatics. 2016; 17(1):390.

[2] https://github.com/Bohdan-Khomtchouk/Microscope
