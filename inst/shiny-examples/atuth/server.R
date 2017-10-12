#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  bf <- callModule(bfabric, "bfabric8",  applicationid = c(168, 204))
  
  observeEvent(input$load, {
    if (input$ssh){
      values$proteinGroups <- bfabricShiny:::.ssh_unzip(file = 'proteinGroups.txt',
                                                        zipfile=file.path('/srv/www/htdocs', input$relativepath))
    }else{
      zipfile <- file.path('/srv/www/htdocs', input$relativepath)
      if (file.exists(zipfile)){
        values$proteinGroups <- bfabricShiny:::.unzip(file = 'proteinGroups.txt', zipfile=zipfile)
      }else{message("zip file canot be found.")}
    }
    
    message(dim(values$proteinGroups))
    message("DONE")
  })
  
  output$load <- renderUI({
    actionButton("load", "load data", icon("upload"))
  })
  
  
  
  output$distPlot <- renderPlot({
    
    # generate bins based on input$bins from ui.R
    x    <- faithful[, 2] 
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white')
    
  })
  
})
