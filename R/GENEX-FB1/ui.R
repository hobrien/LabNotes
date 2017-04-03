#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Gene Expression in the Fetal Brain: Sex Biases"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      textInput("geneID", "geneID"),
      radioButtons("ages", "Post-Conception Weeks", c('12', '13', '14', '15-16', '17-19', '12-19'), selected = '12-19', inline = FALSE,
                   width = NULL),
      plotOutput("sampleSizeHist", height=200)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
))
