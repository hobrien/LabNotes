#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Application title
#titlePanel("Gene Expression in the Fetal Brain: Sex Biases"),

navbarPage("Gene Expression in the Fetal Brain: Sex Biases:",
            tabPanel("Plot",
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
            ),
            tabPanel("Table",
                     sidebarLayout(
                       sidebarPanel(
                         radioButtons("p_type", "Maximum p-value", c('Uncorrected p-values' = 'p', 'FDR corrected p-values (q-values)'= 'q'), selected = 'q', inline = FALSE,
                                      width = NULL),
                         sliderInput("pvalue", "p-value:", 
                                     min = 0, max = 1, value = 0.1, step= 0.01)
                       ),
                       mainPanel(
                         tabsetPanel(
                           id = 'dataset',
                           tabPanel('12-19 PCW', DT::dataTableOutput('mytable1')),
                           tabPanel('12 PCW', DT::dataTableOutput('mytable2')),
                           tabPanel('13 PCW', DT::dataTableOutput('mytable3')),
                           tabPanel('14 PCW', DT::dataTableOutput('mytable4')),
                           tabPanel('15-16 PCW', DT::dataTableOutput('mytable5')),
                           tabPanel('17-19 PCW', DT::dataTableOutput('mytable6'))
                         )
                       )
                     )
            )
)
