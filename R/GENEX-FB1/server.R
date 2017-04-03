#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
source("FormatGGplot.R")

counts <- read_delim("./Data/counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
fitted <- read_delim("./Data/fitted.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
target <- read_delim("./Data/target.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(label=as.character(label))

PlotTimepoint<-function(geneID, counts, fitted, target, ages) {
  ageSplit <- strsplit(ages, '-')[[1]]
  min <- ageSplit[1]
  max <- ageSplit[length(ageSplit)]
  data <- counts %>% filter(SYMBOL == geneID ) %>%  
    dplyr::select(-SYMBOL) %>%
    gather() %>%
    separate(key, into=c('norm', 'label'), sep='[.]') %>%
    dplyr::select(label, value) %>%
    left_join(target) %>%
    filter(PCW >= min & PCW <= max)
  mean <- filter(fitted, SYMBOL == geneID & ageBin==ages)
  fc <- mean$log2FoldChange[1] %>% format(digits=2)
  pval <- mean$padj[1] %>% format(digits=2)
  mean <- mean %>% dplyr::select(Male, Female) %>%
    gather('Sex', 'mean')
  title<-paste0(geneID, ': log2 difference = ', fc, ', adjusted pvalue = ', pval)
  plot<-  ggplot(data, aes(x=Sex, colour=Sex)) + 
    geom_jitter(aes(y=value), height = 0, width=.1, alpha=.75) + 
    geom_errorbar(aes(ymin=mean, ymax=mean), colour='black', size=1, width=.5, data=mean) +
    ylab("normalised counts") +
    xlab('') +
    main_theme() +
    scale_colour_brewer(type = "qual", palette = 6) +
    ggtitle(title) 
  plot
}

PlotSampleSize<-function(target, ages){
  ageSplit <- strsplit(ages, '-')[[1]]
  min <- ageSplit[1]
  max <- ageSplit[length(ageSplit)]
  target2 <- target %>% filter(PCW >= min & PCW <= max) %>% mutate(age_bin =ifelse(PCW == 16, 15, 
                                           ifelse(PCW > 16, 16, PCW)))
  
  plot <- ggplot(target2, aes(x=age_bin, fill=Sex)) +
    geom_bar() +
    facet_grid(Sex ~ .) +
    side_theme() +
    ggtitle("Sample Size") +
    scale_y_continuous(breaks=seq(0,100, 5)) +
    scale_x_continuous(breaks=c(11,12,13,14,15,16),
                       labels=c('11','12','13','14','15-16','17-19')) +
    tufte_theme() +
    xlab("Post Conception Weeks") +
    scale_fill_brewer(type = "qual", palette = 6) 
  plot
}

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  output$distPlot <- renderPlot({
    PlotTimepoint(input$geneID, counts, fitted, target, input$ages)
  })
  output$sampleSizeHist <- renderPlot({
    PlotSampleSize(target, input$ages)
  })
  
})
