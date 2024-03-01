#
# Shiny app build by Lisa Buche - April 2023 

#This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#################################################-
## USER INTERFACE
#################################################-
## Preliminaries ----
#################################################-
library(shiny)
library(shinyMatrix)
#library(rgl)
library(scales)
library(tidyverse)
library(medicaldata)
library(shiny)
library(mlbench)
library(shiny)
library(plotly)
library(ggplot2)
library(ggpubr)
library(rmarkdown)
library(knitr)
library(pander)
library(ggforce)
#################################################-
## Define UI ----
#################################################-
ui <- fluidPage(
  # Application title
  titlePanel(
    h1("Population dynamics in a community of 2 plant species with changing abiotic conditions", h2("With a changing functional form of species interactions"))
  ),
  # HTML tags
  #tags$head(
  #  tags$style(HTML("hr {border-top: 1px solid #000000;}"))
  #),
  # User inputs
  sidebarLayout(
    sidebarPanel(
      h4("Species intrinsic growth rates:"),
      sliderInput("Lambda_i", label ="Lambda_i", min = 1, max = 10,
                  value=1,step = 1),
      sliderInput("Lambda_j", label ="Lambda_j", min = 1, max = 10,
                  value=1,
                  step = 1),
      sliderInput("Ni0",
                  label =withMathJax("Initial abundance of species i: \\(N_{i,t=1}:\\ "),
                  min = 0,
                  max = 10,
                  value = 1,
                  step=1),
      sliderInput("Nj0",
                  label =withMathJax("Initial abundance of species j: \\(N_{j,t=1}:\\ "),
                  min = 0,
                  max = 10,
                  value = 1,
                  step=1),
      sliderInput("RangeN", 
                  label = withMathJax("Range of Neighbour density \\(N_{ }:\\)"),
                  min = -10, max = 50, value = c(0, 20),
                  step = 1),
      sliderInput("time", 
                  label = withMathJax("Number of generations"),
                  min = 1, max = 100, value = c(20),
                  step = 1),
      sliderInput("b.alpha", 
                  label = withMathJax("Impact of abiotic conditions on initial interaction coefficients"),
                  min = 0, max = 1, value = c(0, 0),
                  step = 0.01),
      sliderInput("b.lambda", 
                  label = withMathJax("Impact of abiotic conditions on intrinsic growth rate"),
                  min = 0, max = 10, value = c(0, 0),
                  step = 0.5),
      h4("Initial interaction coefficients:"),
      shinyMatrix::matrixInput(inputId = "alphamat_init",
                               value = matrix(c(-0.2, -0.4,
                                                -0.5, -0.6),
                                              nrow = 2,
                                              dimnames = list(c("α_0(i,_)", "α_0(j,_)"), c("α_0(_,i)", "α_0(_,j)")),
                                              byrow = TRUE),
                               class = "numeric",
                               rows = list(names = TRUE,
                                           editableNames = FALSE),
                               cols = list(names = TRUE,
                                           editableNames = FALSE)),
      hr(),
      h4("Density-dependent interaction coefficients:"),
      shinyMatrix::matrixInput(inputId = "alphamat_slope",
                               value = matrix(c(0, 0,
                                                0, 0),
                                              nrow = 2,
                                              dimnames = list(c("α_(i,_)", "α_(j,_)"), c("α_(_,i)", "α_(_,j)")),
                                              byrow = TRUE),
                               class = "numeric",
                               rows = list(names = TRUE,
                                           editableNames = FALSE),
                               cols = list(names = TRUE,
                                           editableNames = FALSE)),
      hr(),
      h4("Stretching parameter:"),
      shinyMatrix::matrixInput(inputId = "Cmat",
                               value = matrix(c(0, 0,
                                                0, 0),
                                              nrow = 2,
                                              dimnames = list(c("C(i,_)", "C(j,_)"), c("C(_,i)", "C(_,j)")),
                                              byrow = TRUE),
                               class = "numeric",
                               rows = list(names = TRUE,
                                           editableNames = FALSE),
                               cols = list(names = TRUE,
                                           editableNames = FALSE)),
      hr(),
      h4("Optimal density:"),
      shinyMatrix::matrixInput(inputId = "Nomat",
                               value = matrix(c(1, 5,
                                                5, 1),
                                              nrow = 2,
                                              dimnames = list(c("N_opt(i,_)", "N_opt(j,_)"), c("N_opt(_,i)", "N_opt(_,j)")),
                                              byrow = TRUE),
                               class = "numeric",
                               rows = list(names = TRUE,
                                           editableNames = FALSE),
                               cols = list(names = TRUE,
                                           editableNames = FALSE)),
        ),
    # Outputs
  mainPanel(
    fluidRow(
                   plotOutput("Afunction")
                 ),
                
                 fluidRow(
                   plotOutput("fecundity")
                 ),
            
                 fluidRow(
                   plotOutput("timeserie")
                 ),
           
      )
    )
)


server <- function(input, output){
  
  output$Afunction <- renderPlot({
    Ainit =input$alphamat_init
    Aslope = input$alphamat_slope
    C = input$Cmat
    No = input$Nomat
    N = seq(input$RangeN[1],
             input$RangeN[2],0.25)
    b.alpha = c(input$b.alpha[1],
            input$b.alpha[2])
    b.lambda = c(input$b.lambda[1],
                input$b.lambda[2])
    
    functions.plot(Ainit, Aslope,
                   C,N,No,b.alpha)
  }, res = 96)
  
  output$fecundity <- renderPlot({
    lambda =c(input$Lambda_i,input$Lambda_j)
    Ainit =input$alphamat_init
    Aslope = input$alphamat_slope
    C = input$Cmat
    No = input$Nomat
    N = seq(input$RangeN[1],
            input$RangeN[2],1)
    b.alpha = c(input$b.alpha[1],
                input$b.alpha[2])
    b.lambda = c(input$b.lambda[1],
                 input$b.lambda[2])
    
    fecundity.plot(lambda,Ainit, Aslope,
                   C ,N,No,b.alpha, b.lambda)
  }, res = 96)
  
  output$timeserie <- renderPlot({
    lambda =c(input$Lambda_i,input$Lambda_j)
    N0 = c(input$Ni0,input$Nj0)
    Ainit =input$alphamat_init
    Aslope = input$alphamat_slope
    C = input$Cmat
    No = input$Nomat
    t = input$time
    b.alpha = c(input$b.alpha[1],
                input$b.alpha[2])
    b.lambda = c(input$b.lambda[1],
                 input$b.lambda[2])
    
    
    Abundance.plot(lambda,Ainit, Aslope,
                   C ,N0,No, t,b.alpha,b.lambda)
  }, res = 96)
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

