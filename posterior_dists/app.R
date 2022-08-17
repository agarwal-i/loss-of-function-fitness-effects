
library(data.table)
library(stringr)
library(ggplot2)
library(ggrepel)
library(splitstackshape)
library(gridExtra)
library(gridGraphics)
library(grid)
library(gtable)
library(ggpubr)
library(cowplot)
library(scales)
library(RColorBrewer)
library(ggridges)
library(ggsci)
library("shiny")


ui <- fluidPage(
  titlePanel("Sex-averaged strength of selection on the loss of one copy"),
  
  sidebarLayout(
    sidebarPanel(
      
      textInput("caption", "Enter Gene name", "BRCA1"),
      helpText("Some examples are: BRCA1, DMC1,", 
               "PCSK9, APOE, ACE2, DDX3X, etc."),
      verbatimTextOutput("my_text1"),
      verbatimTextOutput("my_text2"),
    ),
    
    mainPanel(
      plotOutput("my_plot")
    )
  )
)

server <- function(input, output) {
  
  posterior_summary = fread("Table_S2.txt", header=T)
  posterior_summary$map = round(10^posterior_summary$log10_map,7)
  posterior_summary$ci = paste0("(", round(10^posterior_summary$log10_ci_low,7),",", round(10^posterior_summary$log10_ci_high,7),")") 
  
  output$my_text1 <- renderText({paste0("MAP estimate = ",posterior_summary$map[posterior_summary$Gene == paste0(req(input$caption))]) })
  output$my_text2 <- renderText({paste0("Credible interval = ",posterior_summary$ci[posterior_summary$Gene == paste0(req(input$caption))]) })
  
  yax = c(readRDS("yax.Rda"))
  xax = c(readRDS("xax.Rda"))
  my_data <- reactive({
    yax1=yax[names(yax) %in% paste0(req(input$caption))]
    test = data.frame(xax, yax1)
    colnames(test) = c("hs","density")
    test
  })
  
  output$my_plot <- renderPlot(
    ggplot(my_data(), aes(x=10^hs, y=density)) + geom_point(size=1.2) + ggtitle(req(input$caption)) +
      theme_bw() + xlab("hs") +
      scale_color_manual(values=c("black")) +
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
  )
}
shinyApp(ui = ui, server = server)
