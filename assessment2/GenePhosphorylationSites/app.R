#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(vegan)
library(dplyr)
library(stringr)
library(seqinr)
library(R.utils)
library(SparkR)
library(R.oo)
library(ggplot2)
library(ggrepel)

gene_data <- read.csv(file = "../plottingdata.csv", stringsAsFactors = FALSE, header = TRUE)


# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("Gene Phosphorylation Sites"),
    
    
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput(inputId = "gene", label = "Gene", choices = gene_data$external_gene_name)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("genePlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    gene_data <- read.csv(file = "../plottingdata.csv", stringsAsFactors = FALSE, header = TRUE)
    
    output$genePlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        gene <- input$gene
        SiteName <- dplyr::filter(gene_data, external_gene_name == gene) %>%
            dplyr::select(Sites) %>%
            str_split(pattern = ",")
        plotdata <- data.frame(Position = as.numeric(str_sub(as.vector(SiteName[[1]]),2)),
                               ones = rep(as.numeric(1), length(as.vector(SiteName[[1]]))),
                               Site = as.vector(SiteName[[1]]))
        
        ggplot(plotdata, aes(Position, ones)) +
            
#            geom_rect(aes(xmin = 1, ymin = 0, 
#                          xmax = length(as.vector(SiteName[[1]])), ymax = 20)) +
            geom_point() +
            geom_label_repel(aes(Position, label = plotdata$Site)) +
            coord_cartesian(xlim = c(0, gene_data[gene_data$external_gene_name == gene, "Seq_length"]))
        
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
