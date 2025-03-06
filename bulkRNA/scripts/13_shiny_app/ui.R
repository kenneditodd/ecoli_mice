library(shiny)
library(plotly)

shinyUI(fluidPage(
  titlePanel("E.coli Bulk RNA-seq Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        condition = "input.tabs == 'CPM'",
        radioButtons("filtering", "Select Filtering Method:", 
                     choices = c("prefiltering", "postfiltering")),
        selectizeInput("gene", "Select Gene:", choices = NULL, selected = "Lcn2"),
        numericInput("pdf_width", "Width:", value = 8, min = 4),
        numericInput("pdf_height", "Height:", value = 8, min = 4),
        downloadButton("download_pdf", "Download PDF"),
        downloadButton("download_cpm", "Download CPM Table")
      ),
      
      conditionalPanel(
        condition = "input.tabs == 'Volcano'",
        selectInput("deg_table", "Select DEG Table:", choices = NULL),
        numericInput("fdrq_cutoff", "FDRq Cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
        numericInput("lfc_cutoff", "log2(FC) Cutoff:", value = 0, min = -5, max = 5, step = 0.1),
        radioButtons("plot_type", "Select Plot Type:", choices = c("Static", "Dynamic")),
        conditionalPanel(
          condition = "input.plot_type == 'Static'",
          selectizeInput("highlight_genes", "Highlight Gene(s):", choices = NULL, multiple = TRUE),
          numericInput("volcano_pdf_width", "Width:", 8, min = 4),
          numericInput("volcano_pdf_height", "Height:", value = 6, min = 4),
          downloadButton("download_volcano_pdf", "Download Plot as PDF")
        )
      ),
      
      conditionalPanel(
        condition = "input.tabs == 'Correlation'",
        selectInput("deg_list1", "Group 1: Select DEG list", choices = NULL),
        numericInput("cor_fdrq_cutoff1", "Group 1: FDRq cutoff", value = 1, min = 0, max = 1, step = 0.01),
        numericInput("cor_lfc_cutoff1", "Group 1: log2(FC) cutoff", value = 0, min = -5, max = 5, step = 0.1),
        selectInput("deg_list2", "Group 2: Select DEG list", choices = NULL),
        numericInput("cor_fdrq_cutoff2", "Group 2: FDRq cutoff", value = 1, min = 0, max = 1, step = 0.01),
        numericInput("cor_lfc_cutoff2", "Group 2: log2(FC) cutoff:", value = 0, min = -5, max = 5, step = 0.1),
        p("You can download the plot as a png in the interactive viewer")
      )
    ),
    
    mainPanel(
      tabsetPanel(id = "tabs",
                  tabPanel("CPM", plotOutput("cpm_plot")),
                  tabPanel("Volcano", uiOutput("volcano_plot")),
                  tabPanel("Correlation", plotlyOutput("correlation_plot"))
      )
    )
  )
))
