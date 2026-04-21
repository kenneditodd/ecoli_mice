library(dplyr)
library(plotly)
library(shiny)
library(ggplot2)
library(readr)
library(ggrepel)
source("functions/plot_boxplot.R")
source("functions/plot_correlation.R")
source("functions/plot_volcano.R")

server <- shinyServer(function(input, output, session) {
  
  # Gene selection
  observe({
    gene_options <- readRDS("data/gene_options.rds")
    updateSelectInput(session, "gene", choices = gene_options, selected = "Lcn2")
  })
  
  # DEG table selection
  observe({
    # DEG file options
    deg_files <- list.files("data/DEGs", pattern = "\\.tsv$", full.names = FALSE)
    # update for volcano tab
    updateSelectInput(session, "deg_table", choices = deg_files)
    # update for correlation tab
    updateSelectInput(session, "corr_deg_list1_or", choices = deg_files, 
                      selected = "E.2D.O.F_vs_S.2D.O.F_FDRq_1.00_LFC_0.00.tsv")
    updateSelectInput(session, "corr_deg_list2_or", choices = deg_files,
                      selected = "E.2D.Y.F_vs_S.2D.Y.F_FDRq_1.00_LFC_0.00.tsv")
    updateSelectInput(session, "corr_deg_list1", choices = deg_files, 
                      selected = "E.2D.O.F_vs_S.2D.O.F_FDRq_1.00_LFC_0.00.tsv")
    updateSelectInput(session, "corr_deg_list2", choices = deg_files,
                      selected = "E.2D.Y.F_vs_S.2D.Y.F_FDRq_1.00_LFC_0.00.tsv")
  })
  
  # Reactive CPM data
  get_cpm_data <- reactive({
    req(input$filtering)
    file_path <- paste0("data/CPM/CPM_", input$filtering, ".tsv")
    read.table(file_path, header = TRUE, row.names = 1, sep = "\t")
  })
  
  # CPM plot
  output$cpm_plot <- renderPlot({
    req(input$gene)
    plot_boxplot(get_cpm_data(), input$gene)
  })
  
  # CPM download
  output$download_cpm <- downloadHandler(
    filename = function() paste0("CPM_", input$filtering, ".tsv"),
    content = function(file) {
      write.table(get_cpm_data(), file, sep="\t", quote = FALSE, col.names = NA)
    }
  )
  
  # CPM PDF download
  output$download_pdf <- downloadHandler(
    filename = function() {
      paste0(input$gene, "_boxplot.pdf")
    },
    content = function(file) {
      pdf(file, width = input$pdf_width, height = input$pdf_height)
      print(plot_boxplot(get_cpm_data(), input$gene))
      dev.off()
    }
  )
  
  # Volcano plot
  output$volcano_plot_dynamic <- renderPlotly({
    plot_volcano(input$fdrq_cutoff, input$lfc_cutoff, input$deg_table,
                 labels = input$highlight_genes, dynamic = TRUE)
  })
  
  output$volcano_plot_static <- renderPlot({
    plot_volcano(input$fdrq_cutoff, input$lfc_cutoff, input$deg_table,
                 labels = input$highlight_genes, dynamic = FALSE)
  })
  
  output$volcano_plot <- renderUI({
    if (input$plot_type == "Dynamic") {
      plotlyOutput("volcano_plot_dynamic")
    } else {
      plotOutput("volcano_plot_static")
    }
  })
  
  # Volcano plot PDF download
  output$download_volcano_pdf <- downloadHandler(
    filename = function() {
      paste0(gsub("1.00_LFC_0.00.tsv", "", input$deg_table),
             sprintf("%.2f", input$fdrq_cutoff), "_LFC_", sprintf("%.2f", input$lfc_cutoff), ".pdf")
    },
    content = function(file) {
      req(input$deg_table)
      pdf(file, width = input$pdf_width, height = input$pdf_height)
      print(plot_volcano(input$fdrq_cutoff, input$lfc_cutoff, input$deg_table,
                         labels = input$highlight_genes, dynamic = FALSE))
      dev.off()
    }
  )
  
  # Correlation plo - detect OR/AND filtering
  output$correlation_plot <- renderUI({
    if (input$corr_filter_and_or == "OR") {
      plotlyOutput("correlation_plot_or")
    } else {
      plotlyOutput("correlation_plot_and")
    }
  })
  
  # correlation plot when filtering by OR
  output$correlation_plot_or <- renderPlotly({
    # required input
    req(input$corr_deg_list1_or, input$corr_deg_list2_or,
        input$corr_fdrq_cutoff, input$corr_lfc_cutoff)
    # plot correlation
    plot_correlation(deg_list1 = input$corr_deg_list1_or,
                     deg_list2 = input$corr_deg_list2_or,
                     fdrq1 = input$corr_fdrq_cutoff,
                     lfc1 = input$corr_lfc_cutoff)     
  })
  
  # correlation plot when filtering by AND
  output$correlation_plot_and <- renderPlotly({
    # required input
    req(input$corr_deg_list1, input$corr_fdrq_cutoff1, input$corr_lfc_cutoff1,
        input$corr_deg_list2, input$corr_fdrq_cutoff2, input$corr_lfc_cutoff2)
    # plot correlation
    plot_correlation(deg_list1 = input$corr_deg_list1,
                     fdrq1 = input$corr_fdrq_cutoff1,
                     lfc1 = input$corr_lfc_cutoff1,
                     deg_list2 = input$corr_deg_list2,
                     fdrq2 = input$corr_fdrq_cutoff2,
                     lfc2 = input$corr_lfc_cutoff2)  
  })
  
})
