library(shiny)
library(SpiecEasi)
library(SPRING)
library(CMiNet)
library(igraph)
library(visNetwork)

server <- function(input, output, session) {
  data <- reactive({
    req(input$file)  # Ensure a file is uploaded
    df <- read.csv(input$file$datapath, row.names = 1, check.names = FALSE)
    isCount <- input$countBased
    if (!isCount) {
      df <- df + 0.0001
    }
    df  # Return the modified data frame
  })

  
  
  taxa_name <- reactive({
    req(input$file)
    df <- read.csv(input$file$datapath, row.names = 1, check.names = FALSE)
    taxa_name_matrix <- matrix(0, nrow = ncol(df), ncol = 2)
    taxa_name_matrix[, 1] <- colnames(df)
    taxa_name_matrix[, 2] <- 1:ncol(df)
    colnames(taxa_name_matrix) <- c("original", "taxa_id")
    taxa_name_matrix
  })
  
  algorithmParams <- list(
    sparcc = list(
      numericInput("imax", "imax", value = 20,min=1,max=100,step=1),
      numericInput("kmax", "kmax", value = 10,min=1,max=100,step=1),
      numericInput("alpha", "alpha", value = 0.1,min=0,max=1,step=0.1),
      numericInput("Vmin", "Vmin", value = 1e-4,min=1e-4,max=1,step=0.01)
    ),
    spiecEasi_mb = list(
      numericInput("lambda_min_ratio_mb", "Lambda Min Ratio", value = 1e-2,min=0,max=1,step=0.01),
      numericInput("nlambda_mb", "Number of Lambdas", value = 15,min=1,max=100,step=1),
      numericInput("rep_num_mb", "Pulsar Rep Num", value = 20,min=1,max=100,step=1),
      numericInput("ncores_mb", "Pulsar Cores", value = 4,min=1,max=20,step=1)
    ),
    spiecEasi_glasso = list(
      numericInput("lambda_min_ratio_glasso", "Lambda Min Ratio", value = 1e-2,min=0,max=1,step=0.01),
      numericInput("nlambda_glasso", "Number of Lambdas", value = 15,min=1,max=100,step=1),
      numericInput("rep_num_glasso", "Pulsar Rep Num", value = 50,min=1,max=100,step=1)
    ),
    spring = list(
      numericInput("ncores_spring", "Number of Cores", value = 5,min=1,max=20,step=1),
      numericInput("nlambda_spring", "Number of Lambdas", value = 15,min=1,max=100,step=1),
      numericInput("rep_num_spring", "Rep Num", value = 20,min=1,max=100,step=1)
    ),
    gcoda = list(
      numericInput("pseudo_gcoda", "Pseudo", value = 0.5,min=0,max=2, step = 0.1),
      numericInput("lambda_min_ratio_gcoda", "Lambda Min Ratio", value = 1e-4),
      numericInput("nlambda_gcoda", "Number of Lambdas", value = 15,min=1,max=100,step=1),
      numericInput("ebic_gamma_gcoda", "EBIC Gamma", value = 0.5,min=0,max=1,step=0.01)
    ),
    c_MI = list(
      numericInput("q1", "Q1", value = 0.7,min=0,max=1, step = 0.1),
      numericInput("q2", "Q2", value = 0.95,min=0,max=1, step = 0.1)
    ),
    cclasso = list(
      numericInput("pseudo_cclasso", "Pseudo", value = 0.5,min=0,max=2, step = 0.1),
      numericInput("k_cv_cclasso", "K CV", value = 3,min=1,max=10, step = 1),
      numericInput("k_max_cclasso", "K Max", value = 70,min=1,max=200, step = 5),
      numericInput("n_boot_cclasso", "Number of Bootstraps", value = 20, min=1,max=200, step = 5)
    ),
    spearman = list(),
    pearson = list(),
    bicor = list()
  )
  
  output$algorithmParams <- renderUI({
    req(input$selectedAlgorithms)
    tagList(
      lapply(input$selectedAlgorithms, function(algorithm) {
        params <- algorithmParams[[algorithm]]
        if (length(params) > 0) {
          tags$fieldset(
            tags$legend(algorithm),
            do.call(tagList, params)
          )
        } else {
          tags$fieldset(
            tags$legend(algorithm),
            p("No parameters to configure for this algorithm.")
          )
        }
      })
    )
  })
  
  output$dataSummary <- renderTable({
    req(data())
    head(data()[, 1:15])
  })
  
  output$dataInfo <- renderUI({
    req(data())
    df <- data()
    num_samples <- nrow(df)
    num_taxa <- ncol(df)
    
    tagList(
      p("Rows: Samples, Columns: Taxa"),
      p(paste("Number of taxa =", num_taxa)),
      p(paste("Number of samples =", num_samples))
    )
  })
  
  
  output$taxaSummary <- renderTable({
    req(taxa_name())
    taxa_name()
  })
  
  results <- eventReactive(input$run, {
    req(data())
    selectedAlgorithms <- input$selectedAlgorithms
    isCount <- input$countBased
    TT <- as.numeric(input$TT)
    dir.create("Network", showWarnings = FALSE)
    dir.create("Binary_Network", showWarnings = FALSE)
    
    algorithmInputs <- lapply(names(algorithmParams), function(algorithm) {
      if (algorithm %in% selectedAlgorithms) {
        params <- lapply(names(algorithmParams[[algorithm]]), function(param) {
          input[[param]]
        })
        names(params) <- names(algorithmParams[[algorithm]])
        list(enabled = TRUE, params = params)
      } else {
        list(enabled = FALSE, params = list())
      }
    })
    names(algorithmInputs) <- names(algorithmParams)
    
    c_MI_params <- if ("c_MI" %in% selectedAlgorithms) {
      list(
        enabled = TRUE,
        params = list(
          quantitative = input$quantitative,
          q1 = input$q1,
          q2 = input$q2
        )
      )
    } else {
      list(enabled = FALSE, params = list())
    }
    
    gcoda_params <- if ("gcoda" %in% selectedAlgorithms) {
      list(
        enabled = TRUE,
        params = list(
          counts = isCount,
          pseudo = input$pseudo_gcoda,
          lambda.min.ratio= input$lambda_min_ratio_gcoda,
          nlambda= input$nlambda_gcoda,
          ebic.gamma = input$ebic_gamma_gcoda
        )
      )
    } else {
      list(enabled = FALSE, params = list())
    } 
    cclasso_params <- if ("cclasso" %in% selectedAlgorithms) {
      list(
        enabled = TRUE,
        params = list(
          counts = isCount,
          pseudo = input$pseudo_cclasso,
          k_cv = input$k_cv_cclasso,
          lam_int= c(1e-4, 1),
          k_max= input$k_max_cclasso,
          n_boot = input$n_boot_cclasso
        )
      )
    } else {
      list(enabled = FALSE, params = list())
    } 

    spiecEasi_mb_params <- if ("spiecEasi_mb" %in% selectedAlgorithms) {
      list(
        enabled = TRUE,
        params = list(
          method = 'mb',
          lambda.min.ratio = input$lambda_min_ratio_mb,
          nlambda = input$nlambda_mb,
          pulsar.params = list(rep.num = input$rep_num_mb, ncores = input$ncores_mb) 
        )
      )
    } else {
      list(enabled = FALSE, params = list())
    } 
    
    spiecEasi_glasso_params <- if ("spiecEasi_glasso" %in% selectedAlgorithms) {
      list(
        enabled = TRUE,
        params = list(
          method = 'glasso',
          lambda.min.ratio  = input$lambda_min_ratio_glasso,
          pulsar.params = list(rep.num =input$rep_num_glasso)
        )
      )
    } else {
      list(enabled = FALSE, params = list())
    }
    
    sparcc_params <- if ("sparcc" %in% selectedAlgorithms) {
      if (isCount) {
        list(
          enabled = TRUE,
          params = list(
            imax = input$imax,
            kmax = input$kmax,
            alpha = input$alpha,
            Vmin = input$Vmin
          )
        )
      } else {
        list(
          enabled = TRUE,
          params = list(
            kmax = input$kmax,
            alpha = input$alpha,
            Vmin = input$Vmin
          )
        )
      }
    } else {
      list(enabled = FALSE, params = list())
    }
    
    spring_params <- if ("spring" %in% selectedAlgorithms) {
      list(
        enabled = TRUE,
        params = list(
          Rmethod = "original",
          quantitative = isCount,
          ncores = input$ncores_spring,
          lambdaseq = "data-specific",
          nlambda = input$nlambda_spring,
          rep.num = input$rep_num_spring
        )
      )
    } else {
      list(enabled = FALSE, params = list())
    }
    
    result <- CMiNet(
      data = data(),
      quantitative = isCount,
      TT = TT,
      pearson = algorithmInputs[["pearson"]],
      spearman = algorithmInputs[["spearman"]],
      bicor = algorithmInputs[["bicor"]],
      sparcc = sparcc_params,
      spiecEasi_mb = spiecEasi_mb_params ,
      spiecEasi_glasso = spiecEasi_glasso_params,
      spring = spring_params,
      gcoda = gcoda_params,
      c_MI = c_MI_params,
      cclasso = cclasso_params
    )
    
    lapply(names(result), function(algorithm) {
      if (is.list(result[[algorithm]]) && !is.null(result[[algorithm]]$binary_matrix)) {
        write.csv(result[[algorithm]]$binary_matrix, file.path("Binary_Network", paste0(algorithm, "_binary.csv")))
        write.csv(result[[algorithm]]$adj_matrix, file.path("Network", paste0(algorithm, ".csv")))
      }
    })
    
    # Save weighted network to CSV
    if (!is.null(result$weighted_network)) {
      write.csv(result$weighted_network, file.path("Network", "weighted_network.csv"))
    }
    result
  })
  
  output$resultsTab <- renderUI({
    req(results())
    tagList(
      selectInput("selectedResultAlgorithm", "Select Algorithm to Download Result:", choices = c(names(results()), "weighted_network")),
      downloadButton("downloadSelectedBinary", "Download Selected Binary Matrix"),
      downloadButton("downloadSelectedNetwork", "Download Selected Network Matrix"),
      downloadButton("downloadWeightedNetwork", "Download Weighted Network"),
      downloadButton("downloadEdgeList", "Download Edge List")
    )
  })
  
  
  # Download handler for the sample data
  output$downloadSampleData <- downloadHandler(
    filename = function() {
      "sample_data.csv"  # The name of the file the user will download
    },
    content = function(file) {
      file.copy("sample_data.csv", file)
    }
  )
  
  
  # Download handlers for selected algorithm results in Results tab
  output$downloadSelectedBinary <- downloadHandler(
    filename = function() {
      paste0(input$selectedResultAlgorithm, "_binary.csv")
    },
    content = function(file) {
      req(results()[[input$selectedResultAlgorithm]])
      write.csv(results()[[input$selectedResultAlgorithm]]$binary_matrix, file, row.names = TRUE)
    }
  )
  output$downloadSelectedNetwork <- downloadHandler(
    filename = function() {
      paste0(input$selectedResultAlgorithm, ".csv")
    },
    content = function(file) {
      req(results()[[input$selectedResultAlgorithm]])
      write.csv(results()[[input$selectedResultAlgorithm]]$adj_matrix, file, row.names = TRUE)
    }
  )
  output$downloadWeightedNetwork <- downloadHandler(
    filename = function() {
      "weighted_network.csv"
    },
    content = function(file) {
      req(results())
      write.csv(results()$weighted_network, file, row.names = TRUE)
    }
  )
  output$downloadEdgeList <- downloadHandler(
    filename = function() {
      "edge_list.csv"
    },
    content = function(file) {
      req(results())
      edge_list <- which(results()$weighted_network > 0, arr.ind = TRUE)
      edges <- data.frame(
        from = rownames(results()$weighted_network)[edge_list[, 1]],
        to = colnames(results()$weighted_network)[edge_list[, 2]]
      )
      write.csv(edges, file, row.names = FALSE)
    }
  )
  
  observeEvent(input$runVisualization, {
    req(input$weightedFile)
    weighted_network <- tryCatch({
      as.matrix(read.csv(input$weightedFile$datapath, row.names = 1, check.names = FALSE))
    }, error = function(e) {
      showNotification("Error reading weighted network file.", type = "error")
      return(NULL)
    })
    req(weighted_network)
    
    if (nrow(weighted_network) != ncol(weighted_network)) {
      showNotification("Error: Weighted network matrix must be square", type = "error")
      return(NULL)
    }
    
    thresholds <- c(input$threshold1, input$threshold2, input$threshold3, input$threshold4)
    show_labels <- c(1 %in% input$showLabels, 2 %in% input$showLabels, 3 %in% input$showLabels, 4 %in% input$showLabels)
    node_colors <- unlist(strsplit(input$nodeColors, ","))
    edge_colors <- unlist(strsplit(input$edgeColors, ","))
    
    output$networkPlot <- renderPlot({
      process_and_visualize_network(
        weighted_network,
        taxa_name()[, 2],
        thresholds,
        show_labels,
        node_colors,
        edge_colors
      )
    })
  })
  
  
  
  observeEvent(input$runFinalNetwork, {
    req(input$finalWeightedFile)
    weighted_network <- tryCatch({
      as.matrix(read.csv(input$finalWeightedFile$datapath, row.names = 1, check.names = FALSE))
    }, error = function(e) {
      showNotification("Error reading weighted network file.", type = "error")
      return(NULL)
    })
    req(weighted_network)
    
    if (nrow(weighted_network) != ncol(weighted_network)) {
      showNotification("Error: Weighted network matrix must be square", type = "error")
      return(NULL)
    }
    
    score_threshold <- input$score
    final_network <- ifelse(weighted_network > score_threshold, 1, 0)
    final_network[lower.tri(final_network)] <- 0
    
    output$finalNetworkPlot <- renderVisNetwork({
      req(final_network)
      edge_list <- which(final_network == 1, arr.ind = TRUE)
      edges <- data.frame(
        from = rownames(final_network)[edge_list[, 1]],
        to = colnames(final_network)[edge_list[, 2]]
      )
      nodes <- data.frame(
        id = unique(c(edges$from, edges$to)),
        label = unique(c(edges$from, edges$to)),
        color = input$finalNodeColor
      )
      visNetwork(nodes, edges, main = "Final Microbiome Network") %>%
        visNodes(color = list(border = "black", highlight = "orange"), font = list(color = input$finalLabelColor)) %>%
        visEdges(color = list(color = input$finalEdgeColor, highlight = "orange")) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visLayout(randomSeed = 123)
    })
    
    output$downloadFinalNetworkPlot <- downloadHandler(
      filename = function() {
        paste("final_network_plot", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        png(file, width = 1200, height = 800, res = 300) 
        # Retrieve and modify `final_network`
        plot_network(
          final_network,
          node_color = input$finalNodeColor,
          edge_color = input$finalEdgeColor,
          label_color = input$finalLabelColor
        )
        dev.off()
      }
    )
  })
  
  output$downloadFinalWeightedNetwork <- downloadHandler(
    filename = function() {
      "final_weighted_network.csv"
    },
    content = function(file) {
      req(input$finalWeightedFile)
      file.copy(input$finalWeightedFile$datapath, file)
    }
  )
}

process_and_visualize_network <- function(weighted_network, taxa_names, thresholds, show_labels, node_colors, edge_colors) {
  par(mfrow = c(2, 2))  # Arrange plots in a 2x2 grid
  for (i in 1:length(thresholds)) {
    threshold <- thresholds[i]
    network <- weighted_network
    network[network <= threshold] <- 0
    graph <- graph_from_adjacency_matrix(as.matrix(network), mode = "undirected", weighted = TRUE, diag = FALSE)
    
    # Set node and edge attributes
    V(graph)$color <- node_colors[i]
    E(graph)$color <- edge_colors[i]
    
    set.seed(123) 
    plot(
      graph,
      vertex.label = ifelse(show_labels[i], taxa_names, NA),
      vertex.size = 5,
      main = paste("Threshold =", threshold))
    mtext(side = 1, line = 4, paste("NN degree >0: ", sum(degree(graph) > 0), ", Edges: ", ecount(graph), ", Max-Degree: ", max(degree(graph))),
          font = 3) # Change font style here)
  }
}

plot_network <- function(network_final, node_color = "skyblue", edge_color = "grey", label_color = "black") {
  edge_list <- which(network_final == 1, arr.ind = TRUE)
  edges <- data.frame(
    from = rownames(network_final)[edge_list[, 1]],
    to = colnames(network_final)[edge_list[, 2]]
  )
  network_graph <- graph_from_data_frame(d = edges, directed = FALSE, vertices = data.frame(name = unique(c(edges$from, edges$to))))
  set.seed(123)  # Set seed for reproducibility of layout
  layout <- layout_with_fr(network_graph, niter = 500, grid = "nogrid") * 3  # Scale layout
  
  plot(network_graph, layout = layout, vertex.color = node_color, vertex.size = 5,
       vertex.label.color = label_color, vertex.label.cex = 0.6, vertex.frame.color = "black",
       edge.color = edge_color, edge.width = 2, main = "Consensus network",
       vertex.label.dist = 0, vertex.label.degree = 0, vertex.label.font = 2, vertex.shape = "circle")
  
  nodes <- data.frame(id = V(network_graph)$name, label = V(network_graph)$name, color = node_color)
  edges <- data.frame(from = edges$from, to = edges$to, color = edge_color)
  
  visNetwork(nodes, edges, main = "Interactive Network Plot") %>%
    visNodes(color = list(border = "black", highlight = "orange"), font = list(color = label_color)) %>%
    visEdges(color = list(color = edge_color, highlight = "orange")) %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visLayout(randomSeed = 123)
}

ui <- navbarPage(
  "CMiNet Shiny App",
  tabPanel(
    "CMiNet",
    sidebarLayout(
      sidebarPanel(
        fileInput("file", "Upload CSV File", accept = c(".csv")),
        downloadButton("downloadSampleData", "Download Sample Data"),
        checkboxInput("countBased", "Is the data count-based?", value = TRUE),
        numericInput("TT", "Quantile Threshold (TT)", value = 0.95,min=0,max=1,step=0.01),
        checkboxGroupInput(
          "selectedAlgorithms", 
          "Select Algorithms", 
          choices = c(
            "SparCC" = "sparcc",
            "SpiecEasi_MB" = "spiecEasi_mb",
            "SpiecEasi_Glasso" = "spiecEasi_glasso",
            "SPRING" = "spring",
            "GCODA" = "gcoda",
            "C_MI" = "c_MI",
            "CCLasso" = "cclasso",
            "Spearman" = "spearman",
            "Pearson" = "pearson",
            "Bicor" = "bicor"
          )
        ),
        uiOutput("algorithmParams"),
        actionButton("run", "Run CMiNet")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Data Head", tableOutput("dataSummary"), uiOutput("dataInfo")),
          tabPanel("Taxa Name", tableOutput("taxaSummary")),
          tabPanel("Results", uiOutput("resultsTab"))
        )
      )
    )
  ),
  tabPanel(
    "Visualization",
    sidebarLayout(
      sidebarPanel(
        fileInput("weightedFile", "Upload Weighted Network (CSV)", accept = c(".csv")),
        numericInput("threshold1", "Threshold 1", value = 6, min=0, max=10,step=1),
        numericInput("threshold2", "Threshold 2", value = 5, min=0, max=10,step=1),
        numericInput("threshold3", "Threshold 3", value = 4, min=0, max=10,step=1),
        numericInput("threshold4", "Threshold 4", value = 3, min=0, max=10,step=1),
        textInput("nodeColors", "Node Colors (comma-separated)", value = "white,lightyellow,lightgreen,lightblue"),
        textInput("edgeColors", "Edge Colors (comma-separated)", value = "blue,#9491D9,#332288,purple"),
        actionButton("runVisualization", "Run Visualization")
      ),
      mainPanel(
        h4("Network Visualization"),
        plotOutput("networkPlot")
      )
    )
  ),
  tabPanel(
    "Final Network",
    sidebarLayout(
      sidebarPanel(
        fileInput("finalWeightedFile", "Upload Weighted Network (CSV)", accept = c(".csv")),
        numericInput("score", "Score Threshold", value = 3, min = 0, max = 10,step=1),
        textInput("finalNodeColor", "Node Color", value = "skyblue"),
        textInput("finalEdgeColor", "Edge Color", value = "grey"),
        textInput("finalLabelColor", "Label Color", value = "black"),
        actionButton("runFinalNetwork", "Run Final Network"),
        downloadButton("downloadFinalNetworkPlot", "Save Final Network")
      ),
      mainPanel(
        h4("Final Network Plot"),
        visNetworkOutput("finalNetworkPlot"),
        #downloadButton("downloadFinalNetworkPlot", "Save Final Network Plot")
      )
    )
  ),
  tabPanel(
    "About",
    fluidPage(
      img(src = "image/logo.png", style = "width:20%; float: center;"),
      
      tags$div(
        align = "justify",
        p(
          tags$b("CMiNet: : Consensus Microbiome Network"),
          "is an R package designed to generate consensus microbiome networks by integrating results from multiple network construction algorithms. ",
          "This tool is specifically tailored for microbiome data, where capturing the intricate relationships between microbial taxa is essential to understanding complex biological systems and their impacts on health and disease."
        ),
        p(
          "The package employs a range of established algorithms, including Pearson and Spearman correlation, Biweight midcorrelation, Sparse Correlations for Compositional data (SparCC), ",
          "Sparse InversE Covariance estimation for Ecological Association and Statistical Inference (SpiecEasi), Semi-Parametric Rank-based Correlation and Partial Correlation Estimation (SPRING), ",
          "Generalized Co-Occurrence Differential Abundance analysis (gCoda), Correlation Inference for Compositional Data through Lasso (CCLasso), and a novel algorithm based on conditional mutual information (c_MI). ",
          "These algorithms construct individual microbial association networks, which CMiNet then combines into a single, weighted consensus network. By leveraging the strengths of each method, ",
          "CMiNet provides a comprehensive and reliable representation of microbial interactions."
        )
      ),
      
      img(src = "image/CMiNet-Page-2.jpg", style = "width:45%; float: right; margin-left: 10px;"),
      
      tags$h3("Algorithms Applied in CMiNet:", style = "font-size: 18px;"),
      tags$ul(
        tags$li("Pearson coefficient (cor() from stats package)"),
        tags$li("Spearman coefficient (cor() from stats package)"),
        tags$li("Biweight Midcorrelation (bicor() from WGCNA package)"),
        tags$li(a("SparCC (R code on GitHub)", href = "https://github.com/huayingfang/CCLasso/blob/master/R/SparCC.R", target = "_blank")),
        tags$li(a("CCLasso (R code on GitHub)", href = "https://github.com/huayingfang/CCLasso/tree/master", target = "_blank")),
        tags$li(a("SpiecEasi (SpiecEasi package)", href = "https://github.com/zdk123/SpiecEasi", target = "_blank")),
        tags$li(a("SPRING (SPRING package)", href = "https://github.com/GraceYoon/SPRING", target = "_blank")),
        tags$li(a("CMIMN (CMIMN package)", href = "https://github.com/solislemuslab/CMIMN", target = "_blank")),
        tags$li(a("gCoda (R code on GitHub)", href = "https://github.com/huayingfang/gCoda", target = "_blank"))
      ),
      
      tags$h3("Running CMiNet Shiny App", style = "font-size: 18px;"),
      tags$div(
        align = "justify",
        p(
          "We put the American Gut data from ",
          a("SpiecEasi package", href = "https://github.com/zdk123/SpiecEasi", target = "_blank"),
          " as an example to run. You can print the original and ID taxa names in the CMiNet page, tab data summary."
        ),
        p(
          "The CMiNet Shiny App contains four main pages:"
        ),
        tags$ol(
          tags$li(tags$b("CMiNet:"), " This page allows users to construct a consensus network from microbiome data using multiple methods. You can select the algorithms you want to use and download the results of each algorithm separately, along with the weighted network."),
          tags$li(tags$b("Visualization:"), " This page processes a weighted microbiome network and visualizes it across different thresholds. Each threshold represents a minimum edge weight required for inclusion in the network plot."),
          tags$li(tags$b("Final Network:"), " This page illustrates the final network produced by CMiNet based on the threshold defined by the user."),
          tags$li(tags$b("About:"), " Contains detailed information about running and interpreting results from CMiNet."),
          tags$li(tags$b("GitHub:"), a("CMiNet Algorihtm", href ="https://github.com/solislemuslab/CMiNet", target = "_blank"))
        ),
        tags$ul(
          tags$li(
            a("Submit your issues and comments on the GitHub page", 
              href = "https://github.com/solislemuslab/CMiNet", target = "_blank")
          )
        )),
      tags$h3("Citation", style = "font-size: 18px; margin-bottom: 10px;"),
      tags$p(
        "If you use CMiNet in your work, we kindly ask that you cite the following papers:",
        style = "font-size: 14px; margin-bottom: 8px;"
      ),
      tags$div(
        align = "justify",
        p(
          "[1] Aghdam R, Solis-Lemus C. CMiNet: R package for learning the Consensus Microbiome Network. arXiv preprint arXiv:2411.08309. 2024. "),
        p(
          "[2] Aghdam R, Tang X, Shan S, Lankau R, Solís-Lemus C. Human limits in Machine Learning: Prediction of plant phenotypes using soil microbiome data. 2024."
        )
      )
      
    )
  )
  
  
)


shinyApp(ui = ui, server = server)
