library(shiny)
library(SpiecEasi)
library(shinyWidgets)
library(shinyBS)
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
  
  algorithm_links <- list(
    SparCC = "https://github.com/huayingfang/CCLasso/blob/master/R/SparCC.R",
    spiecEasi_mb = "https://github.com/zdk123/SpiecEasi",
    spiecEasi_glasso = "https://github.com/zdk123/SpiecEasi",
    spring = "https://github.com/GraceYoon/SPRING",
    gcoda = "https://github.com/huayingfang/gCoda",
    CMIMN = "https://github.com/solislemuslab/CMIMN",
    cclasso = "https://github.com/huayingfang/CCLasso/tree/master",
    spearman = "https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/cor",
    pearson = "https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/cor",
    bicor = "https://www.rdocumentation.org/packages/WGCNA/versions/1.70-3/topics/bicor"
  )
  
  algorithm_display_names <- list(
    SparCC = "SparCC",
    spiecEasi_mb = "SpiecEasi_MB",
    spiecEasi_glasso = "SpiecEasi_Glasso",
    spring = "SPRING",
    gcoda = "GCODA",
    CMIMN = "CMIMN",
    cclasso = "CCLasso",
    spearman = "Spearman",
    pearson = "Pearson",
    bicor = "Bicor"
  )
  algorithmParams <- list(
    SparCC = list(
      numericInput("imax", "Outer Iterations (imax)", value = 20, min = 1, max = 100, step = 1),
      bsTooltip("imax", "Number of iterations in the outer loop (default: 20).", placement = "right", trigger = "hover"),
      
      numericInput("kmax", "Inner Iterations (kmax)", value = 10, min = 1, max = 100, step = 1),
      bsTooltip("kmax", "Number of iterations in the inner loop (default: 10).", placement = "right", trigger = "hover"),
      
      numericInput("alpha", "Threshold (alpha)", value = 0.1, min = 0, max = 1, step = 0.01),
      bsTooltip("alpha", "Absolute value of correlations below this threshold are considered zero (default: 0.1).", placement = "right", trigger = "hover"),
      
      numericInput("Vmin", "Minimum Variance (Vmin)", value = 1e-4, min = 1e-4, max = 1, step = 0.01),
      bsTooltip("Vmin", "Minimum variance threshold (default: 1e-4).", placement = "right", trigger = "hover")
    ),
    spiecEasi_mb = list(
      numericInput("lambda_min_ratio_mb", "Lambda Min Ratio", value = 1e-2, min = 0, max = 1, step = 0.01),
      bsTooltip("lambda_min_ratio_mb", "Minimum ratio of the penalty parameter (default: 1e-2).", placement = "right", trigger = "hover"),
      
      numericInput("nlambda_mb", "Number of Lambdas", value = 15, min = 1, max = 100, step = 1),
      bsTooltip("nlambda_mb", "Number of tuning parameters (default: 15).", placement = "right", trigger = "hover"),
      
      numericInput("rep_num_mb", "Replications (rep_num)", value = 20, min = 1, max = 100, step = 1),
      bsTooltip("rep_num_mb", "Number of replications for stability selection (default: 20).", placement = "right", trigger = "hover"),
      
      numericInput("ncores_mb", "Number of Cores", value = 4, min = 1, max = 20, step = 1),
      bsTooltip("ncores_mb", "Number of CPU cores to use for parallel computation (default: 4).", placement = "right", trigger = "hover")
    ),
    spiecEasi_glasso = list(
      numericInput("lambda_min_ratio_glasso", "Lambda Min Ratio", value = 1e-2, min = 0, max = 1, step = 0.01),
      bsTooltip("lambda_min_ratio_glasso", "Minimum ratio of the penalty parameter (default: 1e-2).", placement = "right", trigger = "hover"),
      
      numericInput("nlambda_glasso", "Number of Lambdas", value = 15, min = 1, max = 100, step = 1),
      bsTooltip("nlambda_glasso", "Number of tuning parameters (default: 15).", placement = "right", trigger = "hover"),
      
      numericInput("rep_num_glasso", "Replications (rep_num)", value = 50, min = 1, max = 100, step = 1),
      bsTooltip("rep_num_glasso", "Number of replications for stability selection (default: 50).", placement = "right", trigger = "hover"),
      numericInput("ncores_glasso", "Number of Cores", value = 4, min = 1, max = 20, step = 1),
      bsTooltip("ncores_glasso", "Number of CPU cores to use for parallel computation (default: 4).", placement = "right", trigger = "hover")
    ),
    spring = list(
      numericInput("ncores_spring", "Number of Cores", value = 5, min = 1, max = 20, step = 1),
      bsTooltip("ncores_spring", "Number of CPU cores to use for parallel computation (default: 5).", placement = "right", trigger = "hover"),
      
      numericInput("nlambda_spring", "Number of Lambdas", value = 15, min = 1, max = 100, step = 1),
      bsTooltip("nlambda_spring", "Number of tuning parameters (default: 15).", placement = "right", trigger = "hover"),
      
      numericInput("rep_num_spring", "Replications (rep_num)", value = 20, min = 1, max = 100, step = 1),
      bsTooltip("rep_num_spring", "Number of replications for stability selection (default: 20).", placement = "right", trigger = "hover")
    ),
    gcoda = list(
      numericInput("pseudo_gcoda", "Pseudo Count", value = 0.5, min = 0, max = 2, step = 0.1),
      bsTooltip("pseudo_gcoda", "Pseudo count to add if the input is count data (default: 0.5).", placement = "right", trigger = "hover"),
      
      numericInput("lambda_min_ratio_gcoda", "Lambda Min Ratio", value = 1e-4),
      bsTooltip("lambda_min_ratio_gcoda", "Minimum ratio of the penalty parameter (default: 1e-4).", placement = "right", trigger = "hover"),
      
      numericInput("nlambda_gcoda", "Number of Lambdas", value = 15, min = 1, max = 100, step = 1),
      bsTooltip("nlambda_gcoda", "Number of tuning parameters (default: 15).", placement = "right", trigger = "hover"),
      
      numericInput("ebic_gamma_gcoda", "EBIC Gamma", value = 0.5, min = 0, max = 1, step = 0.01),
      bsTooltip("ebic_gamma_gcoda", "Gamma parameter for the EBIC criterion (default: 0.5).", placement = "right", trigger = "hover")
    ),
    CMIMN = list(
      numericInput("q1", "Q1 Threshold", value = 0.7, min = 0, max = 1, step = 0.1),
      bsTooltip("q1", "Quantile threshold for the Mutual Information test (default: 0.7).", 
                placement = "right", trigger = "hover"),
      
      numericInput("q2", "Q2 Threshold", value = 0.95, min = 0, max = 1, step = 0.1),
      bsTooltip("q2", "Quantile threshold for the Conditional Mutual Information test (default: 0.95).", 
                placement = "right", trigger = "hover")
    ),
    
    cclasso = list(
      numericInput("pseudo_cclasso", "Pseudo Count", value = 0.5, min = 0, max = 2, step = 0.1),
      bsTooltip("pseudo_cclasso", "Pseudo count added to avoid zeroes in the compositional data (default: 0.5).", 
                placement = "right", trigger = "hover"),
      
      numericInput("k_cv_cclasso", "Cross-Validation Folds (K CV)", value = 3, min = 1, max = 10, step = 1),
      bsTooltip("k_cv_cclasso", "Number of folds for cross-validation to select the optimal tuning parameter (default: 3).", 
                placement = "right", trigger = "hover"),
      
      numericInput("k_max_cclasso", "Maximum Iterations (K Max)", value = 300, min = 1, max = 400, step = 10),
      bsTooltip("k_max_cclasso", "Maximum number of iterations for the optimization method (default: 300).", 
                placement = "right", trigger = "hover"),
      
      numericInput("n_boot_cclasso", "Bootstrap Samples", value = 20, min = 1, max = 200, step = 10),
      bsTooltip("n_boot_cclasso", "Number of bootstrap samples used for model selection (default: 20).", 
                placement = "right", trigger = "hover")
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
        link <- algorithm_links[[algorithm]] # Get the corresponding link
        display_name <- algorithm_display_names[[algorithm]] # Get the display name
        
        if (length(params) > 0) {
          tags$fieldset(
            tags$legend(
              tags$a(
                href = link,
                target = "_blank",
                style = "color: #007bff; text-decoration: underline;",
                display_name # Use the display name
              )
            ),
            do.call(tagList, params)
          )
        } else {
          tags$fieldset(
            tags$legend(
              tags$a(
                href = link,
                target = "_blank",
                style = "color: #007bff; text-decoration: underline;",
                display_name # Use the display name
              )
            ),
            p("No parameters to configure for this algorithm.")
          )
        }
      })
    )
  })
  
  output$dataSummary <- renderTable({
    req(data())
    head(data()[, 1:10])
  })
  
  output$dataInfo <- renderUI({
    req(data())
    df <- data()
    num_samples <- nrow(df)
    num_taxa <- ncol(df)
    
    tagList(
      p("Your data should be formatted as follows:"),
      p("Rows represent samples, and columns represent taxa."),
      p(paste("Number of taxa =", num_taxa)),
      p(paste("Number of samples =", num_samples))
    )
  })
  
  results <- reactiveVal(NULL)
  observeEvent(input$run, {
    withProgress(message = 'Running CMiNet...', value = 0, {
      # Ensure all prerequisites are met
      req(data())  # Data must be uploaded
      req(input$selectedAlgorithms)
      
      # Clear previous results
      unlink("Network/*", recursive = TRUE)
      unlink("Binary_Network/*", recursive = TRUE)
      
      # Increment progress
      incProgress(0.1, detail = "Preparing input data")

      isCount <- input$countBased
      TT <- as.numeric(input$TT)
      dir.create("Network", showWarnings = FALSE)
      dir.create("Binary_Network", showWarnings = FALSE)

      pearson_params <- if ("pearson" %in% input$selectedAlgorithms) {
        list(enabled = TRUE, params = list())
      } else {
        list(enabled = FALSE, params = list())
      }
      
      spearman_params <- if ("spearman" %in% input$selectedAlgorithms) {
        list(enabled = TRUE, params = list())
      } else {
        list(enabled = FALSE, params = list())
      }
      
      bicor_params <- if ("bicor" %in% input$selectedAlgorithms) {
        list(enabled = TRUE, params = list())
      } else {
        list(enabled = FALSE, params = list())
      }
      
      CMIMN_params <- if ("CMIMN" %in% input$selectedAlgorithms) {
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
      
      gcoda_params <- if ("gcoda" %in% input$selectedAlgorithms) {
        list(
          enabled = TRUE,
          params = list(
            counts = isCount,
            pseudo = input$pseudo_gcoda,
            lambda.min.ratio = input$lambda_min_ratio_gcoda,
            nlambda = input$nlambda_gcoda,
            ebic.gamma = input$ebic_gamma_gcoda
          )
        )
      } else {
        list(enabled = FALSE, params = list())
      }
      
      cclasso_params <- if ("cclasso" %in% input$selectedAlgorithms) {
        list(
          enabled = TRUE,
          params = list(
            counts = isCount,
            pseudo = input$pseudo_cclasso,
            k_cv = input$k_cv_cclasso,
            lam_int = c(1e-4, 1),
            k_max = input$k_max_cclasso,
            n_boot = input$n_boot_cclasso
          )
        )
      } else {
        list(enabled = FALSE, params = list())
      }
      
      spiecEasi_mb_params <- if ("spiecEasi_mb" %in% input$selectedAlgorithms) {
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
      
      spiecEasi_glasso_params <- if ("spiecEasi_glasso" %in% input$selectedAlgorithms) {
        list(
          enabled = TRUE,
          params = list(
            method = 'glasso',
            lambda.min.ratio = input$lambda_min_ratio_glasso,
            pulsar.params = list(rep.num = input$rep_num_glasso, ncores = input$ncores_glasso)
          )
        )
      } else {
        list(enabled = FALSE, params = list())
      }
      
      sparcc_params <- if ("SparCC" %in% input$selectedAlgorithms) {
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
      
      spring_params <- if ("spring" %in% input$selectedAlgorithms) {
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
      
      # Increment progress
      incProgress(0.3, detail = "Running CMiNet algorithm")
      
      # Run the CMiNet algorithm
      result <- CMiNet(
        data = data(),
        quantitative = isCount,
        TT = TT,
        pearson = pearson_params,
        spearman = spearman_params,
        bicor = bicor_params,
        sparcc = sparcc_params,
        spiecEasi_mb = spiecEasi_mb_params,
        spiecEasi_glasso = spiecEasi_glasso_params,
        spring = spring_params,
        gcoda = gcoda_params,
        c_MI = CMIMN_params,
        cclasso = cclasso_params
      )
      
      # Increment progress
      incProgress(0.9, detail = "Saving results")
      
      # Save results to files
      lapply(names(result), function(algorithm) {
        if (is.list(result[[algorithm]]) && !is.null(result[[algorithm]]$binary_matrix)) {
          write.csv(result[[algorithm]]$binary_matrix, file.path("Binary_Network", paste0(algorithm, "_binary.csv")))
          write.csv(result[[algorithm]]$adj_matrix, file.path("Network", paste0(algorithm, ".csv")))
        }
      })
      # Save weighted network to CSV
      if (!is.null(result$weighted_network)) {
        weighted_network <- as.matrix(result$weighted_network)
        write.csv(weighted_network, file.path("Network", "weighted_network.csv"))
      }
      
      if (!is.null(result$edge_list) && nrow(result$edge_list) > 0) {
        edge_list <- as.matrix(result$edge_list)
        write.csv(edge_list, file.path("Network", "edge_list.csv"), row.names = FALSE)
        print("Edge list saved successfully.")
      } else {
        print("Edge list is empty or NULL. Skipping save.")
      }
      
      results(result)
      
      # Increment progress
      incProgress(1, detail = "Completed!")
    })
  })
  
  output$resultsTab <- renderUI({
    req(results())  # Ensure results are available
    
    # Check if edge_list exists and is non-empty
    edge_list_exists <- !is.null(results()$edge_list) && nrow(results()$edge_list) > 0
    
    tagList(
      tags$div(
        downloadButton("downloadWeightedNetwork", "Download Weighted Network"),
        tags$p("This CSV file shows the resulting weighted network based on the selected algorithms. The edge values in the network range between 0 and the maximum number of selected algorithms and it represents the number of algorithms that infer the presence of this edge.")
      ),

      if (edge_list_exists) {
        tags$div(
          downloadButton("downloadEdgeList", "Download Edge List"),
          tags$p("This CSV file contains the resulting edge list. It has four columns: the first and second columns represent the names of the connected nodes, the third column indicates the weight of the edge between them (number of algorithms that inferred the presence of the edge), and the fourth column specifies the methods that inferred the edge.")
        )
      },
      tags$div(
        downloadButton("downloadBinaryFolder", "Download Binary Folder"),
        tags$p("This folder contains the binary adjacency matrices of all selected algorithms as CSV files.")
      ),
      tags$div(
        downloadButton("downloadNetworkFolder", "Download Network Folder"),
        tags$p("This folder contains the results of all selected algorithms (the continuous values before binarization) along with an error file and the weighted network CSV file.")
      )
    )
  })
  
  output$downloadBinaryFolder <- downloadHandler(
    filename = function() {
      paste("Binary_Network_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      binary_files <- list.files("Binary_Network", full.names = TRUE)
      lapply(binary_files, function(binary_file) {
        binary_matrix <- as.matrix(read.csv(binary_file, row.names = 1))
        
        # Modify the matrix: zero the diagonal
        diag(binary_matrix) <- 0
        
        # Update row/column names (example: using taxa_name reactive output)
        taxa_names <- taxa_name()[, 1]
        if (length(taxa_names) == nrow(binary_matrix)) {
          rownames(binary_matrix) <- taxa_names
          colnames(binary_matrix) <- taxa_names
        }
        
        write.csv(binary_matrix, binary_file, row.names = TRUE)
      })
      
      # Zip the updated binary matrices
      zip::zip(file, binary_files)
    }
  )

  output$downloadNetworkFolder <- downloadHandler(
    filename = function() {
      paste("Network_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      # Get the list of adjacency matrix files
      network_files <- list.files("Network", full.names = TRUE)
      
      # Exclude files that aren't adjacency matrices (e.g., weighted network)
      network_files <- network_files[!grepl("weighted_network", network_files)]
      network_files <- network_files[!grepl("edge_list", network_files)]
      network_files <- network_files[!grepl("errors_warnings", network_files)]
      # Read, modify, and overwrite each adjacency matrix
      lapply(network_files, function(network_file) {
        adj_matrix <- as.matrix(read.csv(network_file, row.names = 1))
        diag(adj_matrix) <- 0

        taxa_names <- taxa_name()[, 1]
        if (length(taxa_names) == nrow(adj_matrix)) {
          rownames(adj_matrix) <- taxa_names
          colnames(adj_matrix) <- taxa_names
        }
        write.csv(adj_matrix, network_file, row.names = TRUE)
      })
      zip::zip(file, network_files)
    }
  )
  
  
  output$downloadWeightedNetwork <- downloadHandler(
    filename = function() {
      paste("weighted_network_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      file_path <- file.path("Network", "weighted_network.csv")
      if (file.exists(file_path)) {
        weighted_network <- as.matrix(read.csv(file_path, row.names = 1))
        
        diag(weighted_network) <- 0
        taxa_names <- taxa_name()[, 1]
        if (length(taxa_names) == nrow(weighted_network)) {
          rownames(weighted_network) <- taxa_names
          colnames(weighted_network) <- taxa_names
        }
        
        # Save the modified weighted network
        write.csv(weighted_network, file, row.names = TRUE)
      } else {
        stop("The weighted network file does not exist.")
      }
    }
  )
  
  output$downloadEdgeList <- downloadHandler(
    filename = function() {
      paste("edge_list_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      file_path <- file.path("Network", "edge_list.csv")
      if (file.exists(file_path) && file.info(file_path)$size > 0) {
        edge_list <- read.csv(file_path, check.names = FALSE)
        write.csv(edge_list, file, row.names = FALSE)
      } else {
        stop("The edge list file does not exist or is empty.")
      }
    }
  )
  
  output$downloadSampleData <- downloadHandler(
    filename = function() {
      "sample_data.csv"  # The name of the file the user will download
    },
    content = function(file) {
      file.copy("sample/sample_data.csv", file)
    }
  )
  
  output$downloadSampleNet<- downloadHandler(
    filename = function() {
      "weighted_network.csv"  # The name of the file the user will download
    },
    content = function(file) {
      file.copy("sample/weighted_network.csv", file)
    }
  ) 

  observeEvent({
    input$threshold1
    input$threshold2
    input$threshold3
    input$threshold4
    input$runVisualization
  }, {
    #req(input$weightedFile)
    weighted_network <- tryCatch({
      if (input$fileOption == "Result of CMiNet") {
        file_path <- file.path("Network", "weighted_network.csv")
        if (!file.exists(file_path)) {
          showNotification("CMiNet result file (weighted_network.csv) not found.", type = "error")
          return(NULL)
        }
        network <- as.matrix(read.csv(file_path, row.names = 1, check.names = FALSE))
      } else if (input$fileOption == "Upload Weighted Network (CSV file)") {
        req(input$weightedFile)
        network <- as.matrix(read.csv(input$weightedFile$datapath, row.names = 1, check.names = FALSE))
        # Ensure the network is numeric
      } else {
        showNotification("Please select a valid file option.", type = "error")
        return(NULL)
      }
    }, error = function(e) {
      showNotification("Error reading weighted network file.", type = "error")
      return(NULL)
    })
    
    req(weighted_network)
    
    # Ensure the weighted network is square
    if (nrow(weighted_network) != ncol(weighted_network)) {
      showNotification("Error: Weighted network matrix must be square", type = "error")
      return(NULL)
    }
    
    # Get the input parameters
    thresholds <- c(input$threshold1, input$threshold2, input$threshold3, input$threshold4)
    show_labels <- c(1 %in% input$showLabels, 2 %in% input$showLabels, 3 %in% input$showLabels, 4 %in% input$showLabels)
    node_colors <- unlist(strsplit(input$nodeColors, ","))
    edge_colors <- unlist(strsplit(input$edgeColors, ","))
    
    # Render the plot dynamically
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
  
  final_weighted_network <- reactiveVal(NULL)
  observeEvent(input$runFinalNetwork, {
    # Load the weighted network based on the user's selection
    final_weighted_network <- tryCatch({
      if (input$finalFileOption == "Result of CMiNet") {
        # Load CMiNet-generated weighted network
        file_path <- file.path("Network", "weighted_network.csv")
        if (!file.exists(file_path)) {
          showNotification("CMiNet result file (weighted_network.csv) not found.", type = "error")
          return(NULL)
        }
        as.matrix(read.csv(file_path, row.names = 1, check.names = FALSE))
      } else if (input$finalFileOption == "Upload Weighted Network (CSV file)") {
        req(input$finalWeightedFile)
        as.matrix(read.csv(input$finalWeightedFile$datapath, row.names = 1, check.names = FALSE))
      } else {
        showNotification("Please select a valid file option.", type = "error")
        return(NULL)
      }
    }, error = function(e) {
      showNotification("Error reading weighted network file.", type = "error")
      return(NULL)
    })
    
    if (nrow(final_weighted_network) != ncol(final_weighted_network)) {
      showNotification("Error: Weighted network matrix must be square.", type = "error")
      return(NULL)
    }

    output$finalNetworkPlot <- renderVisNetwork({
      score_threshold <- input$score
      final_network <- ifelse(final_weighted_network > score_threshold, 1, 0)
      final_network[lower.tri(final_network)] <- 0
      
      # Create edge and node data for the plot
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

      visNetwork(nodes, edges, main = "") %>%
        visNodes(color = list(background = input$finalNodeColor, border = "black"), font = list(color = input$finalLabelColor)) %>%
        visEdges(color = list(color = input$finalEdgeColor, highlight = "orange")) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visLayout(randomSeed = 123)
    })
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
    
    V(graph)$color <- node_colors[i]
    E(graph)$color <- edge_colors[i]
    
    set.seed(123) 
    plot(
      graph,
      vertex.label = ifelse(show_labels[i], taxa_names, NA),
      vertex.size = 5,
      main = paste("Threshold =", threshold))
    mtext(side = 1, line = 4, paste("NN Deg >0: ", sum(degree(graph) > 0), ", N Edges: ", ecount(graph), ", Max-Deg: ", max(degree(graph))),
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
        fileInput("file", "Upload Microbiome Data", accept = c(".csv")),
        downloadButton("downloadSampleData", "Sample Data"),
        tags$i(
          class = "fa fa-question-circle",
          style = "color: blue; margin-left: 5px; cursor: pointer;",
          title = "Download a Sample Gut Microbiome Dataset from the American Gut Project to explore the data format."
        ),
        checkboxInput("countBased", "Is the data count-based?", value = TRUE),
        numericInputIcon(
          inputId = "TT",
          label = tags$span(
            "Quantile Threshold",
            tags$i(
              class = "fa fa-question-circle",
              style = "color: blue; margin-left: 5px; cursor: pointer;",
              title = "The threshold value is used for threshold-dependent algorithms (Pearson, Spearman, Bicor, SparCC, and CCLasso). By default, it is set to the 0.95 quantile to promote sparsity in the resulted network. For these algorithms, values below this quantile are binarized to 0, while values equal to or above the quantile are binarized to 1."
            )
          ),
          value = 0.95,
          min = 0,
          max = 1,
          step = 0.01,
          icon = icon("adjust", style = "color: blue;")
        ),
        checkboxGroupInput(
          "selectedAlgorithms", 
          "Select Algorithms", 
          choices = c(
            "SparCC" = "SparCC",
            "SpiecEasi_MB" = "spiecEasi_mb",
            "SpiecEasi_Glasso" = "spiecEasi_glasso",
            "SPRING" = "spring",
            "GCODA" = "gcoda",
            "CMIMN" = "CMIMN",
            "CCLasso" = "cclasso",
            "Spearman" = "spearman",
            "Pearson" = "pearson",
            "Bicor" = "bicor"
          )
        ),
        bsTooltip("selectedAlgorithms", "Select the algorithms you want to include in the weighted network construction.", 
                  placement = "right", trigger = "hover"),
        
        uiOutput("algorithmParams"),
        
        actionButton(
          inputId = "run",
          label = "Run CMiNet",
          style = "background-color: #007bff; color: white; font-weight: bold; border-radius: 5px; padding: 10px 20px; border: none;"
        ),
        bsTooltip("run", "Click to run the selected algorithms and generate a weighted network", 
                  placement = "right", trigger = "hover")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Data Head", tableOutput("dataSummary"), uiOutput("dataInfo")),
          #tabPanel("Taxa Name", tableOutput("taxaSummary")),
          tabPanel("Results", uiOutput("resultsTab"))
        )
      )
    )
  ),

  tabPanel(
    "Visualization",
    sidebarLayout(
      
      
      sidebarPanel(
        # Combined Dropdown for File Selection
        tags$div(
          style = "margin-bottom: 10px;",
          selectInput(
            inputId = "fileOption",
            label = "Load Weighted Network",
            choices = c( "Result of CMiNet", "Upload Weighted Network (CSV file)"),
            selected = "Result of CMiNet")
        ),
        
        # File Input for Upload (Initially Hidden)
        conditionalPanel(
          condition = "input.fileOption == 'Upload Weighted Network (CSV file)'",
          fileInput("weightedFile", "Upload Weighted Network", accept = c(".csv"))
        ),
        
        # Tooltip for CMiNet Result
        conditionalPanel(
          condition = "input.fileOption == 'Result of CMiNet'",
          tags$div(
            style = "margin-top: 10px; color: blue;",
            tags$b("The CMiNet result (weighted_network.csv) obtained from the analysis in the CMiNet tab will be loaded automatically.")
          )
        ),
        
        downloadButton("downloadSampleNet", "Sample Weighted Network"),
        tags$i(
          class = "fa fa-question-circle",
          style = "color: blue; margin-left: 5px; cursor: pointer;",
          title = "Download the Weighted Gut Microbiome Network Generated by CMiNet Using Sample Data from the American Gut Project."
        ),
        
        # Threshold 1
        numericInputIcon(
          inputId = "threshold1",
          label = tags$span(
            "Threshold 1",
            tags$i(
              class = "fa fa-question-circle",
              style = "color: blue; margin-left: 5px; cursor: pointer;",
              title = "Edges remain in the network if their weight values are higher than this threshold."
            )
          ),
          value = 7,
          min = 0,
          max = 10,
          step = 1,
          icon = icon("adjust", style = "color: blue;")
        ),

        numericInputIcon(
          inputId = "threshold2",
          label = tags$span(
            "Threshold 2",
            tags$i(
              class = "fa fa-question-circle",
              style = "color: green; margin-left: 5px; cursor: pointer;",
              title = "Edges remain in the network if their weight values are higher than this threshold."
            )
          ),
          value = 5,
          min = 0,
          max = 10,
          step = 1,
          icon = icon("adjust", style = "color: green;")
        ),

        numericInputIcon(
          inputId = "threshold3",
          label = tags$span(
            "Threshold 3",
            tags$i(
              class = "fa fa-question-circle",
              style = "color: orange; margin-left: 5px; cursor: pointer;",
              title = "Edges remain in the network if their weight values are higher than this threshold."
            )
          ),
          value = 4,
          min = 0,
          max = 10,
          step = 1,
          icon = icon("adjust", style = "color: orange;")
        ),
        
        # Threshold 4
        numericInputIcon(
          inputId = "threshold4",
          label = tags$span(
            "Threshold 4",
            tags$i(
              class = "fa fa-question-circle",
              style = "color: purple; margin-left: 5px; cursor: pointer;",
              title = "Edges remain in the network if their weight values are higher than this threshold."
            )
          ),
          value = 3,
          min = 0,
          max = 10,
          step = 1,
          icon = icon("adjust", style = "color: purple;")
        ),
        
        # Node Colors
        tags$div(
          style = "display: flex; align-items: center; margin-bottom: 10px;",
          textInput("nodeColors", "Node Colors (comma-separated)", value = "white,lightyellow,lightgreen,lightblue"),
          tags$i(
            class = "fa fa-palette",
            style = "color: blue; margin-left: 5px; cursor: pointer;",
            title = "Define colors for the nodes in the network, separated by commas (e.g., white, yellow, green, blue)."
          )
        ),

        tags$div(
          style = "display: flex; align-items: center;",
          textInput("edgeColors", "Edge Colors (comma-separated)", value = "blue,#9491D9,#332288,purple"),
          tags$i(
            class = "fa fa-brush",
            style = "color: purple; margin-left: 5px; cursor: pointer;",
            title = "Define colors for the edges in the network, separated by commas (e.g., blue, red, green, purple)."
          )
        ),
        
        # Run Button
        #actionButton("runVisualization", "Run Visualization", style = "background-color: #007bff; color: white;"),
        
        actionButton(
          inputId = "runVisualization",
          label = "Visualization",
          style = "background-color: #007bff; color: white; font-weight: bold; border-radius: 5px; padding: 10px 20px; border: none;"
        ),
        bsTooltip("runVisualization", "Click to generate four consensus networks based on different thresholds and choose the final threshold for the final network.", 
                  placement = "right", trigger = "hover")
      ),
      mainPanel(
        h4("Network Visualization Based on Four Thresholds"),
        p("This visualization compares network structures based on four different thresholds. Below each network depicts the number of nodes with degree > 0 (NN Deg >0), the number of edges (N Edges), and the maximum degree value in the network (Max-Deg)."),
        plotOutput("networkPlot")
      )
    )
  ),

  tabPanel(
    "Final Network",
    sidebarLayout(
      sidebarPanel(
        # File Selection Option
        tags$div(
          style = "margin-bottom: 10px;",
          selectInput(
            inputId = "finalFileOption",
            label = "Load Weighted Network",
            choices = c("Result of CMiNet", "Upload Weighted Network (CSV file)"),
            selected = "Result of CMiNet"
          )
        ),
        conditionalPanel(
          condition = "input.finalFileOption == 'Upload Weighted Network (CSV file)'",
          fileInput("finalWeightedFile", "Upload Weighted Network", accept = c(".csv"))
        ),
        conditionalPanel(
          condition = "input.finalFileOption == 'Result of CMiNet'",
          tags$div(
            style = "margin-top: 10px; color: blue;",
            tags$b("The CMiNet result (weighted_network.csv) obtained from the analysis in the CMiNet tab will be loaded automatically.")
          )
        ),
        numericInputIcon(
          inputId = "score",
          label = tags$span(
            "Threshold",
            tags$i(
              class = "fa fa-question-circle",
              style = "color: blue; margin-left: 5px; cursor: pointer;",
              title = "Edges remain in the network if their weight values are higher than this threshold."
            )
          ),
          value = 6,
          min = 0,
          max = 9,
          step = 1,
          icon = icon("adjust", style = "color: blue;")
        ),
        
        # Node Color Input with Tooltip
        tags$div(
          style = "display: flex; align-items: center; margin-bottom: 10px;",
          textInput("finalNodeColor", "Node Color", value = "skyblue"),
          tags$i(
            class = "fa fa-palette",
            style = "color: blue; margin-left: 5px; cursor: pointer;",
            title = "Select a color for the nodes in the network."
          )
        ),
        
        # Edge Color Input with Tooltip
        tags$div(
          style = "display: flex; align-items: center; margin-bottom: 10px;",
          textInput("finalEdgeColor", "Edge Color", value = "grey"),
          tags$i(
            class = "fa fa-brush",
            style = "color: purple; margin-left: 5px; cursor: pointer;",
            title = "Select a color for the edges in the network."
          )
        ),
        
        # Label Color Input with Tooltip
        tags$div(
          style = "display: flex; align-items: center;",
          textInput("finalLabelColor", "Label Color", value = "black"),
          tags$i(
            class = "fa fa-font",
            style = "color: black; margin-left: 5px; cursor: pointer;",
            title = "Select a color for the labels in the network."
          )
        ),
        
        # Run Button
        #actionButton("runFinalNetwork", "Run Final Network", style = "background-color: #007bff; color: white;")
        actionButton(
          inputId = "runFinalNetwork",
          label = "Generate Final Network",
          style = "background-color: #007bff; color: white; font-weight: bold; border-radius: 5px; padding: 10px 20px; border: none;"
        ),
        bsTooltip("runFinalNetwork", "Click to generate the final network using the selected threshold.", 
                  placement = "right", trigger = "hover")
        
      ),
      mainPanel(
        h4("Final Microbiome Network"),
        visNetworkOutput("finalNetworkPlot")
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
          "Generalized Co-Occurrence Differential Abundance analysis (gCoda), Correlation Inference for Compositional Data through Lasso (CCLasso), and a novel algorithm based on conditional mutual information (CMIMN). ",
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
          "[2] Aghdam, R., Tang, X., Shan, S., Lankau, R., & Solís-Lemus, C. (2024). Human limits in machine learning: prediction of potato yield and disease using soil microbiome data. BMC bioinformatics, 25, 366."
        )
      )
      
    )
  )
)

shinyApp(ui = ui, server = server)
