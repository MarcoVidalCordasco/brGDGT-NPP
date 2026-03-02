#' Diagnose if new brGDGT data falls within the training data space
#'
#' @description
#' Checks whether new brGDGT samples fall within the brGDGT composition space
#' of the training data used to build the prediction model. This helps identify
#' samples that might yield unreliable predictions. If geographic coordinates
#' are provided, it also shows the spatial distribution of similar samples.
#'
#' @param brgdt_data A data frame containing brGDGT measurements. Must include:
#'   * First column: Sample ID or identifier
#'   * Columns with brGDGT compounds: fIa, fIb, fIc, fIIa, fIIa_, fIIb, fIIb_,
#'     fIIc, fIIc_, fIIIa, fIIIa_, fIIIb, fIIIb_, fIIIc, fIIIc_
#'   * Optional: Latitude and Longitude columns for geographic mapping
#' @param reference_data Optional reference data. If NULL, loads from package
#' @param plot Logical; if TRUE, generates diagnostic plots (default = TRUE)
#' @param threshold Numeric; Mahalanobis distance threshold quantile (default = 0.95)
#' @param use_pca Logical; if TRUE, uses PCA to avoid singular matrix (default = TRUE)
#' @param map Logical; if TRUE, generates geographic map (requires Latitude/Longitude)
#' @param n_neighbors Integer; number of nearest neighbors to identify (default = 5)
#'
#' @return A list with diagnostic results:
#'   \item{within_domain}{Logical vector indicating if each sample is within domain}
#'   \item{mahalanobis_dist}{Mahalanobis distances to training data}
#'   \item{p_values}{P-values for each sample}
#'   \item{nearest_neighbor_dist}{Distance to nearest training sample}
#'   \item{nearest_neighbors}{Data frame with nearest neighbor information}
#'   \item{geographic}{Geographic information (if coordinates provided)}
#'   \item{plots}{List of ggplot objects (if plot=TRUE)}
#'
#' @importFrom compositions clr
#' @importFrom ggplot2 ggplot aes geom_point geom_hline labs theme_minimal
#' @importFrom ggplot2 scale_color_manual scale_fill_manual scale_size_continuous
#' @importFrom ggplot2 geom_polygon geom_text coord_fixed element_text
#' @importFrom stats mahalanobis qchisq pchisq quantile cor
#' @importFrom grDevices rgb
#' @importFrom geosphere distHaversine
#' @importFrom readxl read_excel
#' @importFrom plotly plot_ly layout add_trace
#' @export
#'
#' @examples
#' \dontrun{
#' # Load your brGDGT data
#' my_data <- read.csv("my_brgdt_data.csv")
#' 
#' # Check if data falls within model domain
#' diagnostics <- diagnose_applicability(my_data)
#' 
#' # View which samples are safe to predict
#' print(diagnostics$within_domain)
#' 
#' # Generate diagnostic plots with map (if coordinates available)
#' diagnostics <- diagnose_applicability(my_data, plot = TRUE, map = TRUE)
#' }
diagnose_applicability <- function(brgdt_data, reference_data = NULL, 
                                   plot = TRUE, threshold = 0.95,
                                   use_pca = TRUE, map = TRUE,
                                   n_neighbors = 5) {
  
  # Load reference data if not provided
  if(is.null(reference_data)) {
    # Try to load from Excel file first
    excel_path <- system.file("extdata", "Dataset_.xlsx", package = "brGDGTNPP")
    
    if(file.exists(excel_path)) {
      message("Analyses in progress...")
      
      # Read the Dataset sheet
      raw_data <- readxl::read_excel(excel_path, sheet = "Dataset")
      
      # Define brGDGT columns
      brGDGT_cols <- c("fIa", "fIb", "fIc", "fIIa", "fIIa_", "fIIb", "fIIb_",
                       "fIIc", "fIIc_", "fIIIa", "fIIIa_", "fIIIb", "fIIIb_",
                       "fIIIc", "fIIIc_")
      
      # Extract NPP if available
      npp_col <- NULL
      if("NPP" %in% colnames(raw_data)) {
        npp_col <- raw_data$NPP
      } else if("NPP (gC/m²/yr)" %in% colnames(raw_data)) {
        npp_col <- raw_data$`NPP (gC/m²/yr)`
      }
      
      # Extract coordinates if available
      coords <- NULL
      if(all(c("Latitude", "Longitude") %in% colnames(raw_data))) {
        coords <- raw_data[, c("Latitude", "Longitude")]
      }
      
      # Calculate CLR values
      comp_ref <- as.data.frame(raw_data[, brGDGT_cols]) / 100 + 1e-6
      clr_ref <- as.data.frame(compositions::clr(comp_ref))
      
      # Create reference_data list
      reference_data <- list(
        raw_values = raw_data[, brGDGT_cols],
        clr_values = clr_ref,
        npp = npp_col,
        coordinates = coords,
        sites = if("Site" %in% colnames(raw_data)) raw_data$Site else NULL
      )
      
    } else {
      # Fall back to RDS file
      ref_path <- system.file("extdata", "reference_data.rds", 
                              package = "brGDGTNPP")
      if(file.exists(ref_path)) {
        reference_data <- readRDS(ref_path)
      } else {
        # Create reference from model if available
        model <- load_npp_model()
        reference_data <- extract_reference_from_model(model)
      }
    }
  }
  
  # Define brGDGT columns
  brGDGT_cols <- c("fIa", "fIb", "fIc", "fIIa", "fIIa_", "fIIb", "fIIb_",
                   "fIIc", "fIIc_", "fIIIa", "fIIIa_", "fIIIb", "fIIIb_",
                   "fIIIc", "fIIIc_")
  
  # Extract and transform new data
  comp_new <- brgdt_data[, brGDGT_cols] / 100 + 1e-6
  clr_new <- as.data.frame(compositions::clr(comp_new))
  
  # Get reference CLR values
  clr_ref <- reference_data$clr_values
  
  # Ensure same columns
  common_cols <- intersect(colnames(clr_ref), colnames(clr_new))
  clr_ref <- clr_ref[, common_cols]
  clr_new <- clr_new[, common_cols]
  
  # Find nearest neighbors in brGDGT space
  nn_indices <- apply(clr_new, 1, function(x) {
    dists <- sqrt(rowSums(sweep(clr_ref, 2, x)^2))
    order(dists)[1:n_neighbors]
  })
  
  # If multiple samples, transpose
  if(nrow(clr_new) > 1) {
    nn_indices <- t(nn_indices)
  }
  
  # Get neighbor information
  neighbor_list <- list()
  for(i in 1:nrow(clr_new)) {
    neighbor_list[[i]] <- data.frame(
      sample_id = i,
      neighbor_rank = 1:n_neighbors,
      neighbor_index = nn_indices[i, ],
      distance = apply(clr_ref[nn_indices[i, ], ], 1, function(y) {
        sqrt(sum((y - clr_new[i, ])^2))
      })
    )
  }
  neighbors_df <- do.call(rbind, neighbor_list)
  
  # Geographic analysis if reference coordinates available
  geo_info <- NULL
  if(map && !is.null(reference_data$coordinates)) {
    
    geo_info <- list(
      ref_coords = reference_data$coordinates,
      neighbors = neighbors_df,
      n_neighbors = n_neighbors,
      new_sample_ids = 1:nrow(clr_new) 
    )
    
    # Add reference site names if available
    if(!is.null(reference_data$sites)) {
      geo_info$ref_sites <- reference_data$sites
    }
  }
  
  # Calculate Mahalanobis distances (using PCA if requested)
  if(use_pca) {
    # Combine for PCA
    all_data <- rbind(clr_ref, clr_new)
    pca_result <- prcomp(all_data, scale. = TRUE, center = TRUE)
    
    # Get PCs 
    var_explained <- cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2))
    n_pcs <- min(which(var_explained >= 0.95))
    
    pca_ref <- pca_result$x[1:nrow(clr_ref), 1:n_pcs]
    pca_new <- pca_result$x[(nrow(clr_ref)+1):nrow(all_data), 1:n_pcs]
    
    cov_matrix <- cov(pca_ref) + diag(1e-6, n_pcs)
    mean_vector <- colMeans(pca_ref)
    
    mahal_dists <- stats::mahalanobis(pca_new, mean_vector, cov_matrix)
    df_pcs <- n_pcs
    method_used <- "pca"
    
  } else {
    cov_matrix <- cov(clr_ref) + diag(1e-6, ncol(clr_ref))
    mean_vector <- colMeans(clr_ref)
    
    mahal_dists <- stats::mahalanobis(clr_new, mean_vector, cov_matrix)
    df_pcs <- ncol(clr_ref)
    method_used <- "regularized"
  }
  
  # Calculate p-values
  p_values <- 1 - stats::pchisq(mahal_dists, df = df_pcs)
  
  # Determine threshold
  threshold_value <- stats::qchisq(threshold, df = df_pcs)
  within_domain <- mahal_dists <= threshold_value
  
  # Calculate nearest neighbor distances
  nn_dists <- apply(clr_new, 1, function(x) {
    min(sqrt(rowSums(sweep(clr_ref, 2, x)^2)))
  })
  
  # Prepare results
  results <- list(
    within_domain = within_domain,
    mahalanobis_dist = mahal_dists,
    p_values = p_values,
    nearest_neighbor_dist = nn_dists,
    nearest_neighbors = neighbors_df,
    geographic = geo_info,
    threshold_value = threshold_value,
    method = method_used,
    n_pcs = ifelse(use_pca, n_pcs, ncol(clr_ref)),
    # Include raw reference data for plotting
    reference_raw = if(!is.null(reference_data$raw_values)) reference_data$raw_values else NULL,
    reference_npp = if(!is.null(reference_data$npp)) reference_data$npp else NULL
  )
  
  # Add sample IDs
  if(!is.null(brgdt_data[,1])) {
    results$sample_id <- brgdt_data[,1]
    names(results$within_domain) <- brgdt_data[,1]
    names(results$mahalanobis_dist) <- brgdt_data[,1]
  }
  
  # Generate plots
  if(plot) {
    if(use_pca) {
      results$plots <- generate_diagnostic_plots_pca(
        results, pca_ref, pca_new, brgdt_data, reference_data, geo_info
      )
    } else {
      results$plots <- generate_diagnostic_plots(
        results, clr_new, clr_ref, brgdt_data, geo_info
      )
    }
  }
  
  class(results) <- c("brGDGT_diagnostic", "list")
  return(results)
}

#' Generate enhanced diagnostic plots with PCA, ternary plot, and similarity map
#'
#' @param diagnostic_results Output from diagnose_applicability()
#' @param pca_ref PCA scores for reference data
#' @param pca_new PCA scores for new data
#' @param original_data Original brGDGT data (with raw percentages)
#' @param reference_data Reference dataset with raw values and NPP
#' @param geo_info Geographic information from diagnose_applicability()
#'
#' @return List of ggplot objects
#'
#' @keywords internal
#'
generate_diagnostic_plots_pca <- function(diagnostic_results, pca_ref, pca_new,
                                          original_data, reference_data, geo_info) {
  
  plots <- list()
  
  # ============================================
  # PLOT 1: Mahalanobis distances
  # ============================================
  df_mahal <- data.frame(
    Sample = 1:length(diagnostic_results$mahalanobis_dist),
    Distance = diagnostic_results$mahalanobis_dist,
    Within_Domain = diagnostic_results$within_domain
  )
  
  plots$mahalanobis <- ggplot2::ggplot(df_mahal, 
                                       ggplot2::aes(x = Sample, y = Distance, 
                                                    color = Within_Domain)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_hline(yintercept = diagnostic_results$threshold_value, 
                        linetype = "dashed", color = "red", alpha = 0.7) +
    ggplot2::labs(title = "Mahalanobis Distance",
                  subtitle = "Measures how unusual your samples are",
                  x = "Sample Index", y = "Mahalanobis Distance") +
    ggplot2::scale_color_manual(values = c("TRUE" = "#2E8B57", "FALSE" = "#DC143C")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  
  # ============================================
  # PLOT 2: PCA with training data
  # ============================================
  pca_all <- prcomp(rbind(pca_ref, pca_new))
  var_exp <- round(summary(pca_all)$importance[2,1:2] * 100, 1)
  
  pca_df <- data.frame(
    PC1 = c(pca_ref[,1], pca_new[,1]),
    PC2 = c(pca_ref[,2], pca_new[,2]),
    Type = c(rep("Global Training Data", nrow(pca_ref)), 
             rep("Your Samples", nrow(pca_new))),
    Within_Domain = c(rep(NA, nrow(pca_ref)), 
                      diagnostic_results$within_domain)
  )
  
  plots$pca <- ggplot2::ggplot() +
    ggplot2::geom_point(data = subset(pca_df, Type == "Global Training Data"),
                        ggplot2::aes(x = PC1, y = PC2), 
                        color = "gray70", alpha = 0.4, size = 1.5) +
    ggplot2::geom_point(data = subset(pca_df, Type == "Your Samples"),
                        ggplot2::aes(x = PC1, y = PC2, 
                                     color = Within_Domain,
                                     fill = Within_Domain), 
                        size = 4, shape = 21, stroke = 1.2) +
    ggplot2::scale_color_manual(values = c("TRUE" = "#2E8B57", "FALSE" = "#DC143C")) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#90EE90", "FALSE" = "#FF6B6B")) +
    ggplot2::labs(title = "PCA: Your Samples vs Global Training Data",
                  subtitle = paste0("PC1: ", var_exp[1], "%, PC2: ", var_exp[2], "%"),
                  x = "PC1", y = "PC2") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom",
                   plot.title = ggplot2::element_text(face = "bold", size = 14))
  
  # ============================================
  # PLOT 3: TERNARY PLOT - Methylation groups
  # ============================================
  
  # Define methylation groups
  tetramethylated <- c("fIa", "fIb", "fIc")
  pentamethylated <- c("fIIa", "fIIa_", "fIIb", "fIIb_", "fIIc", "fIIc_")
  hexamethylated  <- c("fIIIa", "fIIIa_", "fIIIb", "fIIIb_", "fIIIc", "fIIIc_")
  
  # ===== TRAINING DATA =====
  # Use raw values if available (from the Excel file)
  if(!is.null(diagnostic_results$reference_raw)) {
    message("...")
    
    training_raw <- diagnostic_results$reference_raw
    
    # Calculate methylation groups
    training_methyl <- data.frame(
      Tetra = rowSums(training_raw[, tetramethylated, drop = FALSE], na.rm = TRUE),
      Penta = rowSums(training_raw[, pentamethylated, drop = FALSE], na.rm = TRUE),
      Hexa = rowSums(training_raw[, hexamethylated, drop = FALSE], na.rm = TRUE)
    )
    
    # Convert from percentages to proportions if needed (assuming values are percentages 0-100)
    if(max(training_methyl, na.rm = TRUE) > 1) {
      training_methyl <- training_methyl / 100
    }
    
    # Add NPP
    if(!is.null(diagnostic_results$reference_npp)) {
      training_methyl$NPP <- diagnostic_results$reference_npp
    } else {
      training_methyl$NPP <- 500  # Default middle value
    }
    
  } else {
    # Fall back to CLR back-transformation
    message("Raw reference data not found, back-transforming from CLR")
    
    clr_train <- reference_data$clr_values
    comp_train <- exp(clr_train)
    comp_train <- comp_train / rowSums(comp_train)
    
    train_cols <- colnames(comp_train)
    training_methyl <- data.frame(
      Tetra = rowSums(comp_train[, intersect(tetramethylated, train_cols), drop = FALSE], na.rm = TRUE),
      Penta = rowSums(comp_train[, intersect(pentamethylated, train_cols), drop = FALSE], na.rm = TRUE),
      Hexa = rowSums(comp_train[, intersect(hexamethylated, train_cols), drop = FALSE], na.rm = TRUE)
    )
    
    if(!is.null(reference_data$npp)) {
      training_methyl$NPP <- reference_data$npp
    } else {
      training_methyl$NPP <- 500
    }
  }
  
  # Normalize training data to ensure sum to 1
  training_methyl$Total <- training_methyl$Tetra + training_methyl$Penta + training_methyl$Hexa
  training_methyl$Tetra_prop <- training_methyl$Tetra / training_methyl$Total
  training_methyl$Penta_prop <- training_methyl$Penta / training_methyl$Total
  training_methyl$Hexa_prop <- training_methyl$Hexa / training_methyl$Total
  
  # ===== TESTING DATA =====
  # Testing data comes as percentages (0-100)
  testing_methyl <- data.frame(
    Tetra = rowSums(original_data[, tetramethylated, drop = FALSE], na.rm = TRUE),
    Penta = rowSums(original_data[, pentamethylated, drop = FALSE], na.rm = TRUE),
    Hexa = rowSums(original_data[, hexamethylated, drop = FALSE], na.rm = TRUE)
  )
  
  # Convert from percentages to proportions if needed
  if(max(testing_methyl, na.rm = TRUE) > 1) {
    testing_methyl <- testing_methyl / 100
  }
  
  # Normalize testing data
  testing_methyl$Total <- testing_methyl$Tetra + testing_methyl$Penta + testing_methyl$Hexa
  testing_methyl$Tetra_prop <- testing_methyl$Tetra / testing_methyl$Total
  testing_methyl$Penta_prop <- testing_methyl$Penta / testing_methyl$Total
  testing_methyl$Hexa_prop <- testing_methyl$Hexa / testing_methyl$Total
  
  # Add within domain information
  testing_methyl$Within_Domain <- diagnostic_results$within_domain
  
  # Remove any rows with NA or infinite values
  train_valid <- complete.cases(training_methyl[, c("Tetra_prop", "Penta_prop", "Hexa_prop", "NPP")]) &
    is.finite(training_methyl$Tetra_prop) & 
    is.finite(training_methyl$Penta_prop) & 
    is.finite(training_methyl$Hexa_prop)
  
  test_valid <- complete.cases(testing_methyl[, c("Tetra_prop", "Penta_prop", "Hexa_prop")]) &
    is.finite(testing_methyl$Tetra_prop) & 
    is.finite(testing_methyl$Penta_prop) & 
    is.finite(testing_methyl$Hexa_prop)
  
  if(any(train_valid) && any(test_valid)) {
    
    # Training data for ternary plot
    tern_train <- data.frame(
      Tetra = training_methyl$Tetra_prop[train_valid],
      Penta = training_methyl$Penta_prop[train_valid],
      Hexa = training_methyl$Hexa_prop[train_valid],
      NPP = training_methyl$NPP[train_valid],
      Type = "Training"
    )
    
    # Testing data for ternary plot
    tern_test <- data.frame(
      Tetra = testing_methyl$Tetra_prop[test_valid],
      Penta = testing_methyl$Penta_prop[test_valid],
      Hexa = testing_methyl$Hexa_prop[test_valid],
      NPP = NA,
      Type = "Testing",
      Within_Domain = testing_methyl$Within_Domain[test_valid]
    )
    
    # ============================================
    # PLOT 3: TERNARY PLOT WITH PLOTLY
    # ============================================
    
    # Create ternary plot with plotly (NO CLASS CONFLICTS)
    if(requireNamespace("plotly", quietly = TRUE)) {
      
      # Prepare data for plotly
      tern_train$size_scaled <- 5 + (tern_train$NPP / max(tern_train$NPP)) * 15
      
      # Create ternary plot with plotly
      plots$ternary <- plotly::plot_ly(
        type = 'scatterternary',
        mode = 'markers'
      ) %>%
        # Training data
        plotly::add_trace(
          a = tern_train$Tetra,
          b = tern_train$Penta,
          c = tern_train$Hexa,
          marker = list(
            size = tern_train$size_scaled,
            color = 'gray',
            opacity = 0.5,
            line = list(width = 0)
          ),
          name = "Training data",
          text = ~paste('NPP:', round(tern_train$NPP)),
          hoverinfo = 'text'
        ) %>%
        # Testing data - within domain
        plotly::add_trace(
          a = tern_test$Tetra[tern_test$Within_Domain],
          b = tern_test$Penta[tern_test$Within_Domain],
          c = tern_test$Hexa[tern_test$Within_Domain],
          marker = list(
            size = 12,
            color = '#2E8B57',
            line = list(color = '#90EE90', width = 2)
          ),
          name = "Testing (within domain)"
        ) %>%
        # Testing data - outside domain
        plotly::add_trace(
          a = tern_test$Tetra[!tern_test$Within_Domain],
          b = tern_test$Penta[!tern_test$Within_Domain],
          c = tern_test$Hexa[!tern_test$Within_Domain],
          marker = list(
            size = 12,
            color = '#DC143C',
            line = list(color = '#FF6B6B', width = 2)
          ),
          name = "Testing (outside domain)"
        ) %>%
        plotly::layout(
          title = "brGDGT Methylation Groups",
          ternary = list(
            aaxis = list(title = "Tetra-methylated"),
            baxis = list(title = "Penta-methylated"),
            caxis = list(title = "Hexa-methylated")
          ),
          showlegend = TRUE
        )
      
    } else {
      warning("Package 'plotly' not installed. Skipping ternary plot. Install with: install.packages('plotly')")
      plots$ternary <- NULL
    }
    
  } else {
    warning("No valid data points for ternary plot after filtering")
    plots$ternary <- NULL
  }
  
  # ============================================
  # PLOT 4: Map showing TRAINING SITES most similar in brGDGT COMPOSITION
  # ============================================
  if(!is.null(geo_info) && !is.null(geo_info$ref_coords) && nrow(geo_info$ref_coords) > 0) {
    
    # Get neighbor information from diagnostic_results
    neighbors_df <- diagnostic_results$nearest_neighbors
    
    # Calculate similarity score (based on Euclidean distance in CLR)
    GDGT_similarity <- rep(0, nrow(geo_info$ref_coords))
    
    for(i in 1:nrow(geo_info$ref_coords)) {
      site_neighbors <- subset(neighbors_df, neighbor_index == i)
      if(nrow(site_neighbors) > 0) {
        GDGT_similarity[i] <- 1 / (1 + mean(site_neighbors$distance, na.rm = TRUE))
      }
    }
    
    # Normalize to 0-1 for visualization
    if(max(GDGT_similarity, na.rm = TRUE) > 0) {
      GDGT_similarity <- GDGT_similarity / max(GDGT_similarity, na.rm = TRUE)
    }
    
    # Prepare map data
    world <- ggplot2::map_data("world")
    
    ref_points <- data.frame(
      Longitude = geo_info$ref_coords$Longitude,
      Latitude = geo_info$ref_coords$Latitude,
      brGDGT_Similarity = GDGT_similarity,
      Is_Neighbor = 1:nrow(geo_info$ref_coords) %in% unique(neighbors_df$neighbor_index),
      stringsAsFactors = FALSE
    )
    
    if(!is.null(geo_info$ref_sites)) {
      ref_points$Site <- geo_info$ref_sites
    }
    
    # Find top similar sites
    valid_similarity <- is.finite(GDGT_similarity) & GDGT_similarity > 0
    if(any(valid_similarity)) {
      top_indices <- order(GDGT_similarity, decreasing = TRUE)[1:min(10, sum(valid_similarity))]
      top_sites <- ref_points[top_indices, ]
      top_sites <- top_sites[complete.cases(top_sites), ]
      
      if(nrow(top_sites) > 0) {
        plots$similarity_map <- ggplot2::ggplot() +
          ggplot2::geom_polygon(data = world, 
                                ggplot2::aes(x = long, y = lat, group = group),
                                fill = "gray90", color = "white", size = 0.2) +
          ggplot2::geom_point(data = ref_points,
                              ggplot2::aes(x = Longitude, y = Latitude),
                              color = "black", size = 1, alpha = 0.3) +
          ggplot2::geom_point(data = top_sites,
                              ggplot2::aes(x = Longitude, y = Latitude,
                                           size = brGDGT_Similarity),
                              color = "darkorange", alpha = 0.7) +
          ggplot2::scale_size_continuous(range = c(3, 8), guide = "none") +
          ggplot2::labs(title = "Most similar brGDGTs in training data",
                        subtitle = "Based on Euclidean distance in CLR-transformed space",
                        x = "Longitude", y = "Latitude") +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::coord_fixed(1.3)
      }
    }
  }
  
  return(plots)
}

#' Print method for diagnostic objects
#'
#' @param x brGDGT_diagnostic object
#' @param ... Additional arguments
#'
#' @export
print.brGDGT_diagnostic <- function(x, ...) {
  cat("brGDGT Applicability Domain Diagnostic\n")
  cat("======================================\n\n")
  cat("Samples within domain:", sum(x$within_domain), "of", length(x$within_domain), "\n")
  cat("Percentage:", round(mean(x$within_domain) * 100, 1), "%\n\n")
  
  if(!is.null(x$sample_id)) {
    df_summary <- data.frame(
      Sample = x$sample_id,
      Within_Domain = x$within_domain,
      Mahalanobis_Distance = round(x$mahalanobis_dist, 2),
      p_value = format.pval(x$p_values, digits = 3)
    )
    print(df_summary)
  }
}