#' Predict NPP from brGDGT data
#'
#' @description
#' Predicts Net Primary Productivity (NPP) from brGDGT data using
#' centered log-ratio (CLR) transformation and a pre-trained random forest model.
#'
#' @param brgdt_data A data frame containing brGDGT measurements. Must include:
#'   \itemize{
#'     \item First column: Sample ID or identifier
#'     \item Columns 2-16: brGDGT compounds (fIa, fIb, fIc, fIIa, fIIa_, fIIb,
#'           fIIb_, fIIc, fIIc_, fIIIa, fIIIa_, fIIIb, fIIIb_, fIIIc, fIIIc_)
#'   }
#' @param sample_type Character string specifying sample type. One of:
#'   "Lacustrine Sediment", "Lacustrine SPM", "Low DO Lacustrine SPM",
#'   "Peat", "Riverine Sediment and SPM", or "Soil"
#' @param include_ci Logical; if TRUE, returns 95% confidence intervals.
#'   Default = TRUE.
#'
#' @return A data frame with columns:
#'   \item{ID}{Sample identifier}
#'   \item{NPP_pred_mean}{Predicted NPP value}
#'   \item{NPP_CI95_lower}{Lower 95 percent confidence bound}
#'   \item{NPP_CI95_upper}{Upper 95 percent confidence bound}
#'   \item{NPP_range}{Width of confidence interval}
#'
#' @importFrom compositions clr
#' @importFrom randomForest randomForest
#' @importFrom stats cor quantile
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load your data
#' data(example_data)
#' 
#' # Predict NPP
#' results <- predict_npp(example_data, sample_type = "Soil")
#' head(results)
#' }
#'
#' @seealso
#' \code{\link{diagnose_applicability}} for checking data suitability
#'

predict_npp <- function(brgdt_data, sample_type, include_ci = TRUE) {
  
  # Load the pre-trained model
  model_path <- system.file("extdata", "model_brGDGT-NPP_all1.rds", 
                            package = "brGDGTNPP")
  model_all_data <- readRDS(model_path)
  
  # Define sample types
  sampletypes <- c("Lacustrine Sediment", "Lacustrine SPM", 
                   "Low DO Lacustrine SPM", "Peat", 
                   "Riverine Sediment and SPM", "Soil")
  
  # Extract brGDGT columns (assume first column is ID, next 15 are compounds)
  brGDGT_cols <- colnames(brgdt_data)[2:16]
  comp_new <- brgdt_data[, brGDGT_cols] / 100 + 1e-6
  
  # Apply CLR transformation
  clr_df <- as.data.frame(compositions::clr(comp_new))
  
  # Process sample type
  clr_df$Sampletype_fixed <- sample_type
  clr_df$Sampletype_fixed <- ifelse(
    clr_df$Sampletype_fixed %in% c("Lacustrine SPM", "Low DO Lacustrine SPM",
                                   "Riverine Sediment and SPM", "Lacustrine Meso/Microcosm"),
    "Aquatic_SPM_combined",
    clr_df$Sampletype_fixed
  )
  clr_df$Sampletype_fixed <- factor(clr_df$Sampletype_fixed)
  
  # Prepare model matrix
  clr_cols <- colnames(clr_df)[1:length(brGDGT_cols)]
  newdata_model <- clr_df[, clr_cols, drop = FALSE]
  
  # Create variables for sample types
  for(stype in sampletypes){
    colname <- paste0("Sampletype_fixed", stype)
    newdata_model[[colname]] <- ifelse(clr_df$Sampletype_fixed == stype, 1, 0)
  }
  
  # Ensure columns match model expectations
  model_cols <- c(clr_cols, paste0("Sampletype_fixed", sampletypes))
  model_cols <- model_cols[model_cols %in% colnames(newdata_model)]
  newdata_model <- newdata_model[, model_cols]
  
  # Make predictions
  if(include_ci) {
    # Get individual tree predictions for confidence intervals
    rf_preds_all <- predict(model_all_data$finalModel, 
                            newdata = newdata_model,
                            predict.all = TRUE) 
    
    # Calculate 95% confidence intervals
    ci_95_lower <- apply(rf_preds_all$individual, 1, 
                         stats::quantile, probs = 0.025, na.rm = TRUE)
    ci_95_upper <- apply(rf_preds_all$individual, 1, 
                         stats::quantile, probs = 0.975, na.rm = TRUE)
    
    # Create results data frame
    results <- data.frame(
      NPP_pred_mean = rf_preds_all$aggregate,
      NPP_CI95_lower = ci_95_lower,
      NPP_CI95_upper = ci_95_upper,
      NPP_range = ci_95_upper - ci_95_lower
    )
    
  } else {
    # Get only mean predictions
    rf_preds <- predict(model_all_data$finalModel, newdata = newdata_model)
    
    results <- data.frame(
      NPP_pred_mean = rf_preds
    )
  }
  
  # Add sample IDs if present
  if(!is.null(brgdt_data[,1])) {
    results <- cbind(ID = brgdt_data[,1], results)
  }
  
  return(results)
}