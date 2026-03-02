#' Load the pre-trained NPP prediction model
#'
#' @description
#' Loads the random forest model included with the package. The model was
#' trained on a global compilation of brGDGT data and is used by 
#' \code{\link{predict_npp}}.
#'
#' @return A random forest model object of class "train" from the caret package.
#'   The object contains the fitted model with 1000 trees.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load the global model
#' model <- load_npp_model()
#' 
#' # Examine model structure
#' print(model)
#' 
#' # Get model summary
#' summary(model)
#' }
#'
#' @seealso
#' \code{\link{predict_npp}} for making predictions,
#' \code{\link{diagnose_applicability}} for checking data suitability
#'
#' @importFrom compositions clr
#' @importFrom randomForest randomForest
#' @importFrom stats quantile

load_npp_model <- function() {
  model_path <- system.file("extdata", "model_brGDGT-NPP_all1.rds", 
                            package = "brGDGTNPP")
  
  if(!file.exists(model_path)) {
    stop("Model file not found. Please reinstall the package.")
  }
  
  readRDS(model_path)
}