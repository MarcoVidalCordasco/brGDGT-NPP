#' Predict NPP from mean annual temperature and precipitation with Miami model (Lieth, 1975)
#'
#' @description
#' Predicts potential Net Primary Productivity (NPP) from temperature and precipitation data using
#' the Miami model (Lieth, H. Modeling the Primary Productivity of the World. 237–263 (1975) 
#' doi:10.1007/978-3-642-80913-2_12.
#'
#' @param climate_data A data frame containing temperatures and precipitations. Must include:
#'   \itemize{
#'     \item First column: Sample ID or identifier (e.g., sample depth, archaeological level)
#'     \item Columns 2-3: mean annual temperature (MAT) in °C
#'           and mean annual precipitation (MAP) in mm
#'   }
#' @param npp_observed Optional numeric vector or data frame with observed/predicted NPP values 
#'   from brGDGTs. Required if hanpp = TRUE.
#' @param hanpp Logical; if TRUE, returns Human Appropriation of NPP (HANPP), computed
#'   as the difference between climate-derived NPP (Miami model) and brGDGT-derived NPP.
#'   Default is FALSE.
#'
#' @return A data frame with columns:
#'   \item{ID}{Sample identifier}
#'   \item{NPP_pred_mean}{Predicted potential NPP value from Miami model (g/m²/yr)}
#'   \item{HANPP}{Human Appropriation of NPP (g/m²/yr) - only if hanpp = TRUE}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example climate data
#' climate_data <- data.frame(
#'   ID = 1:5,
#'   MAT = c(10, 15, 20, 25, 30),
#'   MAP = c(500, 1000, 1500, 2000, 2500)
#' )
#' 
#' # Predict potential NPP with Miami model
#' results <- miami_model(climate_data)
#' head(results)
#' 
#' # Example with brGDGT-derived NPP to calculate HANPP
#' brgdt_npp <- c(550, 680, 720, 750, 780)  # example values
#' results_hanpp <- miami_model(climate_data, 
#'                               npp_observed = brgdt_npp, 
#'                               hanpp = TRUE)
#' head(results_hanpp)
#' }
#'
#' @seealso
#' \code{\link{predict_npp}} for computing brGDGT-derived NPP

miami_model <- function(climate_data, npp_observed = NULL, hanpp = FALSE) {
  
  # Check input
  if(ncol(climate_data) < 3) {
    stop("climate_data must have at least 3 columns: ID, MAT, MAP")
  }
  
  # Extract variables
  ID <- climate_data[, 1]
  MAT <- climate_data[, 2]
  MAP <- climate_data[, 3]
  
  # Compute NPP from temperature (g/m²/yr)
  NPP_MAT <- 3000 / (1 + exp(1.315 - 0.119 * MAT))
  
  # Compute NPP from precipitation (g/m²/yr)
  NPP_MAP <- 3000 * (1 - exp(-0.000664 * MAP))
  
  # Potential NPP
  NPP_pot <- pmin(NPP_MAT, NPP_MAP)
  
  # Round potential NPP
  NPP_pot <- round(NPP_pot, 2)
  
  # Create base results data frame
  results <- data.frame(
    ID = ID,
    NPP_pred_mean = NPP_pot
  )
  
  # Calculate HANPP if requested
  if(hanpp) {
    if(is.null(npp_observed)) {
      stop("npp_observed must be provided when hanpp = TRUE")
    }
    
    # Check lengths match
    if(length(npp_observed) != length(NPP_pot)) {
      stop("Length of npp_observed must match number of rows in climate_data")
    }
    
    # Calculate HANPP (positive = appropriation, negative = enhancement)
    HANPP <- NPP_pot - npp_observed
    
    # Add to results
    results$HANPP <- round(HANPP, 2)
  }
  
  return(results)
}