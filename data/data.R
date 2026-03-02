#' Example brGDGT dataset
#'
#' A dataset containing brGDGT measurements from 50 natural archive samples
#' to demonstrate package functionality.
#'
#' @format A data frame with 50 rows and 17 variables:
#' \describe{
#'   \item{ID}{Sample identifier}
#'   \item{fIa}{brGDGT Ia fractional abundance}
#'   \item{fIb}{brGDGT Ib fractional abundance}
#'   \item{fIc}{brGDGT Ic fractional abundance}
#'   \item{fIIa}{brGDGT IIa fractional abundance}
#'   \item{fIIa_}{brGDGT IIa' fractional abundance}
#'   \item{fIIb}{brGDGT IIb fractional abundance}
#'   \item{fIIb_}{brGDGT IIb' fractional abundance}
#'   \item{fIIc}{brGDGT IIc fractional abundance}
#'   \item{fIIc_}{brGDGT IIc' fractional abundance}
#'   \item{fIIIa}{brGDGT IIIa fractional abundance}
#'   \item{fIIIa_}{brGDGT IIIa' fractional abundance}
#'   \item{fIIIb}{brGDGT IIIb fractional abundance}
#'   \item{fIIIb_}{brGDGT IIIb' fractional abundance}
#'   \item{fIIIc}{brGDGT IIIc fractional abundance}
#'   \item{fIIIc_}{brGDGT IIIc' fractional abundance}
#'   \item{Age}{Sample age (years BP)}
#' }
#'
#' @source Natural archive samples from [your study area]
#'
#' @examples
#' data(example_data)
#' head(example_data)
"example_data"

#' Pre-trained random forest model for NPP prediction
#'
#' A random forest model trained on global brGDGT compilation to predict
#' Net Primary Productivity. The model uses CLR-transformed brGDGT abundances
#' and sample type as predictors.
#'
#' @format A randomForest model object with 1000 trees
#' @source Trained on global compilation of 2234 samples
#'
#' @examples
#' \dontrun{
#' model <- load_npp_model()
#' print(model)
#' }
"npp_model"