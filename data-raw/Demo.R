
# Install library
library(devtools)

# Modify directory
install.packages("D:/2026/brGDGT-NPP/brGDGTNPP_0.1.0.tar.gz", 
                 repos = NULL, 
                 type = "source")

# Load library
library(brGDGTNPP)


# Load example data
data(example_data)
head(example_data)

# Check if new brGDGT data falls within the trainig data space

?diagnose_applicability # Function details


check_data <- diagnose_applicability(
  brgdt_data = example_data,
  use_pca = TRUE,
  plot = TRUE,
  map = TRUE,
  n_neighbors = 50
)

# Gráfico 1: Distancias de Mahalanobis (rareza estadística)
check_data$plots$mahalanobis

# Gráfico 2: PCA (visión general de la composición química)
check_data$plots


check_data$plots$top_chemical_sites


par(mfrow = c(1, 3))


?predict_npp
?load_npp_model




results <- predict_npp(example_data, sample_type = "Soil")
head(results)



# ============================================================================
# COMPLETE DEMO: brGDGTNPP Package - All Functions Tutorial
# ============================================================================
# This script demonstrates how to use all functions in the brGDGTNPP package
# with the built-in example dataset.
# ============================================================================

# Load required packages
library(brGDGTNPP)
library(openxlsx)

# ============================================================================
# PART 1: LOAD EXAMPLE DATA
# ============================================================================

cat("\n=== PART 1: Loading Example Data ===\n")

# Load the example dataset included in the package
data(example_data)
head(example_data)

# Check column names (should be ID and brGDGT compounds)
names(example_data)

# ============================================================================
# PART 2: LOAD THE GLOBAL MODEL
# ============================================================================

# Load the pre-trained global model
global_model <- load_npp_model()
global_model

# ============================================================================
# PART 3: PREDICT NPP WITH GLOBAL MODEL
# ============================================================================


# Predict for Lacustrine Sediment samples
cat("\n--- Predicting for Lacustrine Sediment samples ---\n")
predictions_lac <- predict_npp(example_data, 
                               sample_type = "Lacustrine Sediment", 
                               include_ci = TRUE)

cat("\nFirst 5 predictions (Lacustrine Sediment):\n")
print(head(predictions_lac, 5))


# ============================================================================
# PART 4: DIAGNOSE APPLICABILITY DOMAIN
# ============================================================================


# Run diagnostic without plots first
diagnostic <- diagnose_applicability(example_data, 
                                     plot = TRUE, 
                                     use_pca = TRUE)

cat("\nDiagnostic Results:\n")
print(diagnostic)

# Check which samples are within domain
cat("\nSamples within domain:", sum(diagnostic$within_domain), 
    "out of", length(diagnostic$within_domain))
cat("\nPercentage:", round(mean(diagnostic$within_domain) * 100, 1), "%")

# View Mahalanobis distances
cat("\n\nMahalanobis distances (first 10):\n")
print(round(diagnostic$mahalanobis_dist[1:10], 3))

cat("\nThreshold value:", round(diagnostic$threshold_value, 3))


























library(brGDGTNPP)
data(example_data)














# 1. ¿El reference_data tiene coordenadas?
ref_path <- system.file("extdata", "reference_data.rds", package = "brGDGTNPP")
ref_data <- readRDS(ref_path)
names(ref_data)  # Debe incluir "coordinates"

# 2. ¿Las coordenadas tienen los nombres correctos?
head(ref_data$coordinates)

# 3. ¿Tus datos de ejemplo tienen coordenadas?
head(example_data[, c("Longitude", "Latitude")])
