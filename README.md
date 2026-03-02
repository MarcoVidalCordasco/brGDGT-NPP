# brGDGT-NPP
    R package to predict Net Primary Productivity (NPP) from the relative abundance of
    individual brGDGTs using centered log-ratio (CLR) transformation and a pre-trained 
    random forest model. Includes diagnostic tools to assess whether new 
    samples fall within the brGDGT composition space of the training data. When mean annual 
    temperature and precipitation are available, the package also computes 
    potential NPP using the Miami model, allowing comparison 
    between climate-based and brGDGT-based NPP estimates to identify 
    deviations from climatic expectations.

## Installation

# From GitHub
install.packages("devtools")

devtools::install_github("MarcoVidalCordasco/brGDGT-NPP", 
                         build_vignettes = TRUE, 
                         dependencies = TRUE)

# Load the package and check demo
head(example_data)

vignette("brGDGTNPP-demo")


