# Load necessary libraries
library(ggcorrplot)
library(factoextra)
library(dplyr)
library(tidyr)
library(factoextra)
library(RColorBrewer)
library(ggplot2)
library(ellipse)


# Combining Datasets for RStudio PCA ------------------------------------

# Files to Download for Testing:
file_paths <- c(
  "/home/vikentiy/RStudio Data/HWC/HWC1.1",
  "/home/vikentiy/RStudio Data/HWC/HWC1.2",
  "/home/vikentiy/RStudio Data/HWC/HWC1.3",
  "/home/vikentiy/RStudio Data/HWC/HWC2.1",
  "/home/vikentiy/RStudio Data/HWC/HWC2.2",
  "/home/vikentiy/RStudio Data/HWC/HWC2.3",
  "/home/vikentiy/RStudio Data/HWC/HWC3.1",
  "/home/vikentiy/RStudio Data/HWC/HWC3.2",
  "/home/vikentiy/RStudio Data/HWC/HWC3.3",
  "/home/vikentiy/RStudio Data/HWC/HWC4.1",
  "/home/vikentiy/RStudio Data/HWC/HWC4.2",
  "/home/vikentiy/RStudio Data/HWC/HWC4.3",
  "/home/vikentiy/RStudio Data/HWC/HWC5.1",
  "/home/vikentiy/RStudio Data/HWC/HWC5.2",
  "/home/vikentiy/RStudio Data/HWC/HWC5.3",
  "/home/vikentiy/RStudio Data/H13/H131.1",
  "/home/vikentiy/RStudio Data/H13/H131.2",
  "/home/vikentiy/RStudio Data/H13/H131.3",
  "/home/vikentiy/RStudio Data/H13/H132.1",
  "/home/vikentiy/RStudio Data/H13/H132.2",
  "/home/vikentiy/RStudio Data/H13/H132.3",
  "/home/vikentiy/RStudio Data/H13/H133.1",
  "/home/vikentiy/RStudio Data/H13/H133.2",
  "/home/vikentiy/RStudio Data/H13/H133.3",
  "/home/vikentiy/RStudio Data/H13/H134.1",
  "/home/vikentiy/RStudio Data/H13/H134.2",
  "/home/vikentiy/RStudio Data/H13/H134.3",
  "/home/vikentiy/RStudio Data/H13/H135.1",
  "/home/vikentiy/RStudio Data/H13/H135.2",
  "/home/vikentiy/RStudio Data/H13/H135.3",
  "/home/vikentiy/RStudio Data/Crab Biochar (CSB)/CSB1.1",
  "/home/vikentiy/RStudio Data/Crab Biochar (CSB)/CSB1.2",
  "/home/vikentiy/RStudio Data/Crab Biochar (CSB)/CSB1.3",
  "/home/vikentiy/RStudio Data/Crab Biochar (CSB)/CSB2.1",
  "/home/vikentiy/RStudio Data/Crab Biochar (CSB)/CSB2.2",
  "/home/vikentiy/RStudio Data/Crab Biochar (CSB)/CSB2.3",
  "/home/vikentiy/RStudio Data/Crab Biochar (CSB)/CSB3.1",
  "/home/vikentiy/RStudio Data/Crab Biochar (CSB)/CSB3.2",
  "/home/vikentiy/RStudio Data/Crab Biochar (CSB)/CSB3.3",
  "/home/vikentiy/RStudio Data/Crab Biochar (CSB)/CSB4.1",
  "/home/vikentiy/RStudio Data/Crab Biochar (CSB)/CSB4.2",
  "/home/vikentiy/RStudio Data/Crab Biochar (CSB)/CSB4.3",
  "/home/vikentiy/RStudio Data/Crab Biochar (CSB)/CSB5.1",
  "/home/vikentiy/RStudio Data/Crab Biochar (CSB)/CSB5.2",
  "/home/vikentiy/RStudio Data/Crab Biochar (CSB)/CSB5.3"
)

# Read all datasets
all_datasets <- lapply(file_paths, read.csv, header = FALSE)

# Extract column names from file paths
column_names <- sapply(strsplit(file_paths, "/"), tail, 1)

# Create a function to select relevant columns for each dataset
select_columns <- function(dataset, index) {
  if (index == 1) {
    # Keep both columns for subsequent datasets
    return(dataset[, c(1, 2)])
  } else {
    # Keep only the second column for subsequent datasets
    return(dataset[, 2])
  }
}

# Apply the function to each dataset and assign column names
all_datasets <- mapply(select_columns, all_datasets, seq_along(all_datasets), SIMPLIFY = FALSE)

# Rename columns based on file names
names(all_datasets) <- column_names

# Extract sample names without numbers
sample_names <- substr(column_names, 1, 3)

# Add a row with sample names to the merged dataset
merged_dataset <- do.call(cbind, all_datasets)
# Get the wavenumber of the original dataset
wavenumber <- merged_dataset[,1]
merged_dataset <- merged_dataset[,2:ncol(merged_dataset)]

# Transpose (swap rows and columns)
merged_dataset <- t(merged_dataset)

merged_dataset <- data.frame(sample_names = sample_names, merged_dataset)


# Create a data frame with sample names as the first row
PCAcombined <- data.frame(merged_dataset, stringsAsFactors = FALSE)

# Extract the first row as the factor variable
sample_names_factor <- as.factor(PCAcombined[, 1])

# Remove the first row (factor variable) from the data
PCAcombined <- PCAcombined[,-1]

# Clean and preprocess the data
PCAcombined <- PCAcombined[complete.cases(PCAcombined), ]  # Remove missing values
PCAcombined <- as.data.frame(sapply(PCAcombined, as.numeric))  # Ensure all columns are numeric
PCAcombined <- PCAcombined[is.finite(rowSums(PCAcombined)), ]  # Remove infinite values
PCAcombined_normalized <- scale(PCAcombined)  # Standardize the data

# Compute and visualize the correlation matrix
PCAcombined_corr_matrix <- cor(PCAcombined_normalized)

# Apply PCA
PCAcombined.pca <- princomp(PCAcombined_corr_matrix)

#Scree Plot
fviz_eig(PCAcombined.pca, addlabels = TRUE)

# Extract loading vectors for PC1, PC2, and PC3
loadings_pc1 <- PCAcombined.pca$loadings[, 1]
loadings_pc2 <- PCAcombined.pca$loadings[, 2]

# PCAcombined is the data where the rows are the individual samples
pc1_scores <- as.matrix(PCAcombined_normalized) %*% loadings_pc1
pc2_scores <- as.matrix(PCAcombined_normalized) %*% loadings_pc2

# Get unique sample names and assign unique colors
unique_sample_names <- unique(sample_names)
num_unique_samples <- length(unique_sample_names)
sample_colors <- rainbow(num_unique_samples)

# Create a named vector to map sample names to colors
sample_colors_map <- sample_colors[match(sample_names, unique_sample_names)]

scaled_scores_pc1 <- scale(pc1_scores)
scaled_scores_pc2 <- scale(pc2_scores)

# Creating a scatter plot with unique colors for each sample
plot(scaled_scores_pc1, scaled_scores_pc2, 
     xlab = "Principal Component 1 Scores",
     ylab = "Principal Component 2 Scores",
     main = "PCA Projection of Samples",
     col = sample_colors_map)

# Adding labels for each point using sample names
text(pc1_scores, pc2_scores, labels = sample_names, pos = 3, cex = 0.6) # cex controls text size

# Define the file path and name for the CSV file
csv_file <- "/home/vikentiy/Thesis Codes/Fixed codes/Plots/compare_load.csv"
csv1_file <- "/home/vikentiy/Thesis Codes/Fixed codes/Plots/compare_score.csv"

data_matrix <- data.frame(wavenumber, loadings_pc1, loadings_pc2)

# Write the columns of data to CSV files for python visualization
write.csv(data_matrix, file = csv_file, row.names = FALSE, quote = FALSE)
write.csv(cbind(scaled_scores_pc1, scaled_scores_pc2, sample_names), file = csv1_file, row.names = FALSE, quote = FALSE)
     
