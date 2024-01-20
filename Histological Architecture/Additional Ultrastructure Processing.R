## Initialization
library(dplyr)
library(patchwork)
library(DDRTree)
library(ggplot2)
library(stringr)
library(monocle3)
library(RANN)
library(viridis)
library(gridExtra)


## Loading in MATLAB-Processed Ultrastructure Data
# Load data (.csv should have one row per image and one column per matrix parameter, plus headers)
loaded_counts <- read.delim(file = "pdac.csv", header = TRUE, row.names=1, sep = ",")
loaded_counts[is.na(loaded_counts)] <- 0
loaded_counts <- t(loaded_counts)
matrix_object <- new_cell_data_set(loaded_counts)


## Processing of Patient Metadata
# Load in patient metadata
loaded_metadata <- read.delim(file = "Deidentified Patient Data.csv", header = TRUE, row.names=1, sep = ",")
loaded_metadata <- t(loaded_metadata)
# Populate image to patient key
loaded_key <- matrix(nrow = 2, ncol = ncol(loaded_counts))
rownames(loaded_key) <- c("Pic.ID","Patient.ID")
colnames(loaded_key) <- colnames(loaded_counts)
loaded_key["Pic.ID",] <- colnames(loaded_counts)
loaded_key["Patient.ID",] <- as.numeric(str_split(colnames(loaded_key),"-|_|\\.", simplify = TRUE)[,1])
# Set blank values to "NA" for compatibility
loaded_metadata[loaded_metadata == ""] <- NA
# Initialize patient metadata categories
colData(matrix_object)$patient_ID <- NA
colData(matrix_object)$cohort <- NA
colData(matrix_object)$age <- NA
colData(matrix_object)$sex <- NA
colData(matrix_object)$survivorship <- NA
colData(matrix_object)$overall_survival <- NA
colData(matrix_object)$disease_free_survival <- NA
colData(matrix_object)$t <- NA
colData(matrix_object)$n <- NA
colData(matrix_object)$ajcc_stage <- NA
colData(matrix_object)$grade <- NA
colData(matrix_object)$tumor_size <- NA
colData(matrix_object)$lymph_nodes <- NA
colData(matrix_object)$preop_ca199 <- NA
colData(matrix_object)$neoadjuvant_category <- NA
names(colData(matrix_object)$patient_ID) <- colnames(loaded_counts)
names(colData(matrix_object)$cohort) <- colnames(loaded_counts)
names(colData(matrix_object)$age) <- colnames(loaded_counts)
names(colData(matrix_object)$sex) <- colnames(loaded_counts)
names(colData(matrix_object)$survivorship) <- colnames(loaded_counts)
names(colData(matrix_object)$overall_survival) <- colnames(loaded_counts)
names(colData(matrix_object)$disease_free_survival) <- colnames(loaded_counts)
names(colData(matrix_object)$t) <- colnames(loaded_counts)
names(colData(matrix_object)$n) <- colnames(loaded_counts)
names(colData(matrix_object)$ajcc_stage) <- colnames(loaded_counts)
names(colData(matrix_object)$grade) <- colnames(loaded_counts)
names(colData(matrix_object)$tumor_size) <- colnames(loaded_counts)
names(colData(matrix_object)$lymph_nodes) <- colnames(loaded_counts)
names(colData(matrix_object)$preop_ca199) <- colnames(loaded_counts)
names(colData(matrix_object)$neoadjuvant_category) <- colnames(loaded_counts)
# Populate patient metadata based on patient ID lookup
colData(matrix_object)$patient_ID <- as.character(loaded_metadata['Patient.ID',as.character(loaded_key['Patient.ID',])])
colData(matrix_object)$cohort <- as.character(loaded_metadata['Cohort',as.character(loaded_key['Patient.ID',])])
colData(matrix_object)$age <- as.numeric(loaded_metadata['Age',as.character(loaded_key['Patient.ID',])])
colData(matrix_object)$sex <- as.character(loaded_metadata['Sex',as.character(loaded_key['Patient.ID',])])
colData(matrix_object)$survivorship <- as.character(loaded_metadata['Survivorship',as.character(loaded_key['Patient.ID',])])
colData(matrix_object)$overall_survival <- as.numeric(loaded_metadata['OS',as.character(loaded_key['Patient.ID',])])
colData(matrix_object)$disease_free_survival <- as.numeric(loaded_metadata['DFS',as.character(loaded_key['Patient.ID',])])
colData(matrix_object)$t <- as.character(loaded_metadata['T',as.character(loaded_key['Patient.ID',])])
colData(matrix_object)$n <- as.character(loaded_metadata['N',as.character(loaded_key['Patient.ID',])])
colData(matrix_object)$ajcc_stage <- as.character(loaded_metadata['AJCC.Stage',as.character(loaded_key['Patient.ID',])])
colData(matrix_object)$grade <- as.character(loaded_metadata['Grade',as.character(loaded_key['Patient.ID',])])
colData(matrix_object)$tumor_size <- as.numeric(loaded_metadata['Tumor.Size',as.character(loaded_key['Patient.ID',])])
colData(matrix_object)$lymph_nodes <- as.numeric(loaded_metadata['Lymph.Nodes',as.character(loaded_key['Patient.ID',])])
colData(matrix_object)$preop_ca199 <- as.numeric(loaded_metadata['Preop.CA199',as.character(loaded_key['Patient.ID',])])
colData(matrix_object)$neoadjuvant_category <- as.character(loaded_metadata['Neoadjuvant.Category',as.character(loaded_key['Patient.ID',])])
# Store ultrastructural parameter names 
rowData(matrix_object)$gene_name <- rownames(matrix_object)
rowData(matrix_object)$gene_short_name <- rowData(matrix_object)$gene_name


## Trajectory Analysis (DDRTree via Monocle Wrappers)
# PCA
matrix_object <- preprocess_cds(cds = matrix_object, method = "PCA", num_dim=10)
# Run UMAP to generate dummy slots for population 
matrix_object <- reduce_dimension(cds = matrix_object, reduction_method = "UMAP", preprocess_method="PCA", build_nn_index = TRUE) 
# Populate UMAP coordinates with pre-calculated values from MATLAB
loadedUMAP <- t(read.csv("umap.csv", header = FALSE))
reducedDim(matrix_object, "UMAP")[,1] <- loadedUMAP[1,]
reducedDim(matrix_object, "UMAP")[,2] <- loadedUMAP[2,]
matrix_object <- cluster_cells(cds = matrix_object, reduction_method = "UMAP", resolution = 5e-5) 
# Build trajectories (DDRTree via Monocle)
matrix_object <- cluster_cells(cds = matrix_object, reduction_method = "UMAP", resolution = 1e-8) 
matrix_object <- learn_graph(matrix_object, use_partition = TRUE, learn_graph_control = list(minimal_branch_len = 5))
# Order cells in pseudotime 
matrix_object <- order_cells(matrix_object, reduction_method = "UMAP")
# Store pseudotime for each datapoint
pseudotime <- matrix_object@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]


## Mapping of New Data
# Perform entire analysis pipeline up to, but not including PCA (preprocess_cds)
# for new set of images and store object as test_object. 
# All downstream processing of new dataset will be performed by model transformation.
# Set # of PC's (default: # used in original model)
n_pcs <- ncol(reducedDim(matrix_object,"PCA"))
# Make object for NN fitting
nn_object <- matrix_object
# Rebuild NN index
nn_object <- reduce_dimension(cds = nn_object, reduction_method = "UMAP", preprocess_method="PCA", build_nn_index = TRUE) 
# Save original pseudotime values
nn_object@principal_graph_aux@listData[["UMAP"]][["pseudotime"]] <- matrix_object@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
# Import original model into query dataset
save_transform_models(nn_object, "reference_model")
test_object <- load_transform_models(test_object, "reference_model")
# Apply original transform model to query dataset
test_object <- preprocess_transform(test_object)
test_object <- reduce_dimension_transform(test_object)
# Map pseudotime values to new images using nearest neighbors
nn_indices <- as.vector(nn2(data = reducedDim(nn_object, "PCA"), query = reducedDim(test_object, "PCA"), k = 1)$nn.idx)
mapped_pseudotime <- nn_object@principal_graph_aux@listData[["UMAP"]][["pseudotime"]][as.numeric(nn_indices)]
test_object$mapped_pseudotime <- mapped_pseudotime
# Map UMAP coordinates to new images using nearest neighbors
mapped_umap <- t(as.matrix(reducedDim(matrix_object, "UMAP")))[,as.numeric(nn_indices)]
rownames(mapped_umap) <- c("UMAP.X","UMAP.Y")

