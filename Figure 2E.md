# HBV2024_flow-cytometry-data_Figure 2E

library(Matrix)
library(readr)
library(data.table)
library(dplyr)
library(plyr)
library(HDCytoData)
library(flowCore)
library(CATALYST)
library(stringr)
library(ggplot2)

# Specify directories
file_paths <- list(paste0(getwd(),"/Fig2E/"))

# Name file_paths
names(file_paths) <- file_paths

# Search for all files in the specified directories and extract files by a given extension
files_list            <- lapply(file_paths, list.files, recursive=T)
files_by_ext          <- lapply(files_list, function(x){x[endsWith(x, suffix=".fcs")]} )

# Get complete paths to all files
all_file_paths        <- unlist(lapply(seq_along(files_by_ext), function(x) {  paste(names(files_by_ext[x]), files_by_ext[[x]], sep="") } ))
names(all_file_paths) <- lapply(strsplit(all_file_paths,split="/"), function(x) { sub(".fcs","",x[length(x)]) } )
file_names <- unname(unlist(lapply(strsplit(unlist(files_by_ext),split = "/"),tail,1)))

# Read the FCS files in the directory fcs/
data_list <- list()
fcs.par <- list()
for (i in c(1:length(file_names))) {
  data_list[[i]] <- read.FCS(as.character(all_file_paths[i]), alter.names = T)
  fcs.par[[i]] <- as.character(data_list[[i]]@parameters@data$name)
}

# Find all the measured FACS parameters
common.par <- Reduce(intersect, fcs.par)
common.par

# Get all the fluorescent markers measured during FACS analysis
common.par <- common.par[7:21]

# Subset the FCS files for fluorescent markers measured during FACS analysis
dir.create(paste0(getwd(),"/subsetFig2E/"))
sub_data_list <- list()
subset_file_names <- list()

for (i in c(1:length(file_names))) {
  sub_data_list[[i]] <- data_list[[i]][,common.par]
  write.FCS(sub_data_list[[i]] , paste(getwd(),"/subsetFig2E/",paste(sub(".fcs","",file_names[[i]]),"_subset",".fcs", sep = ""), sep = ""))
  subset_file_names[[i]] <- paste(sub(".fcs","",file_names[[i]]),"_subset",".fcs", sep = "")
}

# Read the subset FCS files
fs <- read.flowSet(path =paste0(getwd(),"/subsetFig2E/"))

# Create patient and sample annotations: 
md <- as.data.frame(as.character(unlist(subset_file_names)))
colnames(md) <- c("file_name")
md$patient_id <- str_split_fixed(md$file_name, "_", 5)[,1]
md$sample_id <-  "core"
md$sample_id[str_detect(md$file_name, "pol")] <- "pol"
md$condition <- paste(md$sample_id, md$patient_id, sep="_")

# Create panel annotations
panel <- data.frame(as.character(sub_data_list[[1]]@parameters@data$name),as.character(sub_data_list[[1]]@parameters@data$desc))
rownames(panel) <- c(1:dim(panel)[1])
panel$marker_class <- c("type")
colnames(panel) <- c("fcs_colname",	"antigen", "marker_class")

# Exclude the parameters not to be included in the dimensionality reduction
panel[grep("APC.Cy7.A|Ax700.A|BUV395.A|BUV496.A|PE.A",panel$fcs_colname),]$marker_class <- c("state")

# Create single cell experiment
sce <- prepData(
  fs,
  panel = panel,
  md = md,
  features = panel$fcs_colname,
  transform = TRUE,
  cofactor = 150,
  by_time = F, FACS = T
)

# Number of cells in each cohort
for (i in 1:length(unique(sce$sample_id))) {
  print(paste("The number of cells in",as.character(unique(sce$sample_id)[i]), "are",length(sce$sample_id[sce$sample_id %in% as.character(unique(sce$sample_id)[i])])))
}

# Plot the number of cells in each condition
pdf("HBV CHRONIC NUC VS NAIVEplot1_number of cells in each condition.pdf")
plotCounts(sce, group_by = "condition", color_by = "sample_id")
graphics.off()

# Run t-SNE/UMAP on at most 1000 cells per sample
set.seed(1234)
sce <- runDR(sce, "TSNE", features = "type", cells = 1000)
pdf("PLot_TSNE.pdf")
plotDR(sce, color_by = "sample_id", "TSNE") + geom_point(size=1)
graphics.off()

# Plot scaled expression of the desired markers
cdx <- rownames(sce)[c(1,6:11,13:15)]
pdf("Plot_TSNE scaled expression of markers.pdf")
plotDR(sce, scale = T,color_by = cdx, ncol = 4, a_pal = rev(hcl.colors(10, "Spectral"))) + geom_point(size=0.5) 
graphics.off()

