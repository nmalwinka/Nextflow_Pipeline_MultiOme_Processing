#!/storage/Software/packages/anaconda3/bin/Rscript

library(Signac)
library(Seurat)


Project <- "MBU_jv361_xxx_mgatk"

baseDir <- "/Users/mn367/Documents/MBU-Projects/MBU_groupLeader/MBU_projid_xxx"
setwd(baseDir)

MultiOme_obj <- readRDS("path/to/multiome.Rds")



message("+-------------------------------------------------------------------------------")
message("+        signac compute allele freq of MT variants  - need mgatk results        ")
message("+-------------------------------------------------------------------------------")

# load mgatk output
mito.data <- ReadMGATK(dir="Input/mgatk_input/final")


# create an assay
mito <- CreateAssayObject(counts = mito.data$counts)

# Subset to cell present in the scATAC-seq assat
mito <- subset(mito, cells = colnames(MultiOme_obj))

# add assay and metadata to the seurat object
MultiOme_obj[["mito"]] <- mito
MultiOme_obj <- AddMetaData(MultiOme_obj, metadata = mito.data$depth, col.name = "mtDNA_depth")


rm(mito)
gc()

# Dimension reduction and clustering
# Next we can run a standard dimension reduction and clustering workflow using the scATAC-seq data to identify cell clusters.
DefaultAssay(MultiOme_obj) <- "mito"
MultiOme_obj <- RunTFIDF(MultiOme_obj)
MultiOme_obj <- FindTopFeatures(MultiOme_obj, min.cutoff = 10)
MultiOme_obj <- RunSVD(MultiOme_obj)
MultiOme_obj <- RunUMAP(MultiOme_obj, reduction = "lsi", dims = 2:50)
MultiOme_obj <- FindNeighbors(MultiOme_obj, reduction = "lsi", dims = 2:50)
MultiOme_obj <- FindClusters(MultiOme_obj, resolution = 0.5, algorithm = 3)

#Find informative mtDNA variants
#Next, we can identify sites in the mitochondrial genome that vary across cells, and cluster the cells into clonotypes based on the frequency of these variants in the cells. Signac utilizes the principles established in the original mtscATAC-seq work of identifying high-quality variants.

variable.sites <- IdentifyVariants(MultiOme_obj, assay = "mito", refallele = mito.data$refallele)
VariantPlot(variants = variable.sites)


MultiOme_obj  <- AlleleFreq(MultiOme_obj, variants= rownames(high.conf), slot = "counts", assay = "mito", new.assay.name = "alleles")










