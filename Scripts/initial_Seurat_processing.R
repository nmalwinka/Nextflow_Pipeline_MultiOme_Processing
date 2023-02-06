#!/storage/Software/packages/anaconda3/bin/Rscript

library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)


Project <- "MBU_jv361_xxx_ini"

baseDir <- "/Users/mn367/Documents/MBU-Projects/MBU_groupLeader/MBU_projid_xxx"
setwd(baseDir)


message("+-------------------------------------------------------------------------------")
message("+                     create MultiOme object                                    ")
message("+-------------------------------------------------------------------------------")


MultiOme_data <- Read10X(data.dir = paste0(baseDir,"Input/filtered_feature_bc_matrix"))

fragpath <- paste0(baseDir,"Input/atac_fragments.tsv.gz")


# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"


MultiOme_obj <- CreateSeuratObject(counts = MultiOme_data$`Gene Expression`, assay = "RNA", project = Project)

# create ATAC assay and add it to the object
MultiOme_obj[["ATAC"]] <- CreateChromatinAssay(counts = MultiOme_data$Peaks, sep = c(":", "-"), fragments = fragpath, annotation = annotation)



message("+-------------------------------------------------------------------------------")
message("+                            Quality control                                    ")
message("+-------------------------------------------------------------------------------")

DefaultAssay(MultiOme_obj) <- "ATAC"

MultiOme_obj <- NucleosomeSignal(MultiOme_obj)
MultiOme_obj <- TSSEnrichment(MultiOme_obj)


pdf(paste(Project, "_vlnPlt_", "ATACseq", ".pdf", sep = "_"), onefile=FALSE, width=10, height=6) # "poorQualRm_",
par(bg=NA)
VlnPlot( object = MultiOme_obj, features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), ncol = 4, pt.size = 0)
dev.off()



message("+-------------------------------------------------------------------------------")
message("+                   filter out low quality cells                                ")
message("+-------------------------------------------------------------------------------")


MultiOme_obj <- subset( x = MultiOme_obj, subset = nCount_ATAC < 100000 & nCount_RNA < 60000 &  nCount_ATAC > 1000 &  nCount_RNA > 1000 &  nucleosome_signal < 2 & TSS.enrichment > 1)


pdf(paste(Project, "_vlnPlt_", "ATACseq", "_subsetQC.pdf", sep = "_"), onefile=FALSE, width=10, height=6) # "poorQualRm_",
par(bg=NA)
VlnPlot( object = MultiOme_obj, features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), ncol = 4, pt.size = 0)
dev.off()




message("+-------------------------------------------------------------------------------")
message("+                             Peak calling                                      ")
message("+-------------------------------------------------------------------------------")

# The set of peaks identified using Cellranger often merges distinct peaks that are close together. This can create a problem for certain analyses, particularly motif enrichment analysis and peak-to-gene linkage. To identify a more accurate set of peaks, we can call peaks using MACS2 with the CallPeaks() function. Here we call peaks on all cells together, but we could identify peaks for each group of cells separately by setting the group.by parameter, and this can help identify peaks specific to rare cell populations.


# call peaks using MACS2
peaks <- CallPeaks(MultiOme_obj, macs2.path = "/Users/mn367/miniconda/envs/signac/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE) ### blacklist_hg38_unified
# blacklist for mouse mm10::: http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/


# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(MultiOme_obj),
  features = peaks,
  cells = colnames(MultiOme_obj)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
MultiOme_obj[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)



message("--------------------------------------------------------------------------------")
message("+                           cell cycle scoring                                  ")
message("+-------------------------------------------------------------------------------")

cc.genes.updated.2019
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

s.genes <- tolower(s.genes)
g2m.genes <- tolower(g2m.genes)

s.genes <- firstup(s.genes)
g2m.genes <- firstup(g2m.genes)


DefaultAssay(MultiOme_obj) <- "RNA"
MultiOme_obj <- CellCycleScoring(MultiOme_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

table(MultiOme_obj@meta.data[,c("Phase")])





message("+-------------------------------------------------------------------------------")
message("+                  Gene expression data processing                              ")
message("+-------------------------------------------------------------------------------")


# We can normalize the gene expression data using SCTransform, and reduce the dimensionality using PCA.

DefaultAssay(MultiOme_obj) <- "RNA"
MultiOme_obj <- SCTransform(MultiOme_obj, vars.to.regress =  c("S.Score", "G2M.Score"), return.only.var.genes = FALSE)
MultiOme_obj <- RunPCA(MultiOme_obj)





message("+-------------------------------------------------------------------------------")
message("+                 DNA accessibility data processing                             ")
message("+-------------------------------------------------------------------------------")

# Here we process the DNA accessibility assay the same way we would process a scATAC-seq dataset, by performing latent semantic indexing (LSI).

DefaultAssay(MultiOme_obj) <- "peaks"
MultiOme_obj <- FindTopFeatures(MultiOme_obj, min.cutoff = 5)
MultiOme_obj <- RunTFIDF(MultiOme_obj)
MultiOme_obj <- RunSVD(MultiOme_obj)





message("+-------------------------------------------------------------------------------")
message("+                              joint UMAP                                       ")
message("+-------------------------------------------------------------------------------")

# build a joint neighbor graph using both assays
MultiOme_obj <- FindMultiModalNeighbors(
  object = MultiOme_obj,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
MultiOme_obj <- RunUMAP(
  object = MultiOme_obj,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)











