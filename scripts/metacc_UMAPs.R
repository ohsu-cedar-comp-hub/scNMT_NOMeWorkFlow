args <- commandArgs()

help <- function(){
    cat("metacc_UMAPs.R :
Make UMAPs using accessibility and methylation information.\n")
    cat("Usage: \n")
    cat("--outDir     : output directory            [required]\n")
    cat("--accDir     : path to accessibility files [required]\n")
    cat("--metDir     : path to methylation files   [required]\n")
    cat("--anno       : path to annotation file     [required]\n")
    cat("--dataname   : Sample Name                 [required]\n")
    cat("\n")
    q()
}

io <- list()
opts <- list()

if( !is.na(charmatch("--help",args)) || !is.na(charmatch("--help",args)) ){
    help()
} else {
    io$outDir    <- sub( '--outDir=', '', args[grep('--outDir=', args)] )
    io$met_dir   <- sub( '--metDir=', '', args[grep('--metDir=', args)] )
    io$acc_dir   <- sub( '--accDir=', '', args[grep('--accDir=', args)] )
    io$anno      <- sub( '--anno=', '', args[grep('--anno=', args)] )
    io$data_name <- sub( '--dataname=', '', args[grep('--dataname=', args)] )
}

#io$outDir <- "plots/UMAPs"
#io$met_dir <- "data/met"
#io$acc_dir <- "data/acc"
#io$anno <- "data/anno"
#io$data_name <- "p0_params"

library(remotes)

install_version("Matrix", version = "1.6-5", repos = "http://cran.us.r-project.org")
install.packages("irlba", type = "source", repos = "http://cran.us.r-project.org")

library(Seurat)
library(dplyr)
library(ggplot2)
library(Signac)
library(reshape2)
library(tibble)

data_name <- io$data_name

acc_body <- read.table(paste0(io$acc_dir,"/body.tsv.gz"), sep = "\t", header = T)
met_body <- read.table(paste0(io$met_dir,"/body.tsv.gz"), sep = "\t", header = T)
body_anno <- read.table(paste0(io$anno, "/body1000.bed"), sep = "\t")
colnames(body_anno) <- c("chr", "start", "end", "strand", "id", "anno")

acc_body <- merge(acc_body, body_anno, by = c("id", "anno"))
met_body <- merge(met_body, body_anno, by = c("id", "anno"))

acc_body$count <- "unassigned"
acc_body$count[acc_body$rate > 50] <- 1
acc_body$count[acc_body$rate <= 50] <- 0
acc_body$count <- as.numeric(acc_body$count)

met_body$count <- "unassigned"
met_body$count[met_body$rate > 50] <- 1
met_body$count[met_body$rate <= 50] <- 0
met_body$count <- as.numeric(met_body$count)

acc_body$peak <- paste(acc_body$chr, acc_body$start, acc_body$end, sep = "-")
met_body$peak <- paste(met_body$chr, met_body$start, met_body$end, sep = "-")

acc_body_counts <- dcast(data=acc_body, formula=peak ~ sample, value.var="count") %>% remove_rownames %>% column_to_rownames(var="peak")

met_body_counts <- dcast(data=met_body, formula=peak ~ sample, value.var="count") %>% remove_rownames %>% column_to_rownames(var="peak")


acc_promoter <- read.table(paste0(io$acc_dir, "/promoter.tsv.gz"), sep = "\t", header = T)
met_promoter <- read.table(paste0(io$met_dir, "/promoter.tsv.gz"), sep = "\t", header = T)
promoter_anno <- read.table(paste0(io$anno, "/promoters1000.bed"), sep = "\t")
colnames(promoter_anno) <- c("chr", "start", "end", "strand", "id", "anno")

acc_promoter <- merge(acc_promoter, promoter_anno, by = c("id", "anno"))
met_promoter <- merge(met_promoter, promoter_anno, by = c("id", "anno"))

acc_promoter$count <- "unassigned"
acc_promoter$count[acc_promoter$rate > 50] <- 1
acc_promoter$count[acc_promoter$rate <= 50] <- 0

met_promoter$count <- "unassigned"
met_promoter$count[met_promoter$rate > 50] <- 1
met_promoter$count[met_promoter$rate <= 50] <- 0

acc_promoter$peak <- paste(acc_promoter$chr, acc_promoter$start, acc_promoter$end, sep = "-")
met_promoter$peak <- paste(met_promoter$chr, met_promoter$start, met_promoter$end, sep = "-")

acc_promoter_counts <- dcast(data=acc_promoter, formula=peak ~ sample, value.var="count") %>% remove_rownames %>% column_to_rownames(var="peak")

met_promoter_counts <- dcast(data=met_promoter, formula=peak ~ sample, value.var="count") %>% remove_rownames %>% column_to_rownames(var="peak")

acc_calledPeaks <- read.table(paste0(io$acc_dir, "/", data_name, "_called_peaks.tsv.gz"), sep = "\t", header = T)
met_calledPeaks <- read.table(paste0(io$met_dir, "/", data_name, "_called_peaks.tsv.gz"), sep = "\t", header = T)
calledPeaks_anno <- read.table(paste0(io$anno, "/", data_name, "_called_peaks.bed"), sep = "\t")
colnames(calledPeaks_anno) <- c("chr", "start", "end", "strand", "id", "anno")

acc_calledPeaks <- merge(acc_calledPeaks, calledPeaks_anno, by = c("id", "anno"))
met_calledPeaks <- merge(met_calledPeaks, calledPeaks_anno, by = c("id", "anno"))

acc_calledPeaks$count <- "unassigned"
acc_calledPeaks$count[acc_calledPeaks$rate > 50] <- 1
acc_calledPeaks$count[acc_calledPeaks$rate <= 50] <- 0
acc_calledPeaks$count <- as.numeric(acc_calledPeaks$count)

met_calledPeaks$count <- "unassigned"
met_calledPeaks$count[met_calledPeaks$rate > 50] <- 1
met_calledPeaks$count[met_calledPeaks$rate <= 50] <- 0
met_calledPeaks$count <- as.numeric(met_calledPeaks$count)

acc_calledPeaks$peak <- paste(acc_calledPeaks$chr, acc_calledPeaks$start, acc_calledPeaks$end, sep = "-")
met_calledPeaks$peak <- paste(met_calledPeaks$chr, met_calledPeaks$start, met_calledPeaks$end, sep = "-")

acc_calledPeaks_counts <- dcast(data=acc_calledPeaks, formula=peak ~ sample, value.var="count") %>% remove_rownames %>% column_to_rownames(var="peak")

met_calledPeaks_counts <- dcast(data=met_calledPeaks, formula=peak ~ sample, value.var="count") %>% remove_rownames %>% column_to_rownames(var="peak")

acc_body_counts[which(is.na(acc_body_counts), arr.ind = T)] <- 0

acc_body_chrom_assay <- CreateChromatinAssay(
  counts = acc_body_counts,
  sep = c("-", "-"),
  min.cells = 5,
  min.features = 50
)

acc_body_SO <- CreateSeuratObject(
  counts = acc_body_chrom_assay,
  assay = "peaks")

acc_body_SO$Plate <- sub("_.*", "", colnames(acc_body_SO))
acc_body_SO$Sample <- sub("-.*", "", acc_body_SO$Plate)

#install.packages("Matrix", type = "source")
#install_version("Matrix", version = "1.6-5", repos = "http://cran.us.r-project.org")
#install.packages("irlba", type = "source")

acc_body_SO <- RunTFIDF(acc_body_SO)
acc_body_SO <- FindTopFeatures(acc_body_SO, min.cutoff = 'q0')
acc_body_SO <- RunSVD(acc_body_SO)

acc_body_SO <- RunUMAP(object = acc_body_SO, reduction = 'lsi', dims = 2:20)
acc_body_SO <- FindNeighbors(object = acc_body_SO, reduction = 'lsi', dims = 2:20)
#acc_body_SO <- FindClusters(object = acc_body_SO, verbose = FALSE, algorithm = 3)

pdf(paste0(io$outDir, "/", data_name, "_accBodyUMAP.pdf"))
DimPlot(object = acc_body_SO, group.by = "Plate")
dev.off()

#Idents(acc_body_SO) <- "Sample"

#patient.markers <- FindAllMarkers(acc_body_SO, only.pos = T, test.use = 'LR', latent.vars = 'nCount_peaks')

#top.markers <- patient.markers %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)

#acc_body_SO <- ScaleData(acc_body_SO)

#pdf("allPatients_accDAheatmap.pdf")
#DoHeatmap(acc_body_SO, features = top.markers$gene) + theme(text = element_text(size = 7))
#dev.off()

#write.table(patient.markers, file = "allPatients_accDAlist.tsv", sep = "\t", row.names = F, quote = F)

met_body_counts[which(is.na(met_body_counts), arr.ind = T)] <- 0

met_body_chrom_assay <- CreateChromatinAssay(
  counts = met_body_counts,
  sep = c("-", "-"),
  min.cells = 5,
  min.features = 50
)

met_body_SO <- CreateSeuratObject(
  counts = met_body_chrom_assay,
  assay = "peaks")

met_body_SO$Plate <- sub("_.*", "", colnames(met_body_SO))
met_body_SO$Sample <- sub("-.*", "", met_body_SO$Plate)

#install.packages("Matrix", type = "source")
#install_version("Matrix", version = "1.6-5", repos = "http://cran.us.r-project.org")
#install.packages("irlba", type = "source")

met_body_SO <- RunTFIDF(met_body_SO)
met_body_SO <- FindTopFeatures(met_body_SO, min.cutoff = 'q0')
met_body_SO <- RunSVD(met_body_SO)

met_body_SO <- RunUMAP(object = met_body_SO, reduction = 'lsi', dims = 2:20)
met_body_SO <- FindNeighbors(object = met_body_SO, reduction = 'lsi', dims = 2:20)
#met_body_SO <- FindClusters(object = met_body_SO, verbose = FALSE, algorithm = 3)

pdf(paste0(io$outDir, "/", data_name, "_metBodyUMAP.pdf"))
DimPlot(object = met_body_SO, group.by = "Plate")
dev.off()


acc_promoter_counts[which(is.na(acc_promoter_counts), arr.ind = T)] <- 0

acc_promoter_chrom_assay <- CreateChromatinAssay(
  counts = acc_promoter_counts,
  sep = c("-", "-"),
  min.cells = 5,
  min.features = 50
)

acc_promoter_SO <- CreateSeuratObject(
  counts = acc_promoter_chrom_assay,
  assay = "peaks")

acc_promoter_SO$Plate <- sub("_.*", "", colnames(acc_promoter_SO))
acc_promoter_SO$Sample <- sub("-.*", "", acc_promoter_SO$Plate)

acc_promoter_SO <- RunTFIDF(acc_promoter_SO)
acc_promoter_SO <- FindTopFeatures(acc_promoter_SO, min.cutoff = 'q0')
acc_promoter_SO <- RunSVD(acc_promoter_SO)

acc_promoter_SO <- RunUMAP(object = acc_promoter_SO, reduction = 'lsi', dims = 2:20)
acc_promoter_SO <- FindNeighbors(object = acc_promoter_SO, reduction = 'lsi', dims = 2:20)
#acc_promoter_SO <- FindClusters(object = acc_promoter_SO, verbose = FALSE, algorithm = 3)

pdf(paste0(io$outDir, "/", data_name, "_accPromoterUMAP.pdf"))
DimPlot(object = acc_promoter_SO, group.by = "Plate")
dev.off()

#Idents(acc_promoter_SO) <- "Sample"

#patient.markers <- FindAllMarkers(acc_promoter_SO, only.pos = T, test.use = 'LR', latent.vars = 'nCount_peaks')

#top.markers <- patient.markers %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)

#acc_promoter_SO <- ScaleData(acc_promoter_SO)

#pdf("allPatients_accDAheatmap.pdf")
#DoHeatmap(acc_promoter_SO, features = top.markers$gene) + theme(text = element_text(size = 7))
#dev.off()

#write.table(patient.markers, file = "allPatients_accDAlist.tsv", sep = "\t", row.names = F, quote = F)

met_promoter_counts[which(is.na(met_promoter_counts), arr.ind = T)] <- 0

met_promoter_chrom_assay <- CreateChromatinAssay(
  counts = met_promoter_counts,
  sep = c("-", "-"),
  min.cells = 5,
  min.features = 50
)

met_promoter_SO <- CreateSeuratObject(
  counts = met_promoter_chrom_assay,
  assay = "peaks")

met_promoter_SO$Plate <- sub("_.*", "", colnames(met_promoter_SO))
met_promoter_SO$Sample <- sub("-.*", "", met_promoter_SO$Plate)

#install.packages("Matrix", type = "source")
#install_version("Matrix", version = "1.6-5", repos = "http://cran.us.r-project.org")
#install.packages("irlba", type = "source")

met_promoter_SO <- RunTFIDF(met_promoter_SO)
met_promoter_SO <- FindTopFeatures(met_promoter_SO, min.cutoff = 'q0')
met_promoter_SO <- RunSVD(met_promoter_SO)

met_promoter_SO <- RunUMAP(object = met_promoter_SO, reduction = 'lsi', dims = 2:20)
met_promoter_SO <- FindNeighbors(object = met_promoter_SO, reduction = 'lsi', dims = 2:20)
#met_promoter_SO <- FindClusters(object = met_promoter_SO, verbose = FALSE, algorithm = 3)

pdf(paste0(io$outDir, "/", data_name, "_metPromoterUMAP.pdf"))
DimPlot(object = met_promoter_SO, group.by = "Plate")
dev.off()


acc_calledPeaks_counts[which(is.na(acc_calledPeaks_counts), arr.ind = T)] <- 0

acc_calledPeaks_chrom_assay <- CreateChromatinAssay(
  counts = acc_calledPeaks_counts,
  sep = c("-", "-"),
  min.cells = 5,
  min.features = 50
)

acc_calledPeaks_SO <- CreateSeuratObject(
  counts = acc_calledPeaks_chrom_assay,
  assay = "peaks")

acc_calledPeaks_SO$Plate <- sub("_.*", "", colnames(acc_calledPeaks_SO))
acc_calledPeaks_SO$Sample <- sub("-.*", "", acc_calledPeaks_SO$Plate)

#install.packages("Matrix", type = "source")
#install_version("Matrix", version = "1.6-5", repos = "http://cran.us.r-project.org")
#install.packages("irlba", type = "source")

acc_calledPeaks_SO <- RunTFIDF(acc_calledPeaks_SO)
acc_calledPeaks_SO <- FindTopFeatures(acc_calledPeaks_SO, min.cutoff = 'q0')
acc_calledPeaks_SO <- RunSVD(acc_calledPeaks_SO)

acc_calledPeaks_SO <- RunUMAP(object = acc_calledPeaks_SO, reduction = 'lsi', dims = 2:20)
acc_calledPeaks_SO <- FindNeighbors(object = acc_calledPeaks_SO, reduction = 'lsi', dims = 2:20)
#acc_calledPeaks_SO <- FindClusters(object = acc_calledPeaks_SO, verbose = FALSE, algorithm = 3)

pdf(paste0(io$outDir, "/", data_name, "_accPeaksUMAP.pdf"))
DimPlot(object = acc_calledPeaks_SO, group.by = "Plate")
dev.off()

#Idents(acc_calledPeaks_SO) <- "Sample"

#patient.markers <- FindAllMarkers(acc_calledPeaks_SO, only.pos = T, test.use = 'LR', latent.vars = 'nCount_peaks')

#top.markers <- patient.markers %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)

#acc_calledPeaks_SO <- ScaleData(acc_calledPeaks_SO)

#pdf("allPatients_accDAheatmap.pdf")
#DoHeatmap(acc_calledPeaks_SO, features = top.markers$gene) + theme(text = element_text(size = 7))
#dev.off()

#write.table(patient.markers, file = "allPatients_accDAlist.tsv", sep = "\t", row.names = F, quote = F)

met_calledPeaks_counts[which(is.na(met_calledPeaks_counts), arr.ind = T)] <- 0

met_calledPeaks_chrom_assay <- CreateChromatinAssay(
  counts = met_calledPeaks_counts,
  sep = c("-", "-"),
  min.cells = 5,
  min.features = 50
)

met_calledPeaks_SO <- CreateSeuratObject(
  counts = met_calledPeaks_chrom_assay,
  assay = "peaks")

met_calledPeaks_SO$Plate <- sub("_.*", "", colnames(met_calledPeaks_SO))
met_calledPeaks_SO$Sample <- sub("-.*", "", met_calledPeaks_SO$Plate)

#install.packages("Matrix", type = "source")
#install_version("Matrix", version = "1.6-5", repos = "http://cran.us.r-project.org")
#install.packages("irlba", type = "source")

met_calledPeaks_SO <- RunTFIDF(met_calledPeaks_SO)
met_calledPeaks_SO <- FindTopFeatures(met_calledPeaks_SO, min.cutoff = 'q0')
met_calledPeaks_SO <- RunSVD(met_calledPeaks_SO)

met_calledPeaks_SO <- RunUMAP(object = met_calledPeaks_SO, reduction = 'lsi', dims = 2:20)
met_calledPeaks_SO <- FindNeighbors(object = met_calledPeaks_SO, reduction = 'lsi', dims = 2:20)
#met_calledPeaks_SO <- FindClusters(object = met_calledPeaks_SO, verbose = FALSE, algorithm = 3)

pdf(paste0(io$outDir, "/", data_name, "_metPeaksUMAP.pdf"))
DimPlot(object = met_calledPeaks_SO, group.by = "Plate")
dev.off()
