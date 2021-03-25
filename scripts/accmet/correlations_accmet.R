args <- commandArgs()

help <- function() {
    cat("plotRatesByAnno.R :
Plots rates across the different annotations\n")
    cat("Usage: \n")
    cat("--meta        : stats table with QC info     [required]\n")
    cat("--met         : path to methylation files    [required]\n")
    cat("--acc         : path to accessibility files  [required]\n")
    cat("--plotdir     : location to output plots     [required]\n")
    cat("--anno        : path to annotation files     [required]\n")
    cat("--genes       : path to gene annotation      [required]\n")
    cat("--min_cells_met : min number cells covered (20 or 5)      [required]
                        script dies if you do not meet this min!\n")
    cat("--min_cells_acc : min number cells covered (20 or 5)      [required]\n")
    cat("add option to turn on correlation of cells\n")
    cat("\n")
    q()
}

io <- list()
opts <- list()

print(args)

## Save values of each argument
if( !is.na(charmatch("--help",args)) || !is.na(charmatch("--help",args)) ){
    help()
} else {
    io$nome_meta_data  <- sub( '--meta=', '', args[grep('--meta=', args)] )
    io$met_dir    <- sub( '--met=', '', args[grep('--met=', args)] )
    io$acc_dir    <- sub( '--acc=', '', args[grep('--acc=', args)] )
    io$plot_dir   <- sub( '--plotdir=', '', args[grep('--plotdir=', args)] )
    io$anno_dir   <- sub( '--anno=', '', args[grep('--anno=', args)] )
    io$gene_file  <- sub( '--genes=', '', args[grep('--genes=', args)] )
    opts$min_cells_met <- sub( '--min_cells_met=', '', args[grep('--min_cells_met=', args)] )
    opts$min_cells_acc <- sub( '--min_cells_acc=', '', args[grep('--min_cells_acc', args)] )
}

if ("weights" %in% rownames(installed.packages()) == FALSE) {
    install.packages("weights", repos="https://ftp.osuosl.org/pub/cran/")
}

#library(SingleCellExperiment)
library(scater)
library(data.table)
library(purrr)
library(weights)
library(ggplot2)
library(cowplot)
library(ggrepel)

##### TEST INPUT #####
#setwd("../")
#io <- list()
#opts <- list()
#io$nome_meta_data <- "tables/sample_stats_qcPass.txt"
#io$acc_dir        <- "data/acc/"
#io$met_dir        <- "data/met/"
#io$anno_dir       <- "data/anno/"
#io$plot_dir       <- "plots/corr"
#io$gene_file      <- "data/gene_metadata.tsv"
#opts$min_cells_met <- 5 # loci must have observations in this many cells
#opts$min_cells_acc <- 5 # loci must have observations in this many cells

opts$anno_regex <- "promoter|MCF7_ER_peaks|H3K27ac_peaks|body10000.bed|Repressed|Enhancer|CTCF"
opts$gene_overlap_dist <- 1e5 # overlap annoations with genes within xx bp
opts$min_weight_met <- 1
opts$min_weight_acc <- 1
#opts$min_cells_met <- 20 # loci must have observations in this many cells
#opts$min_cells_acc <- 20 # loci must have observations in this many cells
opts$min.s         <- 5     # minimum number of samples to do the correlation

opts$filt_acc_var <- 0.5 # select the top xx fraction by variance
opts$filt_met_var <- 0.5 # select the top xx fraction by variance

opts$min.coverage   <- 0.5    # minimum coverage per feature across samples (met/acc)
opts$fraction.sites <- 0.5    # fraction of sites (met/acc) to keep based on variance

opts$p_cutoff <- 0.05
opts$cor_samples <- F
opts$weight <- T
opts$method <- "pearson"      # correlation type
#opts$method <- ""      # correlation type

print(io)
print(opts)

### functions ###
fread_gz <- function(path, ...){fread(cmd = paste("zcat", path), ...)}

fread_gz <- function(path, ...){fread(x, ...)}

### load metadata and select cells ####
meta <- fread(io$nome_meta_data)
meta <- meta %>%
  .[pass_accQC == TRUE & pass_metQC == TRUE & pass_CHGQC == TRUE & pass_CHHQC == TRUE]

### load met data ###
met <- dir(io$met_dir, pattern = ".tsv.gz$", full = TRUE) %>%
  .[grep(opts$anno_regex, .)] %>%
  # map(fread_gz) %>%
  map(fread) %>%
  rbindlist()

### load acc data ###
acc <- dir(io$acc_dir, pattern = ".tsv.gz$", full = TRUE) %>%
  .[grep(opts$anno_regex, .)] %>%
  # map(fread_gz) %>%
  map(fread) %>%
  rbindlist()

### merge acc met data ###
metacc_dt <-merge(
  met[,c("sample", "id", "anno", "rate", "N")] %>% setnames(c("rate", "N"), c("met_rate", "Nmet")),
  acc[,c("sample", "id", "anno", "rate", "N")] %>% setnames(c("rate", "N"), c("acc_rate", "Nacc")),
  by=c("sample", "id", "anno")
)

### Remove features with low weight (all types of correlation) ###
metacc_dt <- metacc_dt[Nmet>=opts$min_weight_met & Nacc>=opts$min_weight_acc]
print(head(metacc_dt))

### annotate wiht the nearest gene
genes <- fread(io$gene_file) %>%
  .[, c("start", "end", "chr") := .(start - opts$gene_overlap_dist, 
                                    end + opts$gene_overlap_dist,
                                    gsub("chr", "", chr))] %>%
  setkey(chr, start, end)

anno <- dir(io$anno_dir, full = TRUE, pattern = ".bed$") %>%
  .[grep(opts$anno_regex, .)] %>%
  map(fread) %>%
  rbindlist() %>% setnames(c("chr", "start", "end", "strand", "id", "anno")) %>%
  .[, .(chr = gsub("chr", "", chr),
        start = start,
        end = end,
        strand = strand,
        id = id,
        anno = anno
  )] %>% .[!chr %in% c("MT", "Y", "M")] %>%
  .[, gene_id := grepl("ENSG", id)] %>%
  split(by = "gene_id", keep.by = FALSE)

anno[["TRUE"]] <- anno[["TRUE"]] %>%
  .[, ens_id := id] %>%
  merge(genes[, .(ens_id, gene)], by = "ens_id")

anno[["FALSE"]] <- setkey(anno[["FALSE"]], chr, start, end) %>%
  foverlaps(genes, nomatch = 0L)

anno <- map(anno, ~.[, .(id, anno, gene, ens_id)]) %>%
  rbindlist()

print(head(anno))

metacc_dt <- merge(metacc_dt, anno, by = c("anno", "id"), allow.cartesian = TRUE) # note some loci have >1 gene -> cartesian join
print(head(metacc_dt))
print(unique(metacc_dt$anno))
### filter ###

## remove features with small number of cells
keep_cov_sites <- metacc_dt %>%
    split(.$anno) %>%
    map(~ .[,.(n=.N), by=c("gene","id")] %>%
            .[n>=opts$min_cells_met] %>% .$id)

metacc_dt <- metacc_dt %>% split(.$anno) %>% 
  map2(., names(.), function(x,y) x[id %in% keep_cov_sites[[y]]]) %>% rbindlist

# Remove genomic features with no variability 
keep_var_sites.met <- metacc_dt %>% split(.$anno) %>% map(~ .[,.(var=var(met_rate)), by="id"] %>% .[var>0, id] %>% as.character())
keep_var_sites.acc <- metacc_dt %>% split(.$anno) %>% map(~ .[,.(var=var(acc_rate)), by="id"] %>% .[var>0, id] %>% as.character())
metacc_dt <- metacc_dt %>% split(.$anno) %>%
  map2(., names(.), function(x,y) x[id %in% intersect(keep_var_sites.met[[y]],keep_var_sites.acc[[y]])]) %>% rbindlist

if (opts$cor_samples) {
  
  # Intersect the two data sets
  metacc <- merge(met[,c("sample", "id", "anno", "rate", "N")] %>% setnames(c("rate", "N"), c("met_rate", "met_weight")),
                  acc[,c("sample", "id", "anno", "rate", "N")] %>% setnames(c("rate", "N"), c("acc_rate", "acc_weight")),
                  by=c("sample", "id", "anno"))
  ## Remove features with low weight (all types of correlation)
  metacc <- metacc[met_weight >= opts$min_weight_met & acc_weight >= opts$min_weight_acc]
  # To correlate across samples
  metacc_filt <- copy(metacc)
  
  ## Filter sites with low coverage
  nsamples <- length(unique(metacc$sample))
  metacc_filt <- metacc_filt[, cov := .N / nsamples, by = c("id", "anno")] %>% .[cov >= opts$min.coverage] %>% .[, cov := NULL]
  metacc <- metacc[, cov := .N / nsamples, by = c("id", "anno")] %>% .[cov >= opts$min.coverage] %>% .[, cov := NULL]
  
  
  ## Remove constant sites and filter based on variability (separately for each feature)
  keep_hv_sites <- metacc_filt %>% split(.$anno) %>% map(~ .[,.(met_var = wtd.var(met_rate, met_weight), acc_var = wtd.var(acc_rate, acc_weight)), by = c("id")] %>% .[met_var > 2 | acc_var > 2] %>% .[, var := acc_var * met_var] %>% setorder(-var)  %>% head(n = nrow(.) * opts$fraction.sites) %>% .$id)
  metacc_filt <- metacc_filt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist
  
  ## Filter id pairs with small number of samples to do the correlation
  metacc_filt <- metacc_filt[,n:=.N, by=c("id","anno")] %>% .[n >= opts$min.s] %>% .[,n:=NULL]
}

### compute correlation coeff ###
weight_cor <- function(met_rate, acc_rate, N){
  wtd.cor(met_rate, acc_rate, N) %>%
    as.list() %>%
    set_names(c("r", "std_err", "t", "p"))
}

method_cor <- function(met_rate, acc_rate){
  cor.test(met_rate, acc_rate, alternative = "two.sided", method = opts$method) %>%
    as.list() %>%
    set_names(c("r", "std_err", "t", "p"))
}

# Weighted correlation --> using met weight because more strict
if (opts$weight == TRUE){
  print("weighted approach")
  if (opts$method != "pearson") { print("Weighted correlation only supported for pearson"); stop() }
  if (opts$cor_samples) {
    # Correlate rate across samples
      cor_samples <- metacc_filt[, weight_cor(met_rate, acc_rate, met_weight), .(anno, id)] %>%
      .[, padj := p.adjust(p, method = "fdr"), .(anno)] %>%
      .[, logpadj := -log10(padj)] %>%
      .[, sig := padj < opts$p_cutoff]
  }
  # Correlate rate across genes
  cor_features <- metacc_dt[, weight_cor(met_rate, acc_rate, Nmet), .(anno, id, gene, ens_id)] %>%
    .[, padj := p.adjust(p, method = "fdr"), .(anno)] %>%
    .[, logpadj := -log10(padj)] %>%
    .[, sig := padj < opts$p_cutoff]
}else{
  if (opts$cor_samples) {
    print("non weighted approach used and needs to be fixed")
     # currently pearson method just gives V1 and need to fix
    cor_samples <- metacc_filt[, .(V1 = unlist(cor.test(met_rate, acc_rate, alternative = "two.sided", method = opts$method)[c("estimate", "statistic", "p.value")])), by = c("id", "anno")]
    cor_samples <- metacc_filt[, method_cor(met_rate, acc_rate), .(anno, id)] %>%
      .[, padj := p.adjust(p, method = "fdr"), .(anno)] %>%
      .[, logpadj := -log10(padj)] %>%
      .[, sig := padj < opts$p_cutoff]
  }  
  # Correlate rate across genes
  cor_features <- metacc_dt[, .(V1 = unlist(cor.test(met_rate, acc_rate, alternative = "two.sided", method = opts$method)[c("estimate", "statistic", "p.value")])), by = c("sample", "anno")]
}

print(head(cor_features))
###################
## plot results across loci
###################

labs <- c("NS", "Significant")

p <- ggplot(cor_features, aes(r, logpadj, colour = sig)) +
  geom_point() +
  facet_wrap(~anno) +
  scale_colour_manual(values = c("grey", "navy"), labels = labs, name = NULL) +
  geom_hline(yintercept = -log10(opts$p_cutoff), colour = "blue", linetype = "dashed") +
  geom_vline(aes(xintercept = mean(r)), colour = "blue") +
  #geom_text_repel(data = cors[sig == TRUE]) +
  labs(x = c("Weighted Pearson R"), y = "-log10 q-value") +
  ggtitle("Correlations across loci")+
  guides(label = FALSE) +
  theme_bw()
p

if(!exists(io$plot_dir)){
  dir.create(io$plot_dir)
}

save_plot(paste(io$plot_dir, "acc_met_correlations_loci.pdf", sep="/"), p, base_width = 12, base_height = 8)

fwrite(cor_features
       , paste(io$plot_dir, "acc_met_correlations_loci.tsv", sep="/"), sep="\t", 
       col.names = T, row.names = F, quote=F, na="NA")

###################
## plot results across samples
###################

labs <- c("NS", "Significant")

if (opts$cor_samples) {
    print("make sample cor plot")
    p <- ggplot(cor_samples, aes(r, logpadj, colour = sig)) +
        geom_point() +
        facet_wrap(~anno) +
        scale_colour_manual(values = c("grey", "navy"), labels = labs, name = NULL) +
        geom_hline(yintercept = -log10(opts$p_cutoff), colour = "blue", linetype = "dashed") +
        geom_vline(aes(xintercept = mean(r)), colour = "blue") +
                                        #geom_text_repel(data = cors[sig == TRUE]) +
        labs(x = c("Weighted Pearson R"), y = "-log10 q-value") +
        ggtitle("Correlations across samples")+
        guides(label = FALSE) +
        theme_bw()
    save_plot(paste(io$plot_dir, "acc_met_correlations_samplesAnno.pdf", sep="/"), p, base_width = 12, base_height = 8)
    dev.off()
    fwrite(cor_samples
         , paste(io$plot_dir, "acc_met_correlations_samples.tsv", sep="/"), sep="\t", col.names = T, row.names = F, quote=F, na="NA")
}




