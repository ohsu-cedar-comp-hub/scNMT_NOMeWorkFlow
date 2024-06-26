args <- commandArgs()

help <- function(){
    cat("annotate_arw_met.R :
quantify methylation rates over genomic annotations
input is raw methylation data (1 file per cell)
output is quantified methylation data (1 file per annotation)\n")
    cat("Usage: \n")
    cat("--anno        : path to annotation files                               [required]\n")
    cat("--raw         : path to raw binarised methylation/accessibility files  [required]\n")
    cat("--meta        : stats table with qc information                        [required]\n")
    cat("--outdir      : output directory                                       [required]\n")
    cat("\n")
    q()
}

io   <- list()
opts <- list()

## Save values of each argument
if( !is.na(charmatch("--help",args)) || !is.na(charmatch("--help",args)) ){
    help()
} else {
    io$anno_dir   <- sub( '--anno=', '', args[grep('--anno', args)] )
    io$raw_files  <- sub( '--raw=', '', args[grep('--raw', args)] )
    io$meta_data  <- sub( '--meta=', '', args[grep('--meta', args)] )
    io$out_dir    <- sub( '--outdir=', '', args[grep('--outdir', args)] )
}

library(data.table)
library(purrr)
library(furrr)


# quantify methylation rates over genomic annotations
# input is raw methylation data (1 file per cell)
# output is quantified methylation data (1 file per annotation)


### TEST INPUT ###
#setwd("../")
#io$anno_dir <- "data/anno"
#io$raw_files <- "bismarkSE/CX/coverage2cytosine_1based/filt/binarised"
#io$meta_data <- "tables/sample_stats_qcPass.txt"
#io$out_dir <- "data/met"

dir.create(io$out_dir, recursive = TRUE)

opts$anno_regex <- "promoters1000.bed|MCF7_ER_peaks|H3K27ac_peaks|body1000.bed|Repressed|Enhancer|CTCF|called_peaks" # use to select specific annotations
opts$parallel <- FALSE
opts$cores <- 2
opts$gzip <- TRUE
opts$extend_anno_len <- FALSE
opts$anno_len <- 500


### functions ###
fread_gz = function(filename, ...){
  f <- file(filename)
  type <- summary(f)$class
  close.connection(f)
  if (type == "gzfile") return(fread(cmd = paste("zcat", filename), ...))
  fread(filename, ...)
}


### load metadata and select cells/files ###

meta <- fread(io$meta_data)
meta  <- meta[context == "CG" & pass_accQC == TRUE & pass_CHGQC == TRUE & pass_CHHQC == TRUE]
print(head(meta))

cells <- unique(meta[, sample])
files <- unique(meta[, paste0(io$raw_files, "/", sample, "_CpG", ".gz")])
#} else {
#  meta <- meta[pass_accQC == TRUE]
#  cells <- meta[, sample]
#  files <- meta[, paste0(io$raw_files, "/", sample, "_CpG", ".gz")]
#}

cells <- cells[file.exists(files)]
files <- files[file.exists(files)]
print(head(cells))
print(head(files))

### load annotations ###

anno <- dir(io$anno_dir, pattern = ".bed", full = TRUE) %>%
  .[grep(opts$anno_regex, .)]%>%
#    map2(., sub(".bed", "", basename(.)), ~fread(.x, colClasses = list("factor" = 1L)) %>% 
#           setnames(c("chr", "start", "end", "strand", "id", "anno")) %>% 
#           .[, anno := paste0(anno, "_", .y)]) %>%
    map(fread, colClasses = list("character" = 1)) %>%
    rbindlist() %>% setnames(c("chr", "start", "end", "strand", "id", "anno")) %>%
    .[,chr := gsub("chr", "", chr)] %>% .[!chr == "M" & !chr == "Y"] %>%
setkey(chr, start, end)


if (opts$extend_anno_len) {
  anno <- anno[, mid := start + (end-start)/2] 
  anno[end-start < opts$anno_len, c("start", "end") := .(as.integer(round(mid - opts$anno_len/2), round(mid + opts$anno_len/2)))]
  setkey(anno, chr, start, end)
}

print(unique(anno$anno))

### perform overlap and quantification ###


if (opts$parallel) {
  plan(multiprocess, workers = opts$cores)
} else {
  plan(sequential)
}

acc_dt <- future_map2(files, cells, ~{
#acc_dt <- map2(files, cells, ~{
  if (!file.exists(.x)) return(NULL)
  #fread(.x, select=c(1:2,5), colClasses = list("factor" = 1L)) %>%
  fread_gz(.x, select = c(1:2, 5), colClasses = list("factor" = 1L)) %>%
    .[, .(chr = gsub("chr", "", chr),
    #.[, .(chr = chr,
          start = pos,
          end = pos,
          rate)] %>%
    setkey(chr, start, end) %>%
    foverlaps(anno, nomatch = 0L) %>% 
    .[, .(rate = round(100 * mean(rate)), Nmet = sum(rate == 1), .N, sample = .y, var=var(rate)), .(id, anno)]
}) %>%
  purrr::compact() %>%
  rbindlist()

print(unique(acc_dt$anno))

setkey(acc_dt, anno)
acc_dt <- split(acc_dt, by = "anno")

# save one file per anno
file_names <- paste0(io$out_dir, "/", names(acc_dt), ".tsv")
if (opts$extend_anno_len) file_names <- sub(".tsv", "_", opts$anno_len, "bp.tsv")

# future_map2(acc_dt, file_names, fwrite, sep = "\t")
map2(acc_dt, file_names, fwrite, sep = "\t")

if (opts$gzip) walk(file_names, ~system(paste("gzip -f", .)))
