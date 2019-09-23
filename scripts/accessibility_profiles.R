args <- commandArgs()

help <- function(){
    cat("accessibility_profiles.R :
Make profiles plots around promoters. Output will be save in the parsed directory in the covPath. The promter and body regions (gene body without promoter) will be saved in the outDir.\n")
    cat("Usage: \n")
    cat("--covPath     : path to acc or met data (tsv.gz)                    [required]
                          format: <chr pos met_reads nonmet_reads rate>            \n")
    cat("--outDir      : output directory                                    [required]\n")
    cat("--WinSize     : size of window to average                           [default = 100]\n")
    cat("--StepSize    : slide start of window by step size                  [default = 20]\n")    
    cat("--promUp      : number of nucleotides upstream of the tss           [default = 2000 ]
                           will also be used for downstream                               \n")
    cat("--annoFile    : tsv <ens_id gene chr start end strand >             [required]\n")
    cat("--qcFile      : tsv required feilds:                                [required]
                           <id, context, passmet_QC, pass_accQC>
                           id must be just the name of the sample no context
                           context must either be GC or CG                             \n")
    cat("--context     : either GC or CG                                     [required]
                          assumes files end in _CpG.tsv.gz or _CpG.tsv.gz          \n")
    cat("--covCutOff   : coverage cutoff for high quality cells              [optional]\n")
    cat("\n")
    q()
}

io   <- list()
opts <- list()

## Save values of each argument
if( !is.na(charmatch("--help",args)) || !is.na(charmatch("--help",args)) ){
    help()
} else {
    io$raw_acc   <- sub( '--covPath=', '', args[grep('--covPath=', args)] )
    WinSize      <- sub( '--WinSize=', '', args[grep('--WinSize=', args)] )
    promUp       <- sub( '--promUp=', '', args[grep('--promUp=', args)] )
    StepSize     <- sub( '--StepSize=', '', args[grep('--StepSize=', args)])
    io$anno      <- sub( '--annoFile=', '', args[grep('--annoFile=', args)])
    io$qc_file   <- sub( '--qcFile=', '', args[grep('--qcFile=', args)])
    opts$context <- sub( '--context=', '',args[grep('--context=',args)])
    io$out_dir   <- sub( '--outDir=', '',args[grep('--outDir=',args)])
    covCutOff    <- sub( '--covCutOff=', '',args[grep('--covCutOff=',args)])
    

}

if (identical(StepSize,character(0))){
    opts$covCutOff <- 0
}else{
    opts$covCutOff <- as.numeric(covCutOff)
}

if(!(file.exists( io$out_dir ))) {
    print(paste("mkdir", io$out_dir))
    dir.create(io$out_dir,FALSE,TRUE)  
}

if(!(file.exists( paste(io$out_dir, "anno", sep="/") ))) {
    print(paste("mkdir", paste(io$out_dir, "anno", sep="/")))
    dir.create( paste(io$out_dir, "anno", sep="/") ,FALSE,TRUE)  
}

if (identical(WinSize,character(0))){
   opts$win <- 100
}else{
   opts$win <- as.numeric(WinSize)
}

if (identical(StepSize,character(0))){
   opts$step <- 20
}else{
   opts$step <- as.numeric(StepSize)
}

if (identical(promUp,character(0))){
    opts$promUp   <- 2000
    opts$promDown <- 2000
}else{
    opts$promUp   <- as.numeric(promUp)
    opts$promDown <- as.numeric(promUp)
}

if(opts$context == "GC"){
    opts$outName <- "accessibility"
}else{
    opts$outName <- "methylation"
}

library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)
library(GenomicRanges)
library(dplyr)
library(GenomeInfoDb)

## debugging
#io         <- list()
#io$raw_acc <- "/home/groups/CEDAR/woodfin/projects/NMT-seq/20190523_NM/bismarkSE/CX/coverage2cytosine_1based/filt/binarised"
#io$qc_file <- "/home/groups/CEDAR/woodfin/projects/NMT-seq/20190523_NM/tables/sample_stats_qcPass.txt"
#io$anno    <- "/home/groups/CEDAR/woodfin/projects/NMT-seq/20190523_NM/anno/hg19/gene_hg19.cellRanger_metadata.tsv"
#opts <- list()
#opts$win <- 100
#opts$step <- 20
#opts$context <- "CG"
#opts$promUp <- 2000
#io$out_dir <- "/home/groups/CEDAR/woodfin/projects/NMT-seq/20190523_NM/plots/profiles"
#opts$covCutOff <- 1e7

print(io)
print(opts)

meta <- fread(io$qc_file) %>%
    .[pass_accQC == TRUE & pass_metQC == TRUE]
print(meta)

cells <- unique(meta[, id])

if( opts$context == "GC" ){
    files <- paste0(io$raw_acc,  "/",cells, "_GpC.tsv.gz") %>%
        .[file.exists(.)]
}else{
    files <- paste0(io$raw_acc,  "/",cells, "_CpG.tsv.gz") %>%
        .[file.exists(.)]
}    
head(files)

prom <- fread(io$anno) %>%
  .[strand == "+", tss := start] %>%
  .[strand == "-", tss := end] %>% .[!chr %in% c("MT", "Y")]
print(prom)

gr <- as(prom, "GRanges") %>% .[width(.) > opts$promUp]
print(gr)

# save gene body info for correlations
body      <- resize(gr, fix='end', width=width(gr)-opts$promUp )
body$anno <- "body"
body      <- as.data.table(body)
df        <- body[,.(seqnames, start, end, strand, ens_id, anno)]
print(head(df))

fwrite(df
     , paste0(paste(io$out_dir, "anno", sep="/"), "/", "body.bed"), sep="\t", col.names = F, row.names = F, quote=F, na="NA")

# save promoter info
prom <- as.data.table(promoters(as(prom, "GRanges"), upstream=2000, downstream=2000))

df   <- prom[,.(seqnames, start, end, strand, ens_id, "promoter")]
fwrite(df[!seqnames %in% c("MT", "Y"), ]
     , paste0(paste(io$out_dir, "anno", sep="/"), "/", "promters.bed"), sep="\t", col.names = F, row.names = F, quote=F, na="NA")

##############
# make windows
##############
boo        <- as(prom, "GRanges")
names(boo) <- boo$ens_id

#hmm <- unlist(slidingWindows(boo, 50, step=50))
hmm        <- unlist(slidingWindows(boo, opts$win, step=opts$step))
hmm$ens_id <- sub("\\..*", "", names(hmm))

# get anno for windows
iv       <- match(hmm$ens_id, prom$ens_id)
hmm$gene <- prom[iv,gene]
hmm$tss <- prom[iv,tss]

# 
prom        <- as.data.table(hmm)
names(prom) <- sub("seqnames", "chr", names(prom))

prom[strand == "+", dist := tss - start]
prom[strand == "-", dist := start - tss]

anno <- prom[,c("chr", "start", "end", "gene", "tss", "strand", "dist")]
anno <- anno[!chr %in% c("MT","Y")]
setkey(anno, chr, start, end)

acc <- map2(cells, files, ~{
  cell = .x
  d=fread(.y) %>% 
    .[, .(chr = gsub("^chr", "", chr),
          start = pos,
          end = pos,
          rate)] %>%
    setkey(chr, start, end) %>%
    foverlaps(anno, nomatch = 0L) %>%
    .[, .(rate = mean(rate), cell = cell), dist]
}) %>%
  rbindlist()

acc <- acc[order(acc$cell, acc$dist),]
print(head(acc))
save(acc,file=paste(io$out_dir, paste0(opts$outName, "_at_promoters.RData"), sep="/"))


is.odd <- function(x) x %% 2 != 0
is.even <- function(x) x %% 2 == 0

Avg           <- acc[seq(1,nrow(acc), by=1) %>% is.odd(), ] %>% .[,rate:=0]
Rate <- vector()
for( i in seq(1,nrow(acc)-1,by=2)){
  tmp         <- acc[c(i,i+1)] 
  Rate <- mean(tmp$rate)
}

Rate <- as.data.frame(do.call(rbind,lapply(seq(1,nrow(acc)-1,by=2), function(i){ 
  tmp         <- acc[c(i,i+1)]
  mean(tmp$rate)
})))

Avg$rate <- Rate[,1]

p <- ggplot(Avg, aes(dist, rate, colour = cell)) +
    geom_line() +
    theme_bw() +
    guides(colour = FALSE)
p

save_plot(paste(io$out_dir, paste0(opts$outName, "_at_promoters.pdf"),sep="/"), p)

### repeat with only cells passing stricter threshold ###
if(opts$covCutOff > 0){
    qc <- fread(io$qc_file) %>%
        .[, sample := strsplit(sample, "_") %>% map_chr(1)]
    top_cells <- qc[context ==  opts$context & coverage > opts$covCutOff, id]
    sub <-  Avg[cell %in% top_cells]
    p2 <- ggplot(sub, aes(dist, rate, colour = cell)) +
        theme_bw() +
        geom_line() 
    p2
    pdf(paste(io$out_dir, paste0(opts$outName, "_at_promoters_top_cells.pdf"),sep="/"), 8, 8)
    print(p2)
    dev.off()
}
