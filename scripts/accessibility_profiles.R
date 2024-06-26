args <- commandArgs()

help <- function(){
    cat("accessibility_profiles.R :
Make profiles plots around promoters. Output will be save in the parsed directory in the covPath. The promter and body regions (gene body without promoter) will be saved in the outDir.\n")
    cat("Usage: \n")
    cat("--covPath     : path to acc or met data (tsv.gz)                    [required]
                         format: <chr pos met_reads nonmet_reads rate>                 \n")
    cat("--outDir      : output directory                                    [required]\n")
    cat("--WinSize     : size of window to average                           [default = 100]\n")
    cat("--StepSize    : slide start of window by step size                  [default = 20]\n")    
    cat("--promBed     : path to windows around tss (data/promoters1000.bed) [required]\n")
    cat("--annoFile    : tsv <ens_id gene chr start end strand >             [required]\n")
    cat("--qcFile      : tsv required feilds:                                [required]
                         <id, context, passmet_QC, pass_accQC>
                         id must be just the name of the sample no context
                         context must either be GC or CG                               \n")
    cat("--context     : either GC or CG                                     [required]
                         assumes files end in _CpG.tsv.gz or _CpG.tsv.gz               \n")
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
    io$promBed   <- sub( '--promBed=', '', args[grep('--promBed=', args)] )
    StepSize     <- sub( '--StepSize=', '', args[grep('--StepSize=', args)])
    io$anno      <- sub( '--annoFile=', '', args[grep('--annoFile=', args)])
    io$qc_file   <- sub( '--qcFile=', '', args[grep('--qcFile=', args)])
    opts$context <- sub( '--context=', '',args[grep('--context=',args)])
    io$out_dir   <- sub( '--plotDir=', '',args[grep('--plotDir=',args)])
    covCutOff    <- sub( '--covCutOff=', '',args[grep('--covCutOff=',args)])
    io$data_dir  <- sub( '--dataDir=', '',args[grep('--dataDir=',args)])
    

}

######################
### Test Input Only###
######################
#setwd("../")
#io   <- list()
#opts <- list()
#io$raw_acc <- "bismarkSE/CX/coverage2cytosine_1based/filt/binarised"
#io$out_dir <- "data"
#io$anno <- "data/gene_metadata.tsv"
#io$qc_file <- "tables/sample_stats_qcPass.txt"
#io$promBed <- "data/anno/promoters1000.bed"
#io$bodyBed <- "data/body1000.bed"
#opts$context <- "GC"
#opts$win <- 100
#opts$step <- 20
#opts$covCutOff <- 0
#io$out_dir <- "plots/profiles"
#io$plot_dir <- "plots/profiles"
#io$data_dir <- "data"

##################

if (identical(covCutOff,character(0))){
    opts$covCutOff <- 0
}else{
    opts$covCutOff <- as.numeric(covCutOff)
}

if(!file.exists( io$out_dir )) {
    print(paste("mkdir", io$out_dir))
    dir.create(io$out_dir,FALSE,TRUE)  
}

if(!file.exists( io$data_dir )) {
    print(paste("mkdir", io$data_dir))
    dir.create(io$data_dir,FALSE,TRUE)  
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

#if (identical(promUp,character(0))){
#    opts$promUp   <- 2000
#    opts$promDown <- 2000
#}else{
#    opts$promUp   <- as.numeric(promUp)
#    opts$promDown <- as.numeric(promUp)
#}

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
library(stringr)

print(io)
print(opts)

meta <- fread(io$qc_file) %>%
    .[pass_accQC == TRUE & pass_metQC == TRUE & pass_CHGQC==TRUE & pass_CHHQC==TRUE]
print(meta)

cells <- unique(meta[, id])

if( opts$context == "GC" ){
    files <- paste0(io$raw_acc,  "/",cells, "_GpC.gz") %>%
        .[file.exists(.)]
}else{
    files <- paste0(io$raw_acc,  "/",cells, "_CpG.gz") %>%
        .[file.exists(.)]
}    
head(files)

#prom <- fread(io$anno) %>%
#  .[strand == "+", tss := start] %>%
#  .[strand == "-", tss := end] %>% .[!chr %in% c("MT", "Y")]
#print(prom)

#gr <- as(prom, "GRanges") %>% .[width(.) > opts$promUp]
                                        #print(gr)
prom <- fread(io$promBed) %>% setnames(c("chr","start", "end", "strand", "ens_id", "anno")) %>%
    .[strand == "+", tss := (start+end)/2] %>%                                                    #promoter files were built as windows around promoter, so tss is in the middle 
    .[strand == "-", tss := (start+end)/2] %>% .[!chr %in% c("chrM", "chrY", "M", "Y")]           #promoter files were built as windows around promoter, so tss is in the middle 
prom <- as(prom, "GRanges")
print(prom)

##############
# make windows
##############
boo        <- as(prom, "GRanges")
names(boo) <- boo$ens_id

#hmm <- unlist(slidingWindows(boo, 50, step=50))
hmm        <- unlist(slidingWindows(boo, opts$win, step=opts$step))
hmm$ens_id <- sub("\\..*", "", names(hmm))
print(hmm)

# get anno for windows
anno <- fread(io$anno)
iv   <- match(hmm$ens_id, anno$ens_id)
hmm$gene <- anno[iv,gene]

iv   <- match(hmm$ens_id, prom$ens_id)
hmm$tss <- as.data.frame(prom)[iv,"tss"]
print(hmm)

# 
prom        <- as.data.table(hmm)
names(prom) <- sub("seqnames", "chr", names(prom))

prom[strand == "+", dist := tss - start]
prom[strand == "-", dist := start - tss]

anno <- prom[,c("chr", "start", "end", "gene", "tss", "strand", "dist")]

# make sure neither have chr to avoid format issue
anno$chr <- sub("^chr", "", anno$chr) 
anno <- anno[!chr %in% c("MT","Y")]

setkey(anno, chr, start, end)

print("calculate rates for window around promoters")
acc <- map2(cells, files, ~{
  cell = .x
  d=fread(.y) %>% 
    .[, .(chr = gsub("^chr", "", chr),
          start = pos,
          end = pos,
          rate)] %>%
    setkey(chr, start, end) %>%
    foverlaps(anno, nomatch = 0L) %>%
    .[, .(sd_down = mean(rate)-sd(rate), sd_up = mean(rate)+sd(rate), rate = mean(rate), cell = cell), dist]
}) %>%
  rbindlist()

acc <- acc[order(acc$cell,acc$dist),]
print(head(acc))
print("save rates at promoters")
save(acc,file=paste(io$out_dir, paste0(opts$outName, "_at_promoters.RData"), sep="/"))


is.odd <- function(x) x %% 2 != 0
is.even <- function(x) x %% 2 == 0

### Finding Mean ###

Avg           <- acc[seq(1,nrow(acc), by=1) %>% is.even(), ] %>% .[,rate:=0]
#Rate <- vector()
#for( i in seq(1,nrow(acc)-1,by=2)){
#  tmp         <- acc[c(i,i+1)] 
#  Rate <- mean(tmp$rate)
#}

Rate <- as.data.frame(do.call(rbind,lapply(seq(1,nrow(acc)-1,by=2), function(i){ 
  tmp         <- acc[c(i,i+1)]
  mean(tmp$rate)
})))

Avg$rate <- Rate[,1]
print(head(Avg))

Avg$cell <- gsub("T_", "TD", Avg$cell)
#Avg$condition <- sub("_([^_]*)$", '', Avg$cell)
Avg$condition <- sub("_.*", "", Avg$cell)
###################

Avg_list <- split(Avg, Avg$condition)
Avg_df <- data.frame()

for (elem in Avg_list) {
    cell_Avg <- tapply(elem$rate, elem$dist, mean)
    cell_sd <- tapply(elem$rate, elem$dist, sd)

    test_df <- as.data.frame(cbind(cell_Avg, cell_sd))
    test_df$dist <- rownames(test_df)
    test_df <- as.data.table(test_df)
    test_df$dist <- as.integer(test_df$dist)
    test_df$condition <- rep(unique(elem$condition),length(test_df$dist))

    Avg_df <- rbind(Avg_df, test_df)
}

print("plot average at promoters")
p <- ggplot(Avg_df, aes(dist)) +
    geom_ribbon(aes(y = cell_Avg, ymin = cell_Avg-cell_sd, ymax = cell_Avg+cell_sd), alpha = 0.25) +
    geom_line(aes(y = cell_Avg), color = "red") +
    facet_wrap(~condition, ncol=2) +
    theme_bw() +
    guides(colour = FALSE) +
    labs(y="Rate",x="Distance from Promoter")

save_plot(paste(io$out_dir, paste0(opts$outName, "_average_promoters.pdf"),sep="/"), p)
dev.off()

print("plot all cells at promoters")
p <- ggplot(Avg, aes(dist)) +
    #geom_line(aes(y = sd_down, colour = cell), linetype = "dashed") +
    geom_line(aes(y = rate, colour = cell)) +
    #geom_line(aes(y = sd_up, colour = cell), linetype = "dashed") +
    facet_wrap(~condition, ncol=2) +
    theme_bw() +
    guides(colour = FALSE) +
    labs(y="Rate",x="Distance from Promoter")

save_plot(paste(io$out_dir, paste0(opts$outName, "_at_promoters.pdf"),sep="/"), p)
dev.off()
print("done")
### repeat with only cells passing stricter threshold ###

#if(!is.na(opts$covCutOff)) {
#    if(opts$covCutOff > 0){
#        qc <- fread(io$qc_file) %>%
#            .[, sample := strsplit(sample, "_") %>% map_chr(1)]
#        top_cells <- qc[context ==  opts$context & coverage > opts$covCutOff, id]
#        sub <-  Avg[cell %in% top_cells]
#        p2 <- ggplot(sub, aes(dist, rate, colour = cell)) +
#            theme_bw() +
#            geom_line() 
#        p2
#        pdf(paste(io$out_dir, paste0(opts$outName, "_at_promoters_top_cells.pdf"),sep="/"), 8, 8)
#        print(p2)
#        dev.off()
#    }
#}


