args <- commandArgs()

help <- function(){
    cat("getEnsemblGenes.R :\n")
    cat("Usage: \n")
    cat("--gtf: path to gtf                   [required] \n")
    cat("--met_prom: distance upstream of tss [required] \n")
    cat("-acc_prom: distance upstream of tss  [required] \n")    
     cat("\n")
    q()
}

io   <- list()
opts <- list()

## Save values of each argument
if( !is.na(charmatch("-h",args)) || !is.na(charmatch("--help",args)) ){
    help()
} else {
    io$outfile <- sub( '--gtf=', '', args[grep('--gtf=', args)] )
    opts$met_prom <- as.numeric(sub( '--met_prom=', '', args[grep('--met_prom=', args)] ))
    opts$acc_prom <- as.numeric(sub( '--acc_prom=', '', args[grep('--acc_prom=', args)] ))
}


library(data.table)
library(purrr)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(GenomeInfoDb)

############
#io$outfile <- "/home/groups/CEDAR/anno/gtf/hg19_ens87.chr.gtf"
#opts$acc_prom <- 1000
#opts$met_prom <- 10000

## also generate gene metadata ##
chrs <- sub("^", "chr", c("X", "Y", "MT", 1:22))

print(opts)

if(!(file.exists( "data/anno" ))) {
    print(paste("mkdir", "data/anno"))
    dir.create("data/anno",FALSE,TRUE)  
}

#chrs <- c("X", "Y", "MT", 1:22)

gtf <- rtracklayer::import(io$outfile)


gtf <- as.data.frame(gtf) %>%
  setDT() %>%
  .[type == "gene" & gene_biotype == "protein_coding" & seqnames %in% chrs, 
    .(ens_id = gene_id,
      gene = gene_name, 
      chr = seqnames, 
      start, 
      end, 
      strand)] 

#gtf$chr <- sub("chr", "", gtf$chr)
print("write data/gene_metadata.tsv")
write.table(gtf, file="data/gene_metadata.tsv", sep = "\t", row.names=FALSE, col.names=TRUE)

#######################
## save regions for acc
#######################
prom <- gtf %>%
  .[strand == "+", tss := start] %>%
  .[strand == "-", tss := end] %>% .[!chr %in% c("chrM", "chrY", "M", "Y")]
print(prom)

                                        #gr <- as(prom, "GRanges") %>% .[width(.) > opts$acc_prom]
gr <- as(prom, "GRanges")
gr <- gr[width(gr) > opts$acc_prom]
print(gr)

# save gene body info for correlations
body      <- resize(gr, fix='end', width=width(gr)-opts$acc_prom )
body$anno <- "body"
body      <- as.data.table(body)
df        <- body[,.(seqnames, start, end, strand, ens_id, anno)]
print(head(df))

print("save body without the promoter window")
fwrite(df
     , paste0("data/anno", "body", opts$acc_prom, ".bed"), sep="\t", col.names = F, row.names = F, quote=F, na="NA")

# save promoter info
prom <- as.data.table(promoters(as(prom, "GRanges"), upstream=opts$acc_prom, downstream=opts$acc_prom))

df   <- prom[,.(seqnames, start, end, strand, ens_id, "promoter")]

df <- df[!seqnames %in% c("MT", "Y"), ]

print("save promoter window")
fwrite(df
     , paste0("data/anno/", "promoters", opts$acc_prom, ".bed"), sep="\t", col.names = F, row.names = F, quote=F, na="NA")

#######################
## save regions for met
#######################
prom <- gtf %>%
  .[strand == "+", tss := start] %>%
  .[strand == "-", tss := end] %>% .[!chr %in% c("MT", "Y")]
print(prom)

gr <- as(prom, "GRanges") %>% .[width(.) > opts$met_prom]
print(gr)

# save gene body info for correlations
body      <- resize(gr, fix='end', width=width(gr)-opts$met_prom )
body$anno <- "body"
body      <- as.data.table(body)
df        <- body[,.(seqnames, start, end, strand, ens_id, anno)]
print(head(df))

print("save body without the promoter window")
fwrite(df
     , paste0("data/anno/", "body", opts$met_prom, ".bed"), sep="\t", col.names = F, row.names = F, quote=F, na="NA")

# save promoter info
prom <- as.data.table(promoters(as(prom, "GRanges"), upstream=opts$acc_prom, downstream=opts$met_prom))

df   <- prom[,.(seqnames, start, end, strand, ens_id, "promoter")]

df <- df[!seqnames %in% c("MT", "Y"), ]

print("save promoter window")
fwrite(df
     , paste0("data/anno/", "promoters", opts$met_prom, ".bed"), sep="\t", col.names = F, row.names = F, quote=F, na="NA")
