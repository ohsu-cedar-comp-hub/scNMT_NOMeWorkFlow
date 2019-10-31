args <- commandArgs()

help <- function(){
    cat("getEnsemblGenes.R :\n")
    cat("Usage: \n")
    cat("--gtf: path to gtf [required] \n")
     cat("\n")
    q()
}

io   <- list()
opts <- list()

## Save values of each argument
if( !is.na(charmatch("--help",args)) || !is.na(charmatch("--help",args)) ){
    help()
} else {
    outfile <- sub( '--gtf=', '', args[grep('--gtf=', args)] )
}


library(data.table)
library(purrr)
library(rtracklayer)
############
#outfile <- "/home/groups/CEDAR/anno/gtf/hg19_ens87.chr.gtf"

## also generate gene metadata ##
chrs <- sub("^", "chr", c("X", "Y", "MT", 1:22))

#chrs <- c("X", "Y", "MT", 1:22)

gtf <- rtracklayer::import(outfile)


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

write.table(gtf, file="data/gene_metadata.tsv", sep = "\t", row.names=FALSE, col.names=TRUE)
