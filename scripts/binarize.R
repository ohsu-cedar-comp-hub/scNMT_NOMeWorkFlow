args <- commandArgs()

help <- function(){
    cat("binarize.R :
- (Optional but recommended) Binarise CpG sites and calculate rate.
- output is <chr> <pos> <met_reads> <nonmet_reads> <rate>\n")
    cat("Usage: \n")
    cat("--infile      : Path to input file (.gz)   [ required ]\n")
    cat("--outdir      : Path to output dir         [ required ]\n")
    cat("--input_format: infile format              [ required ]
                           1. chr, position, rate (Bismark default)
                           2. chr, position, met_reads non nomet_reads \n")
    cat("\n")
    q()
}

## Save values of each argument
if(!is.na(charmatch("--help",args)) || !is.na(charmatch("-h",args)) ){
    help()
} else {
    infile       <- sub( '--infile=', '',  args[grep('--infile=', args)] )
    outdir       <- sub( '--outdir=', '',  args[grep('--outdir=',args)] )
    input_format <- as.numeric( sub('--input_format=', '',args[grep('--input_format=',args)]) )

}

#####################################
## Script to filter specific sites ##
#####################################

# Current filters:
# - dinucleotides (non-cg methylation for example)
# - chromosomes

## Input: 
# single-cell methylation files output from Bismark. In either one of the two following formats:
# input_format=1:
# chr     pos     rate
# 1       3019021 0
# 1       3027398 100
# 1       3052955 100

# input_format=2:
# chr     pos     met_reads non nomet_reads
# 1       3019021 0 1
# 1       3027398 1 1
# 1       3052955 1 0

## Output:
# Same format as the input

## Load libraries
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(argparse))

## Define options
opts <- list()

# Define options
opts$input_format <- 2
opts$remove50 <- TRUE # if TRUE, sites with methylation rate of 50% are removed, otherwise they are rounded to 100%
print(opts)

if(!(file.exists( outdir ))) {
    dir.create(outdir,FALSE,TRUE)  
}

cat("Options:\n")
cat(sprintf('- Remove sites with 50 percent methylation rates: %s\n',paste(opts$remove50)))
cat("\n")
    
# Parallelise processing
Sample <- sub(".tsv.gz", "", basename(infile))
outfile <- paste(outdir,Sample,sep="/")
print(outfile)
# Load data
cat(sprintf("Processing %s...\n",Sample))
data <- fread(sprintf("zcat < %s",infile), verbose=F, showProgress=F)

if (opts$input_format == 1) {
    colnames(data) <- c("chr","pos","rate")
    ## data[,rate:=round(rate*100)] # if rate goes from 0 to 1...
} else if (opts$input_format == 2) {
    colnames(data) <- c("chr","pos","met_reads","nonmet_reads")
    data[,rate:=round((met_reads/(met_reads+nonmet_reads))*100)]
}
    
# Sanity check
tmp <- sum((max(data$rate) > 100) | (min(data$rate) < 0))

if (tmp>0) {
    cat(sprintf("%s: There are %d CpG sites that have methylation rate higher than 100 or lower than 0\n",Sample,tmp))
}

# Deal with uncertain sites with rate=50
if (opts$remove50) {
    data <- data[rate!=50,]
} else {
    data[rate==50,rate:=100]
}
    
# Calculate binary methylation status
cat(sprintf("%s: There are %0.03f%% of sites with non-binary methylation rate\n"
          , Sample, 100*mean(!data$rate %in% c(0,100))
            ))
data[,rate:=round(rate/100)]

# Save results
fwrite(data, file=outfile, sep="\t", showProgress=FALSE, verbose=FALSE, col.names=TRUE)
system(sprintf("gzip -f %s",outfile))

print(paste("wrote file:", outfile))
