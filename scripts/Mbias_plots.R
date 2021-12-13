args <- commandArgs()

help <- function(){
    cat("Mbias_plots.R :
- convert Mbias txt files to pdf plots
- output is pdf of Mbias plot")
    cat("Usage: /n")
    cat("--outdir    : Path to ouput dir                        [ required ]\n")
    cat("--infile    : Mbias txt file                           [ required ]\n")
    cat("--sample    : name of sample                           [ required ]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(!is.na(charmatch("--help",args)) || !is.na(charmatch("-h",args)) ){
    help()
} else{
    outdir   <- sub('--outdir=', '', args[grep('--outdir=', args)] )
    infile   <- sub('--infile=', '', args[grep('--infile=', args)] )
    sample_name  <- sub('--sample=', '', args[grep('--sample=', args)] )
}


library(ggplot2)

#infile <- "mega_bismarkSE/CX/FB2-4A_B04_S137_merged.M-bias.txt"
#sample_name <- "FB2-4A_B04_S137_merged"
#outdir <- "plots/mega_Mbias"

if(!(file.exists( outdir ))) {
    dir.create(outdir,FALSE,TRUE)
}

Mbias <- read.table(infile, sep = "\t", fill = T, stringsAsFactors=F)

#cleaned_table <- Mbias[-c(1,2,3,97,98,99,193,194,195),]
cleaned_table <- Mbias[-c(grep("context", Mbias$V1), grep("==", Mbias$V1), grep("position", Mbias$V1)),]
colnames(cleaned_table) <- c("position", "count methylated", "count unmethylated", "% methylated", "coverage")
rep_num <- max(as.numeric(cleaned_table$position))
cleaned_table$context <- c(rep("CpG", rep_num), rep("CHG", rep_num), rep("CHH", rep_num))

cleaned_table$position <- as.numeric(cleaned_table$position)
cleaned_table$`% methylated` <- as.numeric(cleaned_table$`% methylated`)

pdf(paste0(outdir, sample_name, "_Mbias_plot.pdf"))
ggplot(data=cleaned_table, mapping=(aes(x=position,y=`% methylated`, color = context))) + geom_line()
dev.off()
