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
library(dplyr)

#infile <- "bismarkSE/CX/FB2-4A_D02_S159_merged.M-bias.txt"
#sample_name <- "FB2-4A_D02_S159_merged"
outdir <- "plots/Mbias"

if(!(file.exists( outdir ))) {
    dir.create(outdir,FALSE,TRUE)
}

files <- list.files(path="bismarkSE/CX", pattern = "*bias*", full.name = T)

read_and_clean <- function(file) {
    temp_df <- read.table(file, sep = "\t", fill = T, stringsAsFactors=F)
    cleaned_table <- temp_df[-c(grep("context", temp_df$V1), grep("==", temp_df$V1), grep("position", temp_df$V1)),]
    colnames(cleaned_table) <- c("position", "count methylated", "count unmethylated", "% methylated", "coverage")
    cleaned_table$position <- as.numeric(cleaned_table$position)
    cleaned_table$`count methylated` <- as.numeric(cleaned_table$`count methylated`)
    cleaned_table$`count unmethylated` <- as.numeric(cleaned_table$`count unmethylated`)
    cleaned_table$`% methylated` <- as.numeric(cleaned_table$`% methylated`)
    cleaned_table$coverage <- as.numeric(cleaned_table$coverage)
    cleaned_table$`% methylated`[which(is.na(cleaned_table$`% methylated`))] <- 0
    colnames(cleaned_table) <- NULL
    cleaned_table <- as.matrix(cleaned_table)
    print(dim(cleaned_table))
    return(cleaned_table)
}

test <- lapply(files, read_and_clean)
    
avg_mat <- Reduce("+",test)/length(test)

avg_df <- as.data.frame(avg_mat)

colnames(avg_df) <- c("position", "count methylated", "count unmethylated", "% methylated", "coverage")

avg_df$context <- c(rep("CpG", 143), rep("CHG", 143), rep("CHH", 143))

pdf("plots/Mbias/average_Mbias_plot.pdf")
ggplot(data=avg_df, mapping=(aes(x=position,y=`% methylated`, color = context))) + geom_line()
dev.off()

mat_sd <- function(lst) {
    n <- length(lst)
    rc <- dim(lst[[1]])
    ar1 <- array(unlist(lst), c(rc, n))
    round(apply(ar1, c(1, 2), sd), 2)
}

sd_mat <- mat_sd(test)

avg_df$sd <- sd_mat[,4]

pdf("plots/Mbias/average_Mbias_plot_wSD.pdf")
ggplot(data=avg_df, mapping=(aes(x=position,y=`% methylated`, color = context))) + geom_line() + geom_ribbon(aes(y=`% methylated`, ymin= `% methylated` - sd, ymax = `% methylated` + sd, fill = context), alpha = .2)
dev.off()
