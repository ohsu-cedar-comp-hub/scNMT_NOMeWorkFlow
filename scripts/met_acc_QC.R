library(data.table)
library(purrr)
library(tidyr)
library(ggplot2)
library(cowplot)

setwd("/home/groups/CEDAR/woodfin/projects/NMT-seq/20190627_NM_rename/")

Dir <- "plots/QC/filtered"
if(!(file.exists(Dir))) {
    dir.create(Dir,FALSE,TRUE)  
}

barplot_theme <- function() {
  p <- theme(
    plot.title = element_text(size=20, hjust=0.5),
    # axis.title.x = element_text(colour="black", size=25, vjust=1.5),
    axis.title.x = element_text(colour="black", size=18),
    axis.title.y = element_text(colour="black", size=18),
    # axis.text.x = element_text(colour="black",size=rel(1.6)),
    axis.text.y = element_text(colour="black",size=rel(1.5)),
    axis.line = element_line(colour="black", size=rel(0.7)),
    axis.ticks.x = element_line(colour="black", size=rel(0.7)),
    axis.ticks.y = element_line(colour="black", size=rel(0.7)),
    #legend.position="none",
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )
}

## set options for filtering
opts <- list()
opts$acc_coverage_threshold    <- 1e6
opts$acc_mean_threshold.bottom <- 20
opts$acc_mean_threshold.top    <- 45
opts$met_coverage_threshold    <- 1e5

## I/O ##
io                         <- list()
io$sample.metadata         <- "tables/sample_read_report.txt"
io$sample.metadata_updated <- "tables/sample_read_report_qcPASS.txt"
io$in.accdir               <- "bismarkSE/CX/coverage2cytosine_1based/filt/binarised"
io$stats                   <- "tables/sample_stats.txt"
io$stats_updated           <- "tables/sample_stats_qcPass.txt"
io$outdir                  <- Dir
print(io)

# load stats and merge with read report
stats <- fread("tables/sample_stats.txt",sep="\t", verbose=F, showProgress=F)
head(stats)
dim(stats)

metadata    <- fread(io$sample.metadata)
metadata$id <- metadata$sample
metadata$condition <- sub("^sc_.*[0-9]_", "", metadata$sample)

df <- data.table(merge(as.data.frame(stats), as.data.frame(metadata[!duplicated(metadata$sample), ]), by="sample"))
dim(df)
head(df)

# Define which cells to use
opts$cells <- metadata$id_acc
opts$cells_to_drop <- stats$sample[grep("^cells100", stats$sample)]

if(length(opts$cells_to_drop) >0){
    print("filter cells by pattern")
    df <- df[!df$sample %in% opts$cells_to_drop,]
}

###############
## function for each context to make ranked coverage and mean rate plots highlighting cells that do not pass QC
#############
for (Context in c("GC","CG")){
    stats.tmp <- df[context==Context,]
    if(Context=="GC"){
        Title <- "Chromatin accessibility"
        fname <- "qc_accessCoverageBar.pdf"
        col   <- "#F87D42"
        tmp   <- stats.tmp[,c("id","coverage", "mean")] %>%
            setkey(coverage) %>% .[,id:=factor(id,levels=id)]
        tmp$cellcolor <- c("black","red")[as.numeric(tmp$coverage < opts$acc_coverage_threshold)+1]
        Threshold <- opts$acc_coverage_threshold
        tmp.mean <- stats.tmp[,c("id","coverage", "mean")] %>%
            setkey(mean) %>% .[,id:=factor(id,levels=id)]       
    }else{
        Title <- "DNA methylation"
        fname <- "qc_methylCoverageBar.pdf"
        col   <- "#00136C"
        tmp   <- stats.tmp[,c("id","coverage", "mean")] %>%
            setkey(coverage) %>% .[,id:=factor(id,levels=id)]
        tmp$cellcolor <- c("black","red")[as.numeric(tmp$coverage < opts$met_coverage_threshold)+1]
        Threshold     <- opts$met_coverage_threshold
        tmp.mean <- stats.tmp[,c("id","coverage", "mean")] %>%
            setkey(mean) %>% .[,id:=factor(id,levels=id)]
    }
    p1 <- ggplot(tmp, aes(x=id, y=coverage)) +
        geom_bar(stat="identity", position="dodge", fill=col, color=col) +
        labs(title=Title, x=paste("cutoff:",Threshold), y="Number of observed CpG sites") +
        #ggtitle(Title)+
        geom_hline(yintercept=Threshold, colour="black", linetype="dashed") +
        barplot_theme() +
        theme(
            axis.text.x = element_text(angle=90, size=9, vjust=0.5, hjust=1.0, color=tmp$cellcolor)
            #axis.text.x = element_blank(),
            #axis.ticks.x = element_blank()
        )   
    print(p1)
    pdf(paste(Dir, fname, sep="/"), width=8, height=5)
    print(p1)
    dev.off()
    p1 <- ggplot(tmp.mean, aes(x=id, y=mean)) +
        geom_bar(stat="identity", position="dodge", fill=col, color=col) +
        labs(title=Title, x=paste("Coverage cutoff:",Threshold), y="Mean Rate") +
        #ggtitle(Title)+
        #geom_hline(yintercept=Threshold, colour="black", linetype="dashed") +
        barplot_theme() +
        theme(
            axis.text.x = element_text(angle=90, size=9, vjust=0.5, hjust=1.0, color=tmp$cellcolor)
            #axis.text.x = element_blank(),
            #axis.ticks.x = element_blank()
        )   
    print(p1)    
    pdf(paste(Dir, sub("Coverage", "MeanCoverage", fname), sep="/"), width=8, height=5)
    print(p1)
    dev.off()   
    print("Percent failed QC for methylation:")
    failqc <- stats.tmp[coverage<opts$met_coverage_threshold,id]
    print( paste(round(length(failqc)/nrow(stats.tmp)*100, 1), "%") )
    rm(stats.tmp)
    rm(tmp)
}

#############
## Boxplots
#############
p <-
    ggplot(df, aes(x=context, y=mean)) +
    geom_boxplot(aes(fill = context), alpha=1.0, outlier.shape = NA) +
    #geom_jitter(alpha=0.5, color=c("#00136C", "#F87D42")) +
    geom_jitter(aes(color = context), alpha=0.5) +
    scale_color_manual(values=c("#00136C", "#F87D42"))+
    scale_fill_manual("legend", values = c("GC" = "#F87D42", "CG" = "#00136C"),
                      labels=c("CG methylation","GC accessibility"))+
    coord_cartesian(ylim=c(0,1)) +
    ylab("Genome-wide mean rates") +
    theme_bw() +
    theme(
        axis.title.y = element_text(colour="black", size=17, margin=margin(0,20,0,0)),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour="black", angle=90, size=10, vjust=0.5, hjust=1.0),
        axis.text.y = element_text(colour="black", size=11),
        axis.ticks = element_line(colour="black"),
        legend.title=element_blank(),
        legend.text = element_text(size=15)
    )
p

pdf(file=paste(Dir, "qc_meanRateBoxPlot.pdf", sep="/"), width=5, height=5)
print(p)
dev.off()

p <-
    ggplot(df, aes(x=context, y=coverage)) +
    geom_boxplot(aes(fill = context), alpha=1.0, outlier.shape = NA) +
    #geom_jitter(alpha=0.5, color=c("#00136C", "#F87D42")) +
    geom_jitter(aes(color = context), alpha=0.5) +
    scale_color_manual(values=c("#00136C", "#F87D42"))+
    scale_fill_manual("legend", values = c("GC" = "#F87D42", "CG" = "#00136C"),
                      labels=c("CG methylation","GC accessibility"))+
    ylab("Number of observed sites") +
    theme_bw() +
    #scale_y_log10()+
    theme(
        axis.title.y = element_text(colour="black", size=17, margin=margin(0,20,0,0)),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour="black", angle=90, size=10, vjust=0.5, hjust=1.0),
        axis.text.y = element_text(colour="black", size=11),
        axis.ticks = element_line(colour="black"),
        legend.title=element_blank(),
        legend.text = element_text(size=15)
    )
p

pdf(file=paste(Dir, "qc_coverageBoxPlot.pdf", sep="/"), width=5, height=5)
print(p)
dev.off()

######################
## table cells that pass QC
######################
## this table has a line for each context with coverage and rate and pass and fail QC
met_failqc <- df[context=="CG" & coverage < opts$met_coverage_threshold, sample]
acc_failqc <- df[context=="GC" & coverage < opts$acc_coverage_threshold, sample]

df <- df[,c("pass_metQC","pass_accQC"):=TRUE]

df[context=="CG" & sample %in% met_failqc,pass_metQC:=FALSE]
df[context=="GC" & sample %in% acc_failqc,pass_accQC:=FALSE]

fwrite(df, io$stats_updated, sep="\t", col.names = T, row.names = F, quote=F, na="NA")

## this table only has 1 line per sample with the read mapping summaries and pass and fail QC


metadata <- metadata[,c("pass_metQC","pass_accQC"):=FALSE]

metadata[!sample %in% met_failqc,pass_metQC:=TRUE]
metadata[!sample %in% acc_failqc,pass_accQC:=TRUE]

metadata[sample %in% opts$cells_to_drop,c("pass_metQC","pass_accQC"):=FALSE]

metadata <- metadata[!duplicated(sample),]


fwrite(metadata, io$sample.metadata_updated, sep="\t", col.names = T, row.names = F, quote=F, na="NA")

#############
## call the number of cpgs that overlap with the annotate script
#############
