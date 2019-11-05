args <- commandArgs()

help <- function(){
    cat("counts_stats.R :
- Count Stats.
- output is plots found in plots folder")
    cat("Usage: \n")
    cat("--outdir    : Path to output dir  [ required ]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(!is.na(charmatch("--help",args)) || !is.na(charmatch("-h",args)) ){
    help()
} else {
    outdir   <- sub( '--outdir=', '', args[grep('--outdir=', args)] )

}

## Load libraries
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(argparse))

if(!(file.exists( outdir ))) {
    dir.create(outdir,FALSE,TRUE)  
}

Dir <- outdir

cellsToDrop <- "Ignore"

io <- list()
opts <- list()
io$in.metadata <- "data/metadata.txt"
io$in.acc_data <- "bismarkSE/CX/coverage2cytosine_1based/filt/binarised"
io$in.met_data <- "bismarkSE/CX/coverage2cytosine_1based/filt/binarised"
io$file <- "tables/sample_stats.txt"
io$reads <- "tables/sample_read_report.txt"


# Load metadata
metadata <- fread(io$in.metadata)
head(metadata)
#metadata$V8 <- NULL
names(metadata) <- sub("SampleID", "sample", names(metadata))

opts$cells <- metadata[,Cell]

stats <- rbindlist( 
                   list( 
                       data.table(sample=opts$cells, context="CG", coverage=0, mean=0),  
                       data.table(sample=opts$cells, context="GC", coverage=0, rate=0) 
                   ))
head(stats)

for (cell in opts$cells) {
  # Met
    if ( file.exists(sprintf("%s/%s_CpG.gz",io$in.met_data,cell)) ) {
     print( sprintf("Loading %s methylation...",cell) )
     tmp <- fread(sprintf("zcat < %s/%s_CpG.gz",io$in.met_data,cell)
                , sep="\t", verbose=F, showProgress=F)
     stats[sample==cell & context=="CG",coverage:=nrow(tmp)] 
     stats[sample==cell & context=="CG",mean:=mean(tmp$rate)] 
   } else { 
     print(sprintf("Sample %s not found for methylation",cell)) 
   }
  # Acc
    if ( file.exists(sprintf("%s/%s_GpC.gz", io$in.met_data, cell )) ) {
     print(sprintf("Loading %s accessibility...",cell))
     tmp <- fread(sprintf("zcat < %s/%s_GpC.gz",io$in.met_data,cell)
                , sep="\t", verbose=F, showProgress=F)
     stats[sample==cell & context=="GC",coverage:=nrow(tmp)] 
     stats[sample==cell & context=="GC",mean:=mean(tmp$rate)] 
   } else { 
       print(sprintf("Sample %s not found for accessibility",cell))
   }
}
print(head(stats))
print(dim(stats))

fwrite(stats,file=io$file, sep="\t", row.names=F, col.names=T)



#stats <-
#    fread(io$file) %>% .[sample %in% metadata[,sample]] %>% 
#    merge(metadata,by="sample") %>% head()

############################################
## Number of observed CG/GC sites per sample
############################################
stats <- fread(io$file)

opts$cells_to_drop <- NULL#stats$sample[grep(cellsToDrop, stats$sample)]

if(length(opts$cells_to_drop) >0){
    print("filter cells by pattern")
    stats <- stats[!stats$sample %in% opts$cells_to_drop,]
}
print(head(stats))
print(dim(stats))

opts$met_coverage_threshold <- 1e5
opts$acc_coverage_threshold <- 1e6
foo <- stats[,c("sample","coverage","context")]
#foo <- foo[with(foo, order(coverage)), ]

foo$cellcolor  <- rep("black", nrow(foo))
foo[foo$context=="GC" & foo$coverage<opts$acc_coverage_threshold, "cellcolor"]  <- "red"
foo[foo$context=="CG" & foo$coverage<opts$met_coverage_threshold, "cellcolor"]  <- "red"

cellcolor <- foo$cellcolor

p <- ggplot(foo, aes(x=sample, y=coverage, fill=context)) +
    geom_bar(stat='identity', position="dodge") +
    scale_fill_manual("legend", values = c("GC" = "#F87D42", "CG" = "#00136C")
                    , labels=c("GC" = "GC accessibility","CG" = "CG methylation"))+
    scale_y_continuous(expand=c(0,0)) +
    geom_hline(yintercept=opts$met_coverage_threshold, colour="#00136C", linetype="dashed") +
    geom_hline(yintercept=opts$acc_coverage_threshold, colour="#F87D42", linetype="dashed") +
    ylab("Number of observed sites") +
    xlab("Threshold 1e5")+
  theme(
    axis.title.y = element_text(colour="black", size=17, margin=margin(0,20,0,0)),
    axis.title.x = element_text(colour="black", size=17, margin=margin(0,20,0,0)),
    axis.text.x = element_text(colour="black", color=cellcolor, angle=90, size=10, vjust=0.5, hjust=1.0),
    axis.text.y = element_text(colour="black", size=11),
    axis.ticks = element_line(colour="black"),
    axis.line = element_blank(),
    legend.position="top",
    legend.title = element_blank(),
    legend.direction = "horizontal",
    legend.key.width=unit(1.2,"line"),
    legend.key.height=unit(1.0,"line"),
    legend.text = element_text(size=15),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank()
    )
p

pdf(paste(Dir, "coverageBar.pdf", sep="/"), height=5, width=12)
print(p)
dev.off()

#############
# Methylation mean rate bar
#############
stats <- fread("tables/sample_stats.txt")


# foo <- stats[context=="CG" & pass_metQC==T,c("id","mean")]
foo <- stats[context=="CG",c("sample","mean")]
foo %>% setkey(mean) %>% .[,sample:=factor(sample,levels=sample)]

p <- ggplot(foo, aes(x=sample, y=mean)) +
  geom_bar(stat='identity', position="dodge", fill="#00136C") +
  scale_y_continuous(expand=c(0,0)) +
    ylab("Genome-wide mean methylation rate") +
    ylim(0,1)+
  theme(
    axis.title.y = element_text(colour="black", size=17, margin=margin(0,20,0,0)),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour="black",size=rel(1.5)),
    axis.ticks = element_line(colour="black"),
    axis.line = element_line(color="black"),
    legend.position="top",
    legend.title = element_blank(),
    legend.direction = "horizontal",
    legend.key.width=unit(1.2,"line"),
    legend.key.height=unit(1.0,"line"),
    legend.text = element_text(size=15),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )
p

pdf(paste(Dir, "GWmethylRate.pdf", sep="/"), height=5, width=7)
print(p)
dev.off()

##############
# Acessibility mean rate bar
##############
foo <- stats[context=="GC",c("sample","mean")]
foo %>% setkey(mean) %>% .[,sample:=factor(sample,levels=sample)]

p <- ggplot(foo, aes(x=sample, y=mean)) +
    geom_bar(stat='identity', position="dodge", fill="#F87D42") +
    scale_y_continuous(expand=c(0,0)) +
    ylab("Genome-wide mean accessibility rate") + 
    ylim(0,1)+
    theme(
        axis.title.y = element_text(colour="black", size=17, margin=margin(0,20,0,0)),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=rel(1.5)),
        axis.ticks = element_line(colour="black"),
        axis.line = element_line(color="black"),
        legend.position="top",
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.key.width=unit(1.2,"line"),
        legend.key.height=unit(1.0,"line"),
        legend.text = element_text(size=15),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank()
    )
p

pdf(paste(Dir, "GWaccessRate.pdf", sep="/"), height=5, width=7)
print(p)
dev.off()

##########################
# make read mapping report
##########################
Reads       <- read.table("tables/bismarkSE_mapping_report.txt", sep="\t", header = T, fill = T)
sampleReads <- data.table(sample=paste(stats$sample) )

for( x in unique(paste(stats$sample)) ){    
    Tot <- sum(Reads[Reads$Sample==paste0(x,"_R1") | Reads$Sample==paste0(x,"_R2") , which(colnames(Reads)=="TotalReads")] )
    Map <- sum(Reads[Reads$Sample==paste0(x,"_R1") | Reads$Sample==paste0(x,"_R2") , which(colnames(Reads)=="TotalReads")] )
    Per <- round(Map/Tot*100,2)
    sampleReads[sample==x ,total_reads:=Tot]
    sampleReads[sample==x , mapped_reads:=Map]
    sampleReads[sample==x , percent_mapped:=Per]
}

fwrite(sampleReads,file=io$reads, sep="\t", row.names=F, col.names=T)

############################################################
## Mean methylation levels versus Mean accessibility levels
############################################################
foo <- stats[,c("sample","context","mean")] %>% dcast(sample~context, value.var="mean")

p <- ggplot(foo, aes(x=CG, y=GC)) +
  geom_point() +
  labs(x="Mean methylation rate", y="Mean accessibility rate") +
  theme_bw() + theme(
    axis.title.y = element_text(colour="black", size=17, margin=margin(0,10,0,0)),
    axis.title.x = element_text(colour="black", size=17, margin=margin(10,0,0,0)),
    axis.text.x = element_text(colour="black", size=rel(1.5)),
    axis.text.y = element_text(colour="black", size=rel(1.5))
  )
p

pdf(paste(Dir, "meanMethyl_vs_meanAccess.pdf", sep="/"), height=5, width=5)
print(p)
dev.off()

########################
## Mean rate vs coverage
########################

# Methylation
foo <- stats[context=="CG",c("sample","context","mean","coverage")] 

p <- ggplot(foo, aes(x=mean, y=coverage)) +
  geom_point() +
  stat_smooth(method="lm") +
  #facet_wrap(~stage, scales="free", nrow = 3) +
    labs(x="Mean methylation rate", y="Coverage") +
  theme_bw() + theme(
    axis.title.y = element_text(colour="black", size=17, margin=margin(0,10,0,0)),
    axis.title.x = element_text(colour="black", size=17, margin=margin(10,0,0,0)),
    axis.text.x = element_text(colour="black", size=rel(1.5)),
    axis.text.y = element_text(colour="black", size=rel(1.5))
  )
p

# Accessibility
foo <- stats[context=="GC",c("sample","context","mean","coverage")] #%>% merge(metadata[,c("id_acc","stage")]%>%setnames("id_acc","id"))

p1 <- ggplot(foo, aes(x=mean, y=coverage)) +
  geom_point() +
  stat_smooth(method="lm") +
  labs(x="Mean accessibility rate", y="Coverage") +
  theme_bw() + theme(
    axis.title.y = element_text(colour="black", size=17, margin=margin(0,10,0,0)),
    axis.title.x = element_text(colour="black", size=17, margin=margin(10,0,0,0)),
    axis.text.x = element_text(colour="black", size=rel(1.5)),
    axis.text.y = element_text(colour="black", size=rel(1.5))
  )
p1

pdf(file=paste(Dir, "meanRatevsCoverage.pdf", sep="/"), width=19, height=7)
print(cowplot::plot_grid(p, p1, ncol=2, nrow=1))
dev.off()

################################
## histograms
################################
foo <- stats[order(mean),]
foo %>% setkey(mean) %>% .[,sample:=factor(sample,levels=unique(sample))]


p <- ggplot(foo, aes(x=sample, y=mean, group=context)) +
  geom_line(aes(linetype=context))+
    geom_point(aes(shape=context))+
    theme_bw() +
    ylab("Mean rates")+
    theme(
    axis.title.y = element_text(colour="black", size=17, margin=margin(0,20,0,0)),
    axis.title.x = element_text(colour="black", size=17, margin=margin(0,20,0,0)),
    axis.text.x = element_text(colour="black", color=cellcolor, angle=90, size=10, vjust=0.5, hjust=1.0),
    axis.text.y = element_text(colour="black", size=11),
    axis.ticks = element_line(colour="black"),
    legend.position="top",
    legend.direction = "horizontal",
    legend.key.width=unit(1.2,"line"),
    legend.key.height=unit(1.0,"line"),
    legend.text = element_text(size=15),
    )
print(p)

pdf(file=paste(Dir, "meanRateLinePlot.pdf", sep="/"), width=7, height=5)
print(p)
dev.off()

p <- ggplot(foo, aes(x=sample, y=coverage, group=context)) +
  geom_line(aes(linetype=context))+
    geom_point(aes(shape=context))+
    theme_bw() +
    ylab("Number of observed sites") +
    theme(
    axis.title.y = element_text(colour="black", size=17, margin=margin(0,20,0,0)),
    axis.title.x = element_text(colour="black", size=17, margin=margin(0,20,0,0)),
    axis.text.x = element_text(colour="black", color=cellcolor, angle=90, size=10, vjust=0.5, hjust=1.0),
    axis.text.y = element_text(colour="black", size=11),
    axis.ticks = element_line(colour="black"),
    legend.position="top",
    legend.direction = "horizontal",
    legend.key.width=unit(1.2,"line"),
    legend.key.height=unit(1.0,"line"),
    legend.text = element_text(size=15),
    )
print(p)

pdf(file=paste(Dir, "coverageLinePlot.pdf", sep="/"), width=7, height=5)
print(p)
dev.off()


p <-
    ggplot(stats, aes(x=context, y=mean)) +
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

pdf(file=paste(Dir, "meanRateBoxPlot.pdf", sep="/"), width=5, height=5)
print(p)
dev.off()

p <-
    ggplot(stats, aes(x=context, y=coverage)) +
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

pdf(file=paste(Dir, "coverageBoxPlot.pdf", sep="/"), width=5, height=5)
print(p)
dev.off()

p <- ggplot(stats)+
        geom_density(aes(x = coverage, color = context))+
        xlab("Number of observed sites")+
        scale_color_brewer(palette="Dark2")+
        scale_fill_brewer(palette="Dark2")+
        theme_bw()+
        #xlim(0,3000)+
        theme(text = element_text(size=12)#,
              )
print(p)

pdf(file=paste(Dir, "coverageDensityPlot.pdf", sep="/"), width=5, height=5)
print(p)
dev.off()

## coverage distributions
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

for (Context in c("GC","CG")){
    print(Context)
    stats.tmp <- stats[context==Context,]
    names(stats.tmp)  <- sub("sample", "id", names(stats.tmp))
    if(Context=="GC"){
        Title <- "Chromatin accessibility"
        fname <- "qc_accessCoverageBar.pdf"
        col   <- "#F87D42"
        tmp   <- stats.tmp[,c("id","coverage")] %>%
            setkey(coverage) %>% .[,id:=factor(id,levels=id)]
        tmp$cellcolor <- c("black","red")[as.numeric(tmp$coverage < opts$acc_coverage_threshold)+1]
        Threshold <- opts$acc_coverage_threshold
    }else{
        Title <- "DNA methylation"
        fname <- "qc_methylCoverageBar.pdf"
        col   <- "#00136C"
        tmp   <- stats.tmp[,c("id","coverage")] %>%
            setkey(coverage) %>% .[,id:=factor(id,levels=id)]
        tmp$cellcolor <- c("black","red")[as.numeric(tmp$coverage < opts$met_coverage_threshold)+1]
        Threshold     <- opts$met_coverage_threshold
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
    print("Percent failed QC for methylation:")
    failqc <- stats.tmp[coverage<opts$met_coverage_threshold,id]
    print( paste(round(length(failqc)/nrow(stats.tmp)*100, 1), "%") )
    rm(stats.tmp)
    rm(tmp)
}
