args <- commandArgs()

help <- function() {
    cat("plotRatesByAnno.R :
Plots rates across the different annotations\n")
    cat("Usage: \n")
    cat("--meta        : stats table with QC info     [required]\n")
    cat("--met         : path to methylation files    [required]\n")
    cat("--acc         : path to accessibility files  [required]\n")
    cat("--plotdir     : location to output plots     [required]\n")
    cat("--anno        : path to annotation files     [required]\n")
    cat("\n")
    q()
}

io <- list()

## Save values of each argument
if( !is.na(charmatch("--help",args)) || !is.na(charmatch("--help",args)) ){
    help()
} else {
    io$meta_data  <- sub( '--meta=', '', args[grep('--meta', args)] )
    io$met_dir    <- sub( '--met=', '', args[grep('--met', args)] )
    io$acc_dir    <- sub( '--acc=', '', args[grep('--acc', args)] )
    io$plot_dir   <- sub( '--plotdir=', '', args[grep('--plotdir', args)] )
    io$anno_dir   <- sub( '--anno=', '', args[grep('--anno', args)] )
}

if ("weights" %in% rownames(installed.packages()) == FALSE) {
    install.packages("weights", repos="https://ftp.osuosl.org/pub/cran/")
}

library(scater)
library(data.table)
library(purrr)
library(weights)
library(ggplot2)
library(cowplot)
library(ggrepel)



#io <- list()
#io$meta_data <- "tables/sample_stats_qcPass.txt"
#io$met_dir   <- "data/met"
#io$acc_dir   <- "data/acc"
#io$plot_dir  <- "plots"
#io$anno_dir  <- "data/anno"

opts <- list()
opts$anno_regex <- "CGI_promoter|MCF7_ER_peaks|H3K27ac_peaks|body|Repressed|Enhancer|CTCF"


meta <- fread(io$meta_data) %>%
  .[pass_accQC == TRUE & pass_metQC == TRUE]

met <- dir(io$met_dir, pattern = ".tsv.gz$", full = TRUE) %>%
  .[grep(opts$anno_regex, .)] %>%
  # map(fread_gz) %>%
  map(fread) %>%
  rbindlist() %>% .[,Type:="met"]

acc <- dir(io$acc_dir, pattern = ".tsv.gz$", full = TRUE) %>%
  .[grep(opts$anno_regex, .)] %>%
  # map(fread_gz) %>%
  map(fread) %>%
  rbindlist() %>% .[,Type:="acc"]

merg <- rbind(acc, met)

boo <- merg[!is.na(var)]

if( !exists(io$plot_dir) ){
  dir.create(io$plot_dir)
}

##############################
## plot varience across cells
# https://github.com/PMBio/scNMT-seq/blob/master/EB_cells/stats/stats_accmet_features_parsed.Rmd

boxplot_theme <- function() {
  p <- theme(
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(colour="black", size=20, vjust=1.5, margin=margin(10,0,0,0)),
    axis.text.x = element_text(colour="black",size=rel(1.7)),
    axis.text.y = element_text(colour="black",size=rel(1.7)),
    axis.line = element_line(colour="black", size=rel(0.7)),
    axis.ticks.x = element_line(colour="black", size=rel(1.0)),
    axis.ticks.y = element_line(colour="black", size=rel(1.0)),
    # axis.ticks.y = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position="top",
    legend.text=element_text(size=15),
    legend.title=element_blank(),
    legend.background=element_blank(),
    panel.border = element_blank()
  )
}

feature_stats <- merg[,.(mean=mean(rate, na.rm=T), var=var(rate, na.rm=T)), by=c("anno","id","Type")]

p1 <- ggplot(feature_stats[Type=="acc"], aes(x=anno, y=mean)) +
  geom_boxplot(aes(fill=Type), alpha = 1, coef=1.5, outlier.shape=NA, color = "gray50") +
  #facet_grid(Type ~ .) + 
  ggtitle("") + xlab("") + ylab("Mean accessibility rate") +
  scale_fill_manual(values=c("#F87D42")) +
  coord_flip() +
  boxplot_theme()
p1

p2 <- ggplot(feature_stats[Type=="met"], aes(x=anno, y=mean)) +
  geom_boxplot(aes(fill=Type), alpha = 1, coef=1.5, outlier.shape=NA, color="gray50") +
  #facet_grid(.~Type ~ .) + 
  ggtitle("") + xlab("") + ylab("Mean accessibility rate") +
  scale_fill_manual(values=c("#00136C")) +
  coord_flip() +
  boxplot_theme()
p2

p3 <- ggplot(feature_stats[Type=="acc"], aes(x=anno, y=var)) +
  geom_boxplot(aes(fill=Type), alpha = 1, coef=1.5, outlier.shape=NA, color = "gray50") +
  #facet_grid(Type ~ .) + 
  ggtitle("") + xlab("") + ylab("Cell-to-cell variance on the accessibilty rate") +
  scale_fill_manual(values=c("#F87D42")) +
  coord_flip() +
  boxplot_theme()
p3

p4 <- ggplot(feature_stats[Type=="met"], aes(x=anno, y=var)) +
  geom_boxplot(aes(fill=Type), alpha = 1, coef=1.5, outlier.shape=NA, color="gray50") +
  #facet_grid(.~Type ~ .) + 
  ggtitle("") + xlab("") + ylab("Cell-to-cell variance on the methylation rate") +
  scale_fill_manual(values=c("#00136C")) +
  coord_flip() +
  boxplot_theme()
p4

pdf(paste0(io$plot_dir,"/anno_rateVarboxplots.pdf"), height=12, width=17)
print(cowplot::plot_grid(p1,p2,p3,p4, labels = c("Mean Rate","Mean Rate","Varience","Varience"), label_size=20, ncol=2, nrow=2))
dev.off()
