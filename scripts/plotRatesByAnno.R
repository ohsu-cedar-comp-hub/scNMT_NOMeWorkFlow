library(scater)
library(data.table)
library(purrr)
library(weights)
library(ggplot2)
library(cowplot)
library(ggrepel)


io <- list()
io$meta_data <- "Box/external/MCF7_NMTseq_downstream/sample_stats_qcPass.txt"
io$met_dir   <- "Box/external/MCF7_NMTseq_downstream/met/parsed_arw/"
io$acc_dir   <- "Box/external/MCF7_NMTseq_downstream/acc/parsed_arw/"
io$plot_dir  <- "Box/external/MCF7_NMTseq_downstream/plots/corr"
io$anno_dir  <- "Box/external/MCF7_NMTseq_downstream/features/anno/"

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

p <-
  ggplot(boo, aes(x=anno, y=rate, fill=Type)) +
  geom_boxplot(aes(fill = Type), alpha=1.0, outlier.shape = NA) +
  #geom_jitter(alpha=0.5, color=c("#00136C", "#F87D42")) +
  #geom_jitter(aes(color = Type), alpha=0.5) +
  scale_color_manual(values=c("#00136C", "#F87D42"))+
  # scale_fill_manual("legend", values = c("acc" = "#F87D42", "met" = "#00136C"),
  #                  labels=c("CG methylation","GC accessibility"))+
  scale_fill_manual("legend", values = c("acc" = "#F87D42", "met" = "#00136C"))+
  ylab("rate") +
  xlab("") +
  ggtitle("Rates across genomic loci") +
  coord_flip()+
  theme_bw() +
  theme(
    axis.title.y = element_text(colour="black", size=17, margin=margin(0,20,0,0)),
    axis.title.x = element_text(colour="black", size=17, margin=margin(0,20,0,0)),
    axis.text.x = element_text(colour="black", angle=90, size=10, vjust=0.5, hjust=1.0),
    axis.text.y = element_text(colour="black", size=11),
    axis.ticks = element_line(colour="black"),
    legend.title= element_text(size=15),
    legend.text = element_text(size=15)
  )
p

p2 <-
  ggplot(boo, aes(x=anno, y=var, fill=Type)) +
  geom_boxplot(aes(fill = Type), alpha=1.0, outlier.shape = NA) +
  #geom_jitter(alpha=0.5, color=c("#00136C", "#F87D42")) +
  #geom_jitter(aes(color = Type), alpha=0.5) +
  scale_color_manual(values=c("#00136C", "#F87D42"))+
  # scale_fill_manual("legend", values = c("acc" = "#F87D42", "met" = "#00136C"),
  #                  labels=c("CG methylation","GC accessibility"))+
  scale_fill_manual("legend", values = c("acc" = "#F87D42", "met" = "#00136C"))+
  ylab("varience") +
  xlab("") +
  ggtitle("Varience across genomic loci") +
  coord_flip()+
  theme_bw() +
  theme(
    axis.title.y = element_text(colour="black", size=17, margin=margin(0,20,0,0)),
    axis.title.x = element_text(colour="black", size=17, margin=margin(0,20,0,0)),
    axis.text.x = element_text(colour="black", angle=90, size=10, vjust=0.5, hjust=1.0),
    axis.text.y = element_text(colour="black", size=11),
    axis.ticks = element_line(colour="black"),
    legend.title= element_text(size=15),
    legend.text = element_text(size=15)
  )
p2


if( !exists(io$plot_dir) ){
  dir.create(io$plot_dir)
}

save_plot(paste(io$plot_dir, "accmet_rates_anno.boxplot.pdf", sep="/"), p, base_width = 6, base_height = 6)
save_plot(paste(io$plot_dir, "accmet_variance_anno.boxplot.pdf", sep="/"), p2, base_width = 6, base_height = 6)
