rule get_genes:
     input: rules.create_reports.output
     output:
        "data/gene_hg19.cellRanger_metadata.tsv"
     shell:
        "Rscript scripts/getEnsemblGenes.R"
rule accessibility_profile:
     input:
        "tables/sample_stats_qcPass.txt",
        "data/gene_hg19.cellRanger_metadata.tsv"
     output:
        "data/accessibility_at_promoters.pdf",
        "data/accessibility_at_promoters.RData",
        "data/anno/body.bed",
        "data/anno/promoters.bed"
     shell:
        "Rscript scripts/accessibility_profiles.R --covPath=bismarkSE/CX/coverage2cytosine_1based/filt/binarised --outDir=data --annoFile={input[1]} --qcFile={input[0]} --context=GC --promUp=10000"

rule methylation_profile:
     input:
        "tables/sample_stats_qcPass.txt",
        "data/gene_hg19.cellRanger_metadata.tsv"
     output:
        "data/methylation_at_promoters.pdf",
        "data/methylation_at_promoters.RData",
     shell:
        "Rscript scripts/accessibility_profiles.R --covPath=bismarkSE/CX/coverage2cytosine_1based/filt/binarised --outDir=data --annoFile={input[1]} --qcFile={input[0]} --context=CG --promUp=10000"

rule annotate_acc:
     input:
        "tables/sample_stats_qcPass.txt"
     output:
        "data/acc/body.tsv.gz",
        "data/acc/CGI_promoter.tsv.gz",
        "data/acc/CTCF.tsv.gz",
        "data/acc/Enhancer.tsv.gz",
        "data/acc/MCF7_ER_peaks.tsv.gz",
        "data/acc/MCF7_H3K27ac_peaks.tsv.gz",
        "data/acc/nonCGI_promoter.tsv.gz",
        "data/acc/Repressed.tsv.gz"
     shell:
        "Rscript scripts/accmet/annotate_arw_acc.R --anno=data/anno --raw=bismarkSE/CX/coverage2cytosine_1based/filt/binarised --meta={input} --outdir=data/acc"

rule annotate_met:
     input:
        "tables/sample_stats_qcPass.txt"
     output:
        "data/met/body.tsv.gz",
        "data/met/CGI_promoter.tsv.gz",
        "data/met/CTCF.tsv.gz",
        "data/met/Enhancer.tsv.gz",
        "data/met/MCF7_ER_peaks.tsv.gz",
        "data/met/MCF7_H3K27ac_peaks.tsv.gz",
        "data/met/nonCGI_promoter.tsv.gz",
        "data/met/Repressed.tsv.gz"
     shell:
        "Rscript scripts/accmet/annotate_arw_met.R --anno=data/anno --raw=bismarkSE/CX/coverage2cytosine_1based/filt/binarised --meta={input} --outdir=data/met"

rule plot_anno:
     input:
        "tables/sample_stats_qcPass.txt",
        "data/acc/body.tsv.gz",
        "data/acc/CGI_promoter.tsv.gz",
        "data/acc/CTCF.tsv.gz",
        "data/acc/Enhancer.tsv.gz",
        "data/acc/MCF7_ER_peaks.tsv.gz",
        "data/acc/MCF7_H3K27ac_peaks.tsv.gz",
        "data/acc/nonCGI_promoter.tsv.gz",
        "data/acc/Repressed.tsv.gz",
        "data/met/body.tsv.gz",
        "data/met/CGI_promoter.tsv.gz",
        "data/met/CTCF.tsv.gz",
        "data/met/Enhancer.tsv.gz",
        "data/met/MCF7_ER_peaks.tsv.gz",
        "data/met/MCF7_H3K27ac_peaks.tsv.gz",
        "data/met/nonCGI_promoter.tsv.gz",
        "data/met/Repressed.tsv.gz"
     output:
        "plots/accmet_rates_anno.boxplot.pdf",
        "plots/accmet_variance_anno.boxplot.pdf"
     shell:
        "Rscript scripts/plotRatesByAnno.R --meta={input[0]} --met=data/met --acc=data/acc --plotdir=plots --anno=data/anno"

rule accmet_corr:
     input:
        "tables/sample_stats_qcPass.txt",
	"data/acc/body.tsv.gz",
        "data/acc/CGI_promoter.tsv.gz",
        "data/acc/CTCF.tsv.gz",
        "data/acc/Enhancer.tsv.gz",
        "data/acc/MCF7_ER_peaks.tsv.gz",
        "data/acc/MCF7_H3K27ac_peaks.tsv.gz",
        "data/acc/nonCGI_promoter.tsv.gz",
        "data/acc/Repressed.tsv.gz",
        "data/met/body.tsv.gz",
        "data/met/CGI_promoter.tsv.gz",
        "data/met/CTCF.tsv.gz",
        "data/met/Enhancer.tsv.gz",
        "data/met/MCF7_ER_peaks.tsv.gz",
        "data/met/MCF7_H3K27ac_peaks.tsv.gz",
        "data/met/nonCGI_promoter.tsv.gz",
        "data/met/Repressed.tsv.gz"
     output:
        "plots/corr/acc_met_correlations_loci.pdf",
        "plots/corr/acc_met_correlations_samplesAnno.pdf"
     shell:
        "Rscript scripts/accmet/correlations_accmet.R --meta={input[0]} --met=data/met --acc=data/acc --plotdir=plots/corr --anno=data/anno"
