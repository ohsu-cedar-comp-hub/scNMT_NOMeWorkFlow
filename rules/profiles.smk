rule get_genes:
     input: rules.create_reports.output
     output:
        a="data/gene_metadata.tsv",
	b=expand("data/anno/body{acc_prom}.bed", acc_prom = acc_prom),
	c=expand("data/anno/promoters{acc_prom}.bed", acc_prom = acc_prom),
        d=expand("data/anno/body{met_prom}.bed", met_prom = met_prom),
        e=expand("data/anno/promoters{met_prom}.bed", met_prom = met_prom)
     params:
        gtf = config["gtf_file"],
	acc = config['acc_prom'],
	met = config['met_prom']
     conda:
        "../envs/GetEnsemblGenes.yaml"
     shell:
        "Rscript scripts/getEnsemblGenes.R  --gtf={params.gtf} --met_prom={params.met} --acc_prom={params.acc}"

rule accessibility_profile:
     input:
        "tables/sample_stats_qcPass.txt",
        "data/gene_metadata.tsv"
     output:
        "plots/profiles/accessibility_at_promoters.pdf",
        "plots/profiles/accessibility_average_promoters.pdf",	
        "plots/profiles/accessibility_at_promoters.RData"
     params:
        accBed = config['accBed']
     conda:
        "../envs/NMT_Profiles.yaml"
     shell:
        "Rscript scripts/accessibility_profiles.R --covPath=bismarkSE/CX/coverage2cytosine_1based/filt/binarised --dataDir=data --plotDir=plots/profiles --annoFile={input[1]} --qcFile={input[0]} --context=GC --promBed={params.accBed}"

rule methylation_profile:
     input:
        "tables/sample_stats_qcPass.txt",
        "data/gene_metadata.tsv"
     output:
        "plots/profiles/methylation_at_promoters.RData",
        "plots/profiles/methylation_average_promoters.pdf",
        "plots/profiles/methylation_at_promoters.pdf"
     params:
        metBed = config['metBed']
     conda:
        "../envs/NMT_Profiles.yaml"	 
     shell:
        "Rscript scripts/accessibility_profiles.R --covPath=bismarkSE/CX/coverage2cytosine_1based/filt/binarised --dataDir=data --plotDir=plots/profiles --annoFile={input[1]} --qcFile={input[0]} --context=CG --promBed={params.metBed}"

#rule define_cpgProm:
#     input: rules.get_genes.output.c #"data/anno/promoters1000.bed"
#     output:
#        "data/anno/CGI_promoter.bed",
#        "data/anno/nonCGI_promoter.bed"
#     params:
#        bedtools = config['bedtools'],
#        cpg = config['cpg']
#     shell:
#        """
#        {params.bedtools} -u -a {input} -b {params.cpg} > data/anno/CGI_promoter.bed
#        {params.bedtools} -v -a {input} -b {params.cpg} > data/anno/nonCGI_promoter.bed
#        """

rule annotate_acc:
     input:
        "tables/sample_stats_qcPass.txt",
        "data/gene_metadata.tsv"
     output:
        "data/acc/body.tsv.gz",
        "data/acc/CGI_promoter.tsv.gz",
        "data/acc/CTCF.tsv.gz",
        "data/acc/Enhancer.tsv.gz",
        "data/acc/MCF7_ER_peaks.tsv.gz",
        "data/acc/MCF7_H3K27ac_peaks.tsv.gz",
        "data/acc/nonCGI_promoter.tsv.gz",
        "data/acc/Repressed.tsv.gz"
     conda:
        "../envs/NMT_annotate.yaml"        	
     shell:
        "Rscript scripts/accmet/annotate_arw_acc.R --anno=data/anno --raw=bismarkSE/CX/coverage2cytosine_1based/filt/binarised --meta={input[0]} --outdir=data/acc --genes={input[1]}"

rule annotate_met:
     input:
        "tables/sample_stats_qcPass.txt",
	"data/gene_metadata.tsv"
     output:
        "data/met/body.tsv.gz",
        "data/met/CGI_promoter.tsv.gz",
        "data/met/CTCF.tsv.gz",
        "data/met/Enhancer.tsv.gz",
        "data/met/MCF7_ER_peaks.tsv.gz",
        "data/met/MCF7_H3K27ac_peaks.tsv.gz",
        "data/met/nonCGI_promoter.tsv.gz",
        "data/met/Repressed.tsv.gz"
     conda:
        "../envs/NMT_annotate.yaml"        
     shell:
        "Rscript scripts/accmet/annotate_arw_met.R --anno=data/anno --raw=bismarkSE/CX/coverage2cytosine_1based/filt/binarised --meta={input[0]} --outdir=data/met --genes={input[1]}"

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
        "plots/anno_rateVarboxplots.pdf"
     conda:
        "../envs/NMT_plotRatesAnno.yaml"
     shell:
        "Rscript scripts/plotRatesByAnno.R --meta={input[0]} --met=data/met --acc=data/acc --plotdir=plots --anno=data/anno"

rule accmet_corr:
     input:
        "tables/sample_stats_qcPass.txt",
	"data/gene_metadata.tsv",
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
        "plots/corr/acc_met_correlations_loci.tsv"
     params:
        acc = config["min_cells_acc"],
        met = config["min_cells_met"]
     conda:
        "../envs/NMT_plotRatesAnno.yaml"	
     shell:
        "Rscript scripts/accmet/correlations_accmet.R --meta={input[0]} --genes={input[1]} --met=data/met --acc=data/acc --plotdir=plots/corr --anno=data/anno --min_cells_met={params.met} --min_cells_acc={params.acc}"
