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
        "Rscript scripts/getEnsemblGenes.R --gtf={params.gtf} --met_prom={params.met} --acc_prom={params.acc}"

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
        "Rscript scripts/accessibility_profiles.R --covPath=bismarkSE/CX/coverage2cytosine_1based/filt/binarised --dataDir=data --plotDir=plot/profiles --annoFile={input[1]} --qcFile={input[0]} --context=CG --promBed={params.metBed}"

rule define_cpgProm:
     input: rules.get_genes.output.c #"data/anno/promoters1000.bed"
     output:
        "data/anno/CGI_promoter.bed",
        "data/anno/nonCGI_promoter.bed"
     params:
        bedtools = config['bedtools'],
        cpg = config['cpg']
     shell:
        """
        {params.bedtools} -u -a {input} -b {params.cpg} > data/anno/CGI_promoter.bed
        {params.bedtools} -v -a {input} -b {params.cpg} > data/anno/nonCGI_promoter.bed
        """

