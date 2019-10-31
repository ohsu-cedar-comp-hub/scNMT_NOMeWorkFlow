rule get_genes:
     input: rules.create_reports.output
     output:
        "data/gene_metadata.tsv"
     params:
        gtf = config["gtf_file"]
     conda:
        "../envs/GetEnsemblGenes.yaml"
     shell:
        "Rscript scripts/getEnsemblGenes.R --gtf={params.gtf}"

rule accessibility_profile:
     input:
        "tables/sample_stats_qcPass.txt",
        "data/gene_metadata.tsv"
     output:
        expand("data/{regions}{acc_prom}.bed", acc_prom = acc_prom, regions = regions),
        "plots/profiles/accessibility_at_promoters.pdf",
        "plots/profiles/accessibility_at_promoters.RData"
     params:
        promUp = config['acc_prom']
     conda:
        "../envs/NMT_Profiles.yaml"
     shell:
        "Rscript scripts/accessibility_profiles.R --covPath=bismarkSE/CX/coverage2cytosine_1based/filt/binarised --dataDir=data --plotDir=plots/profiles --annoFile={input[1]} --qcFile={input[0]} --context=GC --promUp={params.promUp}"

rule methylation_profile:
     input:
        "tables/sample_stats_qcPass.txt",
        "data/gene_metadata.tsv"
     output:
        expand("data/{regions}{met_prom}.bed", met_prom = met_prom, regions = regions),
        "plots/profiles/methylation_at_promoters.RData",
        "plots/profiles/methylation_at_promoters.pdf"
     params:
        promUp = config['met_prom']	
     conda:
        "../envs/NMT_Profiles.yaml"	 
     shell:
        "Rscript scripts/accessibility_profiles.R --covPath=bismarkSE/CX/coverage2cytosine_1based/filt/binarised --dataDir=data --plotDir=plot/profiles --annoFile={input[1]} --qcFile={input[0]} --context=CG --promUp={params.promUp}"