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