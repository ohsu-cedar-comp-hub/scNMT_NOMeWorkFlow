rule trimming:
    input:
        read1 = "samples/raw/{sample}_R1.fastq.gz",
        read2 = "samples/raw/{sample}_R2.fastq.gz"
    output:
        "samples/trim/{sample}_R1.fastq.gz_trimming_report.txt",
        "samples/trim/{sample}_R1_val_1_fastqc.html",
        "samples/trim/{sample}_R1_val_1_fastqc.zip",
        "samples/trim/{sample}_R1_val_1.fq.gz",
        "samples/trim/{sample}_R2.fastq.gz_trimming_report.txt",
        "samples/trim/{sample}_R2_val_2_fastqc.html",
        "samples/trim/{sample}_R2_val_2_fastqc.zip",
        "samples/trim/{sample}_R2_val_2.fq.gz"
    conda:
        "../envs/trimG.yaml"
    shell:
        """trim_galore --paired --clip_R1 6 --clip_R2 6 --trim1 --gzip --fastqc --output_dir samples/trim {input.read1} {input.read2}"""

rule mapping_R1:
    input:
        fwd = "samples/trim/{sample}_R1_val_1.fq.gz",
    output:
        "bismarkSE/{sample}_R1.{sample}_R1_val_1_bismark_bt2.bam",
        "bismarkSE/{sample}_R1.{sample}_R1_val_1_bismark_bt2_SE_report.txt",
    params:	
        bismark_index = config["bismark_index"],
    conda:
        "../envs/methylome.yaml"
    shell:
        """bismark --prefix {wildcards.sample}_R1 --output_dir bismarkSE --non_directional --parallel=2 --score_min=L,-50,-0.2 --gzip -n 1 {params.bismark_index} {input.fwd}"""
        
rule mapping_R2:
    input:
        rev = "samples/trim/{sample}_R2_val_2.fq.gz",
    output:
        "bismarkSE/{sample}_R2.{sample}_R2_val_2_bismark_bt2.bam",
        "bismarkSE/{sample}_R2.{sample}_R2_val_2_bismark_bt2_SE_report.txt",
    params:	
        bismark_index = config["bismark_index"],
    conda:
        "../envs/methylome.yaml"
    shell:
        """bismark --prefix {wildcards.sample}_R2 --output_dir bismarkSE --non_directional --parallel=2 --score_min=L,-50,-0.2 --gzip -n 1 {params.bismark_index} {input.rev}"""

rule deduplcate_bam_R1:
    input:
        fwd = "bismarkSE/{sample}_R1.{sample}_R1_val_1_bismark_bt2.bam",
    output:
        "bismarkSE/dedup/{sample}_R1.{sample}_R1_val_1_bismark_bt2.deduplicated.bam",
	"bismarkSE/dedup/{sample}_R1.{sample}_R1_val_1_bismark_bt2.deduplication_report.txt"
    conda:
        "../envs/methylome.yaml"
    shell:
        """deduplicate_bismark --single --output_dir bismarkSE/dedup --bam {input.fwd}"""

rule deduplcate_bam_R2:
    input:
        rev = "bismarkSE/{sample}_R2.{sample}_R2_val_2_bismark_bt2.bam",
    output:
        "bismarkSE/dedup/{sample}_R2.{sample}_R2_val_2_bismark_bt2.deduplicated.bam",
	"bismarkSE/dedup/{sample}_R2.{sample}_R2_val_2_bismark_bt2.deduplication_report.txt"
    conda:
        "../envs/methylome.yaml"
    shell:
        """deduplicate_bismark --single --output_dir bismarkSE/dedup --bam {input.rev}"""

rule methylation_extractor:
    input:
        fwd = "bismarkSE/dedup/{sample}_R1.{sample}_R1_val_1_bismark_bt2.deduplicated.bam",
        rev = "bismarkSE/dedup/{sample}_R2.{sample}_R2_val_2_bismark_bt2.deduplicated.bam",
    output:
        "bismarkSE/dedup/{sample}_merged.bam",
        "bismarkSE/CX/CHG_CTOT_{sample}_merged.txt.gz",
        "bismarkSE/CX/CHG_OT_{sample}_merged.txt.gz",
        "bismarkSE/CX/CHH_CTOT_{sample}_merged.txt.gz",
        "bismarkSE/CX/CHH_OT_{sample}_merged.txt.gz",
        "bismarkSE/CX/CpG_CTOT_{sample}_merged.txt.gz",
        "bismarkSE/CX/CpG_OT_{sample}_merged.txt.gz",
        "bismarkSE/CX/{sample}_merged.bismark.cov.gz",
        "bismarkSE/CX/{sample}_merged_splitting_report.txt",
        "bismarkSE/CX/CHG_CTOB_{sample}_merged.txt.gz",
        "bismarkSE/CX/CHG_OB_{sample}_merged.txt.gz",
        "bismarkSE/CX/CHH_CTOB_{sample}_merged.txt.gz",
        "bismarkSE/CX/CHH_OB_{sample}_merged.txt.gz",
        "bismarkSE/CX/CpG_CTOB_{sample}_merged.txt.gz",
        "bismarkSE/CX/CpG_OB_{sample}_merged.txt.gz",
        "bismarkSE/CX/{sample}_merged.bedGraph.gz",
        "bismarkSE/CX/{sample}_merged.M-bias.txt"
    conda:
        "../envs/methylome.yaml"
    shell:
        """
	samtools merge bismarkSE/dedup/{wildcards.sample}_merged.bam {input.fwd} {input.rev}
	bismark_methylation_extractor -s --output bismarkSE/CX --CX --parallel 12 --gzip --bedGraph bismarkSE/dedup/{wildcards.sample}_merged.bam
	"""

rule coverage2cytosine:
    input:
        "bismarkSE/CX/{sample}_merged.bismark.cov.gz" 
    output:
        "bismarkSE/CX/coverage2cytosine_1based/{sample}_merged.NOMe.CpG.cov.gz",
        "bismarkSE/CX/coverage2cytosine_1based/{sample}_merged.NOMe.CpG_report.txt.gz",
        "bismarkSE/CX/coverage2cytosine_1based/{sample}_merged.NOMe.GpC.cov.gz",
        "bismarkSE/CX/coverage2cytosine_1based/{sample}_merged.NOMe.GpC_report.txt.gz"	
    params:
        bismark_index = config["bismark_index"]
    conda:
        "../envs/methylome.yaml"
    shell:
        """
        coverage2cytosine --nome-seq --gzip --output {wildcards.sample}_merged --dir bismarkSE/CX/coverage2cytosine_1based --genome_folder {params.bismark_index} {input}
        """

#rule c2c_cov_filter:
#     input:
#        "bismarkSE/CX/coverage2cytosine_1based/{sample}_merged.NOMe.CpG.cov.gz",
#        "bismarkSE/CX/coverage2cytosine_1based/{sample}_merged.NOMe.GpC.cov.gz"
#     output:
#        "bismarkSE/CX/coverage2cytosine_1based/filt/{sample}_CpG.tsv.gz",
#        "bismarkSE/CX/coverage2cytosine_1based/filt/{sample}_GpC.tsv.gz"
#     shell:
#        """
#        zmore {input[0]} | awk  -vOFS='\\t' '{{ if ($5>0 || $6>0)  {{print $1,$2,$5,$6}}}}' | sort -k 1,1 -k2,2n -S 1G --parallel 6 | gzip > bismarkSE/CX/coverage2cytosine_1based/filt/{wildcards.sample}_CpG.tsv.gz
#        zmore {input[1]} | awk  -vOFS='\\t' '{{ if ($5>0 || $6>0)  {{print $1,$2,$5,$6}}}}' | sort -k 1,1 -k2,2n -S 1G --parallel 6 | gzip > bismarkSE/CX/coverage2cytosine_1based/filt/{wildcards.sample}_GpC.tsv.gz
#        """

#rule binarize:
#     input:
#        expand("bismarkSE/CX/coverage2cytosine_1based/filt/{{sample}}_{filtbin}",filtbin=filtbin)
#        #cpg="bismarkSE/CX/coverage2cytosine_1based/filt/{sample}_CpG.tsv.gz",
#        #gpc="bismarkSE/CX/coverage2cytosine_1based/filt/{sample}_GpC.tsv.gz"
#     output:
#        "bismarkSE/CX/coverage2cytosine_1based/filt/binarised/{sample}_CpG.gz",
#        "bismarkSE/CX/coverage2cytosine_1based/filt/binarised/{sample}_GpC.gz"
#     conda:
#        "../envs/NMT_Binarize.yaml"
#     shell:
#        """
#        Rscript scripts/binarize.R --infile={input[0]} --outdir=bismarkSE/CX/coverage2cytosine_1based/filt/binarised --input_format=2
#        Rscript scripts/binarize.R --infile={input[1]} --outdir=bismarkSE/CX/coverage2cytosine_1based/filt/binarised --input_format=2
#        """

rule c2c_cov_filter_met:
     input:
        "bismarkSE/CX/coverage2cytosine_1based/{sample}_merged.NOMe.CpG.cov.gz"
     output:
        "bismarkSE/CX/coverage2cytosine_1based/filt/{sample}_CpG.tsv.gz"
     shell:
        """
        bash scripts/coverage2cytosine_cov_filter.sh -i {input} -o {output} -c 6 
        """

rule c2c_cov_filter_acc:
     input:
        "bismarkSE/CX/coverage2cytosine_1based/{sample}_merged.NOMe.GpC.cov.gz"
     output:
        "bismarkSE/CX/coverage2cytosine_1based/filt/{sample}_GpC.tsv.gz"
     shell:
        """
        bash scripts/coverage2cytosine_cov_filter.sh -i {input} -o {output} -c 6 
        """

rule binarize_met:
     input:
        "bismarkSE/CX/coverage2cytosine_1based/filt/{sample}_CpG.tsv.gz"
     output:
        "bismarkSE/CX/coverage2cytosine_1based/filt/binarised/{sample}_CpG.gz"
     conda:
        "../envs/NMT_Binarize.yaml"
     shell:
        """
        Rscript scripts/binarize.R --infile={input} --outdir=bismarkSE/CX/coverage2cytosine_1based/filt/binarised --input_format=2
        """

rule binarize_acc:
     input:
        "bismarkSE/CX/coverage2cytosine_1based/filt/{sample}_GpC.tsv.gz"
     output:
        "bismarkSE/CX/coverage2cytosine_1based/filt/binarised/{sample}_GpC.gz"
     conda:
        "../envs/NMT_Binarize.yaml"
     shell:
        """
        Rscript scripts/binarize.R --infile={input} --outdir=bismarkSE/CX/coverage2cytosine_1based/filt/binarised --input_format=2
        """

rule create_reports:
     input:
        expand("bismarkSE/CX/coverage2cytosine_1based/filt/binarised/{sample}_CpG.gz", sample = SAMPLES)
     output:
        "tables/bismarkSE_mapping_report.txt"
     shell:
        "bash scripts/reports.sh"
