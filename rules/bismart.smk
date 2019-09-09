rule trimming:
    input:
        "samples/raw/{sample}_R1.fastq.gz",
        "samples/raw/{sample}_R2.fastq.gz"
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
        """trim_galore --paired --clip_R1 6 --clip_R2 6 --trim1 --gzip --fastqc --output_dir samples/trim {input[0]} {input[1]}"""

rule mapping_R1:
    input:
        fwd = "samples/trim/{sample}_R1_val_1.fq.gz",
    output:
        "bismarkSE/{sample}_R1.bam",
        "bismarkSE/{sample}_R1_SE_report.txt",
    params:	
        bismark_index = config["bismark_index"],
    conda:
        "../envs/methylome.yaml"
    shell:
        """bismark --basename {wildcards.sample}_R1 --output_dir bismarkSE --non_directional --score_min=L,-50,-0.2 --gzip -n 1 {params.bismark_index} {input.fwd}"""
        
rule mapping_R2:
    input:
        rev = "samples/trim/{sample}_R2_val_2.fq.gz",
    output:
        "bismarkSE/{sample}_R2.bam",
        "bismarkSE/{sample}_R2_SE_report.txt",
    params:	
        bismark_index = config["bismark_index"],
    conda:
        "../envs/methylome.yaml"
    shell:
        """bismark --basename {wildcards.sample}_R2 --output_dir bismarkSE --non_directional --score_min=L,-50,-0.2 --gzip -n 1 {params.bismark_index} {input.rev}"""

rule deduplcate_bam_R1:
    input:
        fwd = "bismarkSE/{sample}_R1.bam",
    output:
        "bismarkSE/dedup/{sample}_R1.deduplicated.bam",
    conda:
        "../envs/methylome.yaml"
    shell:
        """deduplicate_bismark --single --output_dir bismarkSE/dedup --bam {input.fwd}"""

rule deduplcate_bam_R2:
    input:
        rev = "bismarkSE/{sample}_R2.bam",
    output:
        "bismarkSE/dedup/{sample}_R2.deduplicated.bam"
    conda:
        "../envs/methylome.yaml"
    shell:
        """deduplicate_bismark --single --output_dir bismarkSE/dedup --bam {input.rev}"""

rule methylation_extractor:
    input:
        fwd = "bismarkSE/dedup/{sample}_R1.deduplicated.bam",
        rev = "bismarkSE/dedup/{sample}_R2.deduplicated.bam",
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
	bismark_methylation_extractor -s --output bismarkSE/CX --CX --parallel 4 --gzip --bedGraph bismarkSE/dedup/{wildcards.sample}_merged.bam
	"""

rule coverage2cytosine:
    input:
        cov = "bismarkSE/CX/{sample}_merged.bismark.cov.gz", 
    output:
        "bismarkSE/CX/coverage2cytosine_1based/{sample}_merged.NOMe.CpG.cov.gz",
	"bismarkSE/CX/coverage2cytosine_1based/{sample}_merged.NOMe.CpG_report.txt.gz",
        "bismarkSE/CX/coverage2cytosine_1based/{sample}_merged.NOMe.GpC.cov.gz",
	"bismarkSE/CX/coverage2cytosine_1based/{sample}_merged.NOMe.GpC_report.txt.gz",	
    params:	
        bismark_index = config["bismark_index"]
    conda:
        "../envs/methylome.yaml"
    shell:
        """coverage2cytosine --nome-seq --gzip --output {wildcards.sample}_merged --dir bismarkSE/CX/coverage2cytosine_1based --genome_folder {params.bismark_index} {input.cov}"""
