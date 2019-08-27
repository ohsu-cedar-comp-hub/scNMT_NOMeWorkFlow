rule trimming:
    input:
        fwd = "samples/raw/{sample}_R1.fastq",
        rev = "samples/raw/{sample}_R2.fastq"
    output:
        "samples/trim/{sample}.fastq_trimming_report.txt",
        "samples/trim/{sample}_R1_val_1_fastqc.html",
        "samples/trim/{sample}_R1_val_1_fastqc.zip",
        "samples/trim/{sample}_R1_val_1.fq",
        "samples/trim/{sample}_R2.fastq_trimming_report.txt",
        "samples/trim/{sample}_R2_val_2_fastqc.html",
        "samples/trim/{sample}_R2_val_2_fastqc.zip",
        "samples/trim/{sample}_R2_val_2.fq"
    conda:
        "../envs/align.yaml"
    shell:
        """trim_galore --paired --clip_R1 6 --clip_R2 6 --trim1 --dont_gzip --fastqc --output_dir samples/trim {input.fwd} {input.rev}"""

rule mapping
    input:
        fwd = "samples/trim/{sample}_R1_val_2.fq",
        rev = "samples/trim/{sample}_R2_val_2.fq",
	name = {sample}
    output:
        "bismark/{wildcards.sample}_R1.bam",
        "bismark/{wildcards.sample}_R1_report.txt",
        "bismark/{wildcards.sample}_R2.bam",
        "bismark/{wildcards.sample}_report.txt"
    params:	
        bismark_index = config["bismark_index"]
	bismark_outdir = config["map_dir"]
    conda:
        "../envs/align.yaml"
    shell:
        """bismark --basename {wildcards.sample}_R1 --output_dir {params.bismark_outdir} --non_directional --bowtie1 --gzip -n 1 {params.bismark_index} {input.fwd}"""
        """bismark --basename {wildcards.sample}_R2 --output_dir {params.bismark_outdir} --non_directional --bowtie1 --gzip -n 1 {params.bismark_index} {input.rev}"""

rule deduplcate_bam
    input:
        fwd = "bismarkSE/{sample}_R1.bam",
        rev = "bismarkSE/{sample}_R2.bam",
	name = {sample}
    output:
        "bismarkSE/dedup/{sample}_R1.deduplicated.bam",
        "bismarkSE/dedup/{sample}_R2.deduplicated.bam"
    conda:
        "../envs/align.yaml"
    shell:
        """deduplicate_bismark --single --output_dir bismarkSE/dedup --bam {input.fwd}"""
        """deduplicate_bismark --single --output_dir bismarkSE/dedup --bam {input.rev}"""

rule merge_bams
    input:
        fwd = "bismarkSE/{sample}_R1.bam",
        rev = "bismarkSE/{sample}_R2.bam",
	name = {sample}
    output:
        "bismarkSE/dedup/{sample}_R1.deduplicated.bam",
        "bismarkSE/dedup/{sample}_R2.deduplicated.bam"
    conda:
        "../envs/align.yaml"
    shell:
        """deduplicate_bismark --single --output_dir bismarkSE/dedup --bam {input.fwd}"""
        """deduplicate_bismark --single --output_dir bismarkSE/dedup --bam {input.rev}"""


rule methylation_extractor
    input:
        fwd = "bismarkSE/dedup/{sample}_R1.deduplicated.bam",
        rev = "bismarkSE/dedup/{sample}_R2.deduplicated.bam"
	name = {sample}
    output:
        "bismarkSE/dedup/{sample}_R1.deduplicated.header.sam",
	"bismarkSE/dedup/{sample}.bam",
        "bismarkSE/dedup/{sample}_R2.deduplicated.bam"
    conda:
        "../envs/align.yaml"
    shell:
            """samtools view -H {input.fwd} > bismarkSE/dedup/{input.name}_R1.deduplicated.header.sam"""
	    """samtools merge -h bismarkSE/dedup/{input.name}_R1.deduplicated.header.sam bismarkSE/dedup/{input.name}.bam {input.fwd} {input.rev}"""
	    """bismark_methylation_extractor --output bismarkSE/CX --CX --parallel 4 --gzip --bedGraph bismarkSE/dedup/{input.name}.bam"""

rule coverage2cytosine
    input:
        cov = bismarkSE/CX/{sample}.bismark.cov.gz,
	name = {sample}
    output:
        "bismarkSE/CX/coverage2cytosine_1based/{sample}.NOMe.CpG.cov.gz",
	"bismarkSE/CX/coverage2cytosine_1based/{sample}.NOMe.CpG_report.txt.gz",
        "bismarkSE/CX/coverage2cytosine_1based/{sample}.NOMe.GpC.cov.gz"
	"bismarkSE/CX/coverage2cytosine_1based/{sample}.NOMe.GpC_report.txt.gz",	
   params:	
        bismark_index = config["bismark_index"]
    conda:
        "../envs/align.yaml"
    shell:
            """coverage2cytosine --nome-seq --gzip --output {input.name} --dir bismarkSE/CX/coverage2cytosine_1based --genome_folder {params.bismark_index} {input.cov}"""
