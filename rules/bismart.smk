rule merge_fastqs:
    input:
        read1 = "samples/raw/{sample}_R1.fastq.gz",
        read2 = "samples/raw/{sample}_R2.fastq.gz"
    output:
        temp("samples/merged/{sample}_merged_R1.fq.gz"),
        temp("samples/merged/{sample}_merged_R2.fq.gz")
    shell:
        """
        mkdir -p samples/merged
        cat samples/raw/{wildcards.sample}*R1*.fastq.gz > samples/merged/{wildcards.sample}_merged_R1.fq.gz
        cat samples/raw/{wildcards.sample}*R2*.fastq.gz > samples/merged/{wildcards.sample}_merged_R2.fq.gz
        """

rule trimming:
    input:
        read1 = "samples/merged/{sample}_merged_R1.fq.gz",
        read2 = "samples/merged/{sample}_merged_R2.fq.gz"
    output:
        "samples/trim/{sample}_merged_R1.fq.gz_trimming_report.txt",
        "samples/trim/{sample}_merged_R1_val_1_fastqc.html",
        "samples/trim/{sample}_merged_R1_val_1_fastqc.zip",
        temp("samples/trim/{sample}_merged_R1_val_1.fq.gz"),
        "samples/trim/{sample}_merged_R2.fq.gz_trimming_report.txt",
        "samples/trim/{sample}_merged_R2_val_2_fastqc.html",
        "samples/trim/{sample}_merged_R2_val_2_fastqc.zip",
        temp("samples/trim/{sample}_merged_R2_val_2.fq.gz")
    conda:
        "../envs/trimG.yaml"
    shell:
        """trim_galore --paired --clip_R1 16 --clip_R2 10 --trim1 --gzip --fastqc --output_dir samples/trim {input.read1} {input.read2}"""

rule split_fastqs:
    input:
        read1 = "samples/trim/{sample}_merged_R1_val_1.fq.gz",
        read2 = "samples/trim/{sample}_merged_R2_val_2.fq.gz"
    output:
        read1_1 = temp("samples/split/{sample}_merged_R1_val_1_1.fq.gz"),
        read1_2 = temp("samples/split/{sample}_merged_R1_val_1_2.fq.gz"),
        read1_3 = temp("samples/split/{sample}_merged_R1_val_1_3.fq.gz"),
        read1_4 = temp("samples/split/{sample}_merged_R1_val_1_4.fq.gz"),
        read2_1 = temp("samples/split/{sample}_merged_R2_val_2_1.fq.gz"),
        read2_2 = temp("samples/split/{sample}_merged_R2_val_2_2.fq.gz"),
        read2_3 = temp("samples/split/{sample}_merged_R2_val_2_3.fq.gz"),
        read2_4 = temp("samples/split/{sample}_merged_R2_val_2_4.fq.gz")
    conda:
        "../envs/split_fastqs.yaml"
    shell:
        """
        mkdir -p samples/split
        fastqsplitter -i {input.read1} -o {output.read1_1} -o {output.read1_2} -o {output.read1_3} -o {output.read1_4}
        fastqsplitter -i {input.read2} -o {output.read2_1} -o {output.read2_2} -o {output.read2_3} -o {output.read2_4}
        """
        
rule mapping_R1_1:
    input:
        fwd = "samples/split/{sample}_merged_R1_val_1_1.fq.gz",
    output:
        temp("bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_1_bismark_bt2.bam"),
        "bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_1_bismark_bt2_SE_report.txt"
    params:	
        bismark_index = config["bismark_index"],
    conda:
        "../envs/methylome.yaml"
    shell:
        """bismark --prefix {wildcards.sample}_R1 --output_dir bismarkSE --temp_dir /mnt/scratch --non_directional --parallel=2 --score_min=L,0,-0.2 --gzip -n 1 {params.bismark_index} {input.fwd}"""
        
rule mapping_R2_1:
    input:
        rev = "samples/split/{sample}_merged_R2_val_2_1.fq.gz",
    output:
        temp("bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_1_bismark_bt2.bam"),
        "bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_1_bismark_bt2_SE_report.txt"
    params:	
        bismark_index = config["bismark_index"],
    conda:
        "../envs/methylome.yaml"
    shell:
        """bismark --prefix {wildcards.sample}_R2 --output_dir bismarkSE --temp_dir /mnt/scratch --non_directional --parallel=2 --score_min=L,0,-0.2 --gzip -n 1 {params.bismark_index} {input.rev}"""

rule mapping_R1_2:
    input:
        fwd = "samples/split/{sample}_merged_R1_val_1_2.fq.gz",
    output:
        temp("bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_2_bismark_bt2.bam"),
        "bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_2_bismark_bt2_SE_report.txt"
    params:	
        bismark_index = config["bismark_index"],
    conda:
        "../envs/methylome.yaml"
    shell:
        """bismark --prefix {wildcards.sample}_R1 --output_dir bismarkSE --temp_dir /mnt/scratch --non_directional --parallel=2 --score_min=L,0,-0.2 --gzip -n 1 {params.bismark_index} {input.fwd}"""
        
rule mapping_R2_2:
    input:
        rev = "samples/split/{sample}_merged_R2_val_2_2.fq.gz",
    output:
        temp("bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_2_bismark_bt2.bam"),
        "bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_2_bismark_bt2_SE_report.txt"
    params:	
        bismark_index = config["bismark_index"],
    conda:
        "../envs/methylome.yaml"
    shell:
        """bismark --prefix {wildcards.sample}_R2 --output_dir bismarkSE --temp_dir /mnt/scratch --non_directional --parallel=2 --score_min=L,0,-0.2 --gzip -n 1 {params.bismark_index} {input.rev}"""

rule mapping_R1_3:
    input:
        fwd = "samples/split/{sample}_merged_R1_val_1_3.fq.gz",
    output:
        temp("bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_3_bismark_bt2.bam"),
        "bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_3_bismark_bt2_SE_report.txt"
    params:	
        bismark_index = config["bismark_index"],
    conda:
        "../envs/methylome.yaml"
    shell:
        """bismark --prefix {wildcards.sample}_R1 --output_dir bismarkSE --temp_dir /mnt/scratch --non_directional --parallel=2 --score_min=L,0,-0.2 --gzip -n 1 {params.bismark_index} {input.fwd}"""
        
rule mapping_R2_3:
    input:
        rev = "samples/split/{sample}_merged_R2_val_2_3.fq.gz",
    output:
        temp("bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_3_bismark_bt2.bam"),
        "bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_3_bismark_bt2_SE_report.txt"
    params:	
        bismark_index = config["bismark_index"],
    conda:
        "../envs/methylome.yaml"
    shell:
        """bismark --prefix {wildcards.sample}_R2 --output_dir bismarkSE --temp_dir /mnt/scratch --non_directional --parallel=2 --score_min=L,0,-0.2 --gzip -n 1 {params.bismark_index} {input.rev}"""

rule mapping_R1_4:
    input:
        fwd = "samples/split/{sample}_merged_R1_val_1_4.fq.gz",
    output:
        temp("bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_4_bismark_bt2.bam"),
        "bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_4_bismark_bt2_SE_report.txt"
    params:	
        bismark_index = config["bismark_index"],
    conda:
        "../envs/methylome.yaml"
    shell:
        """bismark --prefix {wildcards.sample}_R1 --output_dir bismarkSE --temp_dir /mnt/scratch --non_directional --parallel=2 --score_min=L,0,-0.2 --gzip -n 1 {params.bismark_index} {input.fwd}"""
        
rule mapping_R2_4:
    input:
        rev = "samples/split/{sample}_merged_R2_val_2_4.fq.gz",
    output:
        temp("bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_4_bismark_bt2.bam"),
        "bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_4_bismark_bt2_SE_report.txt"
    params:	
        bismark_index = config["bismark_index"],
    conda:
        "../envs/methylome.yaml"
    shell:
        """bismark --prefix {wildcards.sample}_R2 --output_dir bismarkSE --temp_dir /mnt/scratch --non_directional --parallel=2 --score_min=L,0,-0.2 --gzip -n 1 {params.bismark_index} {input.rev}"""

rule merge_bismark_bams:
    input:
        bam1_1 = "bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_1_bismark_bt2.bam",
        bam1_2 = "bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_2_bismark_bt2.bam",
        bam1_3 = "bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_3_bismark_bt2.bam",
        bam1_4 = "bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_4_bismark_bt2.bam",
        bam2_1 = "bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_1_bismark_bt2.bam",
        bam2_2 = "bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_2_bismark_bt2.bam",
        bam2_3 = "bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_3_bismark_bt2.bam",
        bam2_4 = "bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_4_bismark_bt2.bam",
    output:
        fwd = "bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_bismark_bt2.bam",
        rev = "bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_bismark_bt2.bam",
    shell:
        """
        samtools merge {output.fwd} {input.bam1_1} {input.bam1_2} {input.bam1_3} {input.bam1_4}
        samtools merge {output.rev} {input.bam2_1} {input.bam2_2} {input.bam2_3} {input.bam2_4}
        """

rule merge_bismark_reports:
    input:
        rep1_1 = "bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_1_bismark_bt2_SE_report.txt",
        rep1_2 = "bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_2_bismark_bt2_SE_report.txt",
        rep1_3 = "bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_3_bismark_bt2_SE_report.txt",
        rep1_4 = "bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_4_bismark_bt2_SE_report.txt",
        rep2_1 = "bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_1_bismark_bt2_SE_report.txt",
        rep2_2 = "bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_2_bismark_bt2_SE_report.txt",
        rep2_3 = "bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_3_bismark_bt2_SE_report.txt",
        rep2_4 = "bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_4_bismark_bt2_SE_report.txt"
    output:
        rep1 = "bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_bismark_bt2_SE_report_final.txt",
        rep2 = "bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_bismark_bt2_SE_report_final.txt"
    shell:
        """
        bash scripts/merge_reports.sh -a {input.rep1_1} -b {input.rep1_2} -c {input.rep1_3} -d {input.rep1_4} -o {output.rep1} -p {wildcards.sample}_R1
        bash scripts/merge_reports.sh -a {input.rep2_1} -b {input.rep2_2} -c {input.rep2_3} -d {input.rep2_4} -o {output.rep2} -p {wildcards.sample}_R2
        """
        
rule deduplcate_bam_R1:
    input:
        fwd = "bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_bismark_bt2.bam",
    output:
        "bismarkSE/dedup/{sample}_R1.{sample}_merged_R1_val_1_bismark_bt2.deduplicated.bam",
	"bismarkSE/dedup/{sample}_R1.{sample}_merged_R1_val_1_bismark_bt2.deduplication_report.txt"
    conda:
        "../envs/methylome.yaml"
    shell:
        """deduplicate_bismark --single --output_dir bismarkSE/dedup --bam {input.fwd}"""

rule deduplcate_bam_R2:
    input:
        rev = "bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_bismark_bt2.bam",
    output:
        "bismarkSE/dedup/{sample}_R2.{sample}_merged_R2_val_2_bismark_bt2.deduplicated.bam",
	"bismarkSE/dedup/{sample}_R2.{sample}_merged_R2_val_2_bismark_bt2.deduplication_report.txt"
    conda:
        "../envs/methylome.yaml"
    shell:
        """deduplicate_bismark --single --output_dir bismarkSE/dedup --bam {input.rev}"""

rule deduplcate_reports:
    input:
       fwd = "bismarkSE/dedup/{sample}_R1.{sample}_merged_R1_val_1_bismark_bt2.deduplication_report.txt",
       rev = "bismarkSE/dedup/{sample}_R2.{sample}_merged_R2_val_2_bismark_bt2.deduplication_report.txt"
    output:
       fwd = "bismarkSE/dedup/{sample}_R1_deduplication_report_final.txt",
       rev = "bismarkSE/dedup/{sample}_R2_deduplication_report_final.txt"
    shell:
       """
       bash scripts/deduplication_reports.sh -a {input.fwd} -o {output.fwd} -p {wildcards.sample}_R1
       bash scripts/deduplication_reports.sh -a {input.rev} -o {output.rev} -p {wildcards.sample}_R2
       """

rule methylation_extractor:
    input:
        fwd = "bismarkSE/dedup/{sample}_R1.{sample}_merged_R1_val_1_bismark_bt2.deduplicated.bam",
        rev = "bismarkSE/dedup/{sample}_R2.{sample}_merged_R2_val_2_bismark_bt2.deduplicated.bam",
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

rule Mbias_plot:
      input:
         "bismarkSE/CX/{sample}_merged.M-bias.txt"
      output:
         "plots/Mbias/{sample}_merged_Mbias_plot.pdf"
      conda:
         "../envs/NMT_NOMe_QC.yaml"
      shell:
         "Rscript scripts/Mbias_plots.R --infile={input} --sample={wildcards.sample}_merged --outdir=plots/Mbias/"

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

#rule create_reports:
#     input:
#        expand("bismarkSE/CX/coverage2cytosine_1based/filt/binarised/{sample}_CpG.gz", sample = SAMPLES)
#     output:
#        "tables/bismarkSE_mapping_report.txt"
#     params:	
#        name = config["project_id"],
#        seq_type = config["seq_type"]
#     shell:
#        """
#        echo -e "ProjectID\tSampleID\tReadType" > data/metadata.txt
#        ls -1 bismarkSE/ | awk -F. '{{print $1}}' | awk -F_ '{{ print {params.name}"\t"$1"_"$2"_"$3"\t"{params.seq_type} }}' | sort | uniq >> data/metadata.txt
#        bash scripts/reports.sh
#        """

rule create_reports:
     input:
        expand("bismarkSE/CX/coverage2cytosine_1based/filt/binarised/{sample}_CpG.gz", sample = SAMPLES)
     output:
        "tables/bismarkSE_mapping_report.txt"
     params:	
        name = config["project_id"],
        seq_type = config["seq_type"]
     shell:
        """
        Rscript scripts/aggregate_reports.R
        echo -e "ProjectID\tSampleID\tReadType" > data/metadata.txt
        ls -1 bismarkSE/ | awk -F. '{{print $1}}' | awk -F_ '{{ print {params.name}"\t"$1"_"$2"\t"{params.seq_type} }}' | sort | uniq >> data/metadata.txt
        """


