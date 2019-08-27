__author__ = "Ashley Woodfin"
__email__ = "woodfin@ohsu.edu"
__license__ = "MIT"

"""Methylation data processing pipeline"""


import datetime
import sys
import os
import pandas as pd
import json


timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

configfile:"config.yaml"
project_id = config["project_id"]


SAMPLES, = glob_wildcards("samples/raw/{sample}_R1.fastq.gz")

ext = ['r','R1.pdf','R2.pdf','xls']
fastq_ext = ['1','2']
trim_ext = ['.fastq_trimming_report.txt', '_val_1_fastqc.html', '_val_2_fastqc.html', '_val_1_fastqc.zip', '_val_2_fastqc.zip', '_val_1.fq', '_val_2.fq']
bismark_ext = ['pe.bam', 'PE_report.txt']
fastqscreen_ext = ['html','png','txt']
insertion_and_clipping_prof_ext = ['r','R1.pdf','R2.pdf','xls']
inner_distance_ext = ['_freq.txt','_plot.pdf','_plot.r','.txt']
read_dist_ext = ['txt']
read_gc_ext = ['.xls','_plot.r','_plot.pdf']

## generates log dirs from json file
with open('cluster.json') as json_file:
    json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())
rule_dirs.pop(rule_dirs.index('__default__'))


for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
        log_out = os.path.join(os.getcwd(), 'logs', rule)
        os.makedirs(log_out)
        print(log_out)

## output dirs that are created
result_dirs = ['diffexp','tables']
for rule in result_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'results',rule)):
        log_out = os.path.join(os.getcwd(), 'results', rule)
        os.makedirs(log_out)
        print(log_out)


def message(mes):
    sys.stderr.write("|--- " + mes + "\n")

## used for omic rules for python that writes an rscript wont need
def format_plot_columns():
    factors = config['meta_columns_to_plot'].keys()
    reformat_factors = '"' + '","'.join(factors) + '"'
    return 'c({})'.format(reformat_factors)

## dont need unless using DEseq
def get_deseq2_threads(wildcards=None):
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(config["omic_meta_data"]) < 100 or few_coeffs else 6

## used to define contrasts
def get_contrast(wildcards):
    """Return each contrast provided in the configuration file"""
    return config["diffexp"]["contrasts"][wildcards.contrast]

## just print statement
for sample in SAMPLES:
    message("Sample " + sample + " will be processed")

### just need the one at the end not all of them but good practice
### everytime output has variable need to expand if no single file no expansion
rule all:
    input:
        expand("samples/trim/{sample}_R{fastq_ext}.fastq_trimming_report.txt", sample = SAMPLES, fastq_ext = fastq_ext),    
        expand("samples/trim/{sample}_R{fastq_ext}_val_{fastq_ext}_fastqc.html", sample = SAMPLES, fastq_ext = fastq_ext),
        expand("samples/trim/{sample}_R{fastq_ext}_val_{fastq_ext}_fastqc.zip", sample = SAMPLES, fastq_ext = fastq_ext),
        expand("samples/trim/{sample}_R{fastq_ext}_val_{fastq_ext}.fq", sample = SAMPLES, fastq_ext = fastq_ext),	
        expand("bismark/{sample}_{bismark_ext}", sample = SAMPLES, bismark_ext = bismark_ext)
	
rule all:
    input:
        expand("results/tables/{project_id}_STAR_mapping_statistics.txt", project_id = config['project_id']),
        expand("samples/fastqc/{sample}/{sample}_{fastq_ext}_t_fastqc.zip", sample = SAMPLES, fastq_ext = fastq_ext),
        expand("samples/fastqscreen/{sample}/{sample}_{fastq_ext}_t.good_screen.{fastqscreen_ext}", sample=SAMPLES, fastq_ext=fastq_ext, fastqscreen_ext=fastqscreen_ext),
        expand("samples/samtools_stats/{sample}.txt",sample=SAMPLES),
        expand("rseqc/insertion_profile/{sample}/{sample}.insertion_profile.{ext}",sample=SAMPLES, ext=insertion_and_clipping_prof_ext),
        expand("rseqc/inner_distance/{sample}/{sample}.inner_distance{ext}", sample = SAMPLES, ext = inner_distance_ext),
        expand("rseqc/clipping_profile/{sample}/{sample}.clipping_profile.{ext}", sample = SAMPLES, ext = insertion_and_clipping_prof_ext),
        expand("rseqc/read_distribution/{sample}/{sample}.read_distribution.{ext}", sample = SAMPLES, ext = read_dist_ext),
        expand("rseqc/read_GC/{sample}/{sample}.GC{ext}", sample = SAMPLES, ext = read_gc_ext),
        expand("samples/samtools_stats/{sample}.txt",sample=SAMPLES),
        "data/{project_id}_coverage.txt".format(project_id=config["project_id"]),
        "results/tables/{}_Normed_with_Ratio_and_Abundance.txt".format(config['project_id']),
        "results/diffexp/pca.pdf",
        expand(["results/diffexp/GOterms/{contrast}.diffexp.downFC.2.adjp.0.01_BP_GO.txt", "results/diffexp/GOterms/{contrast}.diffexp.upFC.2.adjp.0.01_BP_GO.txt", "results/diffexp/GOterms/{contrast}.diffexp.downFC.2.adjp.0.01.BP.pdf", "results/diffexp/GOterms/{contrast}.diffexp.upFC.2.adjp.0.01.BP.pdf","results/diffexp/GOterms/{contrast}.diffexp.downFC.2.adjp.0.01_BP_classic_5_all.pdf","results/diffexp/GOterms/{contrast}.diffexp.upFC.2.adjp.0.01_BP_classic_5_all.pdf"], contrast = config["diffexp"]["contrasts"]),
        expand("results/diffexp/{contrast}.diffexp.01.VolcanoPlot.pdf", contrast = config["diffexp"]["contrasts"]),
        expand(["results/diffexp/glimma-plots/{contrast}.ma_plot.html","results/diffexp/glimma-plots/{contrast}.volcano_plot.html"],contrast = config["diffexp"]["contrasts"]),
        "results/diffexp/glimma-plots/{project_id}.mds_plot.html".format(project_id=project_id,

## must include all rules
include: "rules/align_rmdp.smk"
include: "rules/omic_qc.smk"
include: "rules/deseq.smk"

## results ask Joey about result dirs not use if it is required