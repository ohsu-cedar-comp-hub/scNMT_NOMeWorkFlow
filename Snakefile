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

bismark_ext = ['pe.bam', 'PE_report.txt']
CX_ext1 = ['CHG','CHH','CpG']
CX_ext2 = ['CTOT','OT','CTOB','OB']
CX_report = ['.bismark.cov.gz','_splitting_report.txt','.bedGraph.gz','.M-bias.txt']
c2c = ['CpG_report.txt.gz','GpC_report.txt.gz']

## generates log dirs from json file
with open('cluster.json') as json_file:
    json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())


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

## just print statement
for sample in SAMPLES:
    message("Sample " + sample + " will be processed")

### just need the one at the end not all of them but good practice
### everytime output has variable need to expand if no single file no expansion
rule all:
    input:
        expand("samples/trim/{sample}_R1.fastq.gz_trimming_report.txt", sample = SAMPLES),
	expand("samples/trim/{sample}_R2.fastq.gz_trimming_report.txt", sample = SAMPLES),
        expand("samples/trim/{sample}_R1_val_1_fastqc.html", sample = SAMPLES),
        expand("samples/trim/{sample}_R1_val_1_fastqc.zip", sample = SAMPLES),
	expand("samples/trim/{sample}_R2_val_2_fastqc.html", sample = SAMPLES),
        expand("samples/trim/{sample}_R2_val_2_fastqc.zip", sample = SAMPLES),
        expand("samples/trim/{sample}_R1_val_1.fq.gz", sample = SAMPLES),
	expand("samples/trim/{sample}_R2_val_2.fq.gz", sample = SAMPLES),
	expand("bismarkSE/{sample}_R1.bam", sample = SAMPLES),
	expand("bismarkSE/{sample}_R2.bam", sample = SAMPLES),
	expand("bismarkSE/dedup/{sample}_R1.deduplicated.bam", sample = SAMPLES),
	expand("bismarkSE/dedup/{sample}_R2.deduplicated.bam", sample = SAMPLES),
	expand("bismarkSE/dedup/{sample}_merged.bam", sample = SAMPLES),
	expand("bismarkSE/{sample}_R1_SE_report.txt", sample = SAMPLES),
	expand("bismarkSE/{sample}_R2_SE_report.txt", sample = SAMPLES),
	expand("bismarkSE/dedup/{sample}_R1.deduplication_report.txt", sample = SAMPLES),
	expand("bismarkSE/dedup/{sample}_R2.deduplication_report.txt", sample = SAMPLES),
	expand("bismarkSE/CX/{CX_ext1}_{CX_ext2}_{sample}_merged.txt.gz", sample = SAMPLES, CX_ext1 = CX_ext1, CX_ext2 = CX_ext2),
	expand("bismarkSE/CX/{sample}_merged{CX_report}", sample = SAMPLES, CX_report = CX_report),
	expand("bismarkSE/CX/coverage2cytosine_1based/{sample}_merged.NOMe.{c2c}", sample = SAMPLES, c2c = c2c),
	expand("bismarkSE/CX/coverage2cytosine_1based/filt/binarised/{sample}_CpG.gz", sample = SAMPLES),
	expand("bismarkSE/CX/coverage2cytosine_1based/filt/binarised/{sample}_GpC.gz", sample = SAMPLES)
	

## must include all rules
include: "rules/bismart.smk"

## results ask Joey about result dirs not use if it is required
