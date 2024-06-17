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
acc_prom = config["acc_prom"]
met_prom = config["met_prom"]

SAMPLES, = glob_wildcards("samples/raw/{sample}_R1.fastq.gz")

CX_ext1 = ['CHG','CHH','CpG']
CX_ext2 = ['CTOT','OT','CTOB','OB']
CX_report = ['_splitting_report.txt','.bedGraph.gz','.M-bias.txt']
c2c = ['CpG_report.txt.gz','GpC_report.txt.gz']
filt = ['CpG.tsv.gz','GpC.tsv.gz']
filtbin = ['CpG.gz','GpC.gz']
filtbindir = ['bismarkSE/CX/coverage2cytosine_1based/filt/binarised']
metORacc = ['met','acc']
regions = ['body', 'promoters']

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
#        expand("samples/trim/{sample}_R{replicate}.fastq.gz_trimming_report.txt", sample = SAMPLES, replicate=[1, 2]),
#        expand("samples/trim/{sample}_R{replicate}_val_{replicate}_fastqc.html", sample = SAMPLES, replicate=[1, 2]),
#        expand("samples/trim/{sample}_R{replicate}_val_{replicate}_fastqc.zip", sample = SAMPLES, replicate=[1, 2]),
        expand("samples/trim/{sample}_merged_R2_val_2_fastqc.zip", sample = SAMPLES),
	expand("bismarkSE/{sample}_R1.{sample}_merged_R1_val_1_bismark_bt2_SE_report_final.txt", sample = SAMPLES),
	expand("bismarkSE/{sample}_R2.{sample}_merged_R2_val_2_bismark_bt2_SE_report_final.txt", sample = SAMPLES),
#	expand("bismarkSE/dedup/{sample}_R1.{sample}_merged_R1_val_1_bismark_bt2.deduplication_report.txt", sample = SAMPLES),
#	expand("bismarkSE/dedup/{sample}_R2.{sample}_merged_R2_val_2_bismark_bt2.deduplication_report.txt", sample = SAMPLES),
	expand("bismarkSE/dedup/{sample}_R1_deduplication_report_final.txt", sample = SAMPLES),
	expand("bismarkSE/dedup/{sample}_R2_deduplication_report_final.txt", sample = SAMPLES),
	expand("bismarkSE/dedup/{sample}_merged.bam", sample = SAMPLES),
	expand("bismarkSE/CX/{CX_ext1}_{CX_ext2}_{sample}_merged.txt.gz", sample = SAMPLES, CX_ext1 = CX_ext1, CX_ext2 = CX_ext2),
	expand("bismarkSE/CX/{sample}_merged{CX_report}", sample = SAMPLES, CX_report = CX_report),
	expand("bismarkSE/CX/coverage2cytosine_1based/{sample}_merged.NOMe.{c2c}", sample = SAMPLES, c2c = c2c),
	expand("bismarkSE/CX/coverage2cytosine_1based/filt/binarised/{sample}_GpC.gz", sample = SAMPLES),
        expand("plots/Mbias/{sample}_merged_Mbias_plot.pdf", sample = SAMPLES),
        expand("data/anno/{sample_name}_called_peaks.bed", sample_name = config["project_id"]),
        "tables/bismarkSE_mapping_report.txt",
        "plots/counts_stats/coverageBar.pdf",
        "plots/counts_stats/GWmethylRate.pdf",
        "plots/counts_stats/GWaccessRate.pdf",
        "tables/sample_read_report.txt",
        "plots/counts_stats/meanMethyl_vs_meanAccess.pdf",
        "plots/counts_stats/meanRatevsCoverage.pdf",
        "plots/counts_stats/meanRateLinePlot.pdf",
        "plots/counts_stats/coverageLinePlot.pdf",
        "plots/counts_stats/meanRateBoxPlot.pdf",
        "plots/counts_stats/coverageBoxPlot.pdf",
        "plots/counts_stats/coverageDensityPlot.pdf",
        "plots/counts_stats/qc_accessCoverageBar.pdf",
        "plots/counts_stats/qc_methylCoverageBar.pdf",
        "plots/met_acc_QC/qc_accessCoverageBar.pdf",
        "plots/met_acc_QC/qc_methylCoverageBar.pdf",
        "plots/met_acc_QC/qc_methylMeanCoverageBar.pdf",
        "plots/met_acc_QC/qc_accessMeanCoverageBar.pdf",
        "plots/met_acc_QC/qc_meanRateBoxPlot.pdf",
        "plots/met_acc_QC/qc_coverageBoxPlot.pdf",
        "tables/sample_stats_qcPass.txt",
        "tables/sample_read_report_qcPASS.txt",
        "plots/profiles/accessibility_at_promoters.RData",	
	"plots/profiles/accessibility_at_promoters.pdf",
	"plots/profiles/accessibility_average_promoters.pdf",
#	"data/anno/body.bed",
#	"data/anno/promoters.bed",
	expand("data/anno/{regions}{acc_prom}.bed", acc_prom = acc_prom, regions = regions),
	expand("data/anno/{regions}{met_prom}.bed", met_prom = met_prom, regions = regions),
	"plots/profiles/methylation_at_promoters.RData",
        "plots/profiles/methylation_at_promoters.pdf",
	"plots/profiles/methylation_average_promoters.pdf",
#	expand("data/{metORacc}/body.tsv.gz", metORacc = metORacc),
#        expand("data/{metORacc}/CTCF.tsv.gz", metORacc = metORacc),
#        expand("data/{metORacc}/Enhancer.tsv.gz", metORacc = metORacc),
#        expand("data/{metORacc}/MCF7_ER_peaks.tsv.gz", metORacc = metORacc),
#        expand("data/{metORacc}/MCF7_H3K27ac_peaks.tsv.gz", metORacc = metORacc),
#        expand("data/{metORacc}/Repressed.tsv.gz", metORacc = metORacc),
	"plots/anno_rateVarboxplots.pdf",
	"plots/corr/acc_met_correlations_loci.pdf",
        "plots/corr/acc_met_correlations_loci.tsv",
        "plots/UMAPs/allUMAPsDone.txt"

## must include all rules
include: "rules/bismart.smk"
include: "rules/QC.smk"
include: "rules/profiles.smk"

## results ask Joey about result dirs not use if it is required
