rule counts_stats:
     input: rules.create_reports.output
     output:
        z="plots/counts_stats/coverageBar.pdf",
        y="plots/counts_stats/GWmethylRate.pdf",
        x="plots/counts_stats/GWaccessRate.pdf",
        a="tables/sample_stats.txt",
        b="tables/sample_read_report.txt",
        c="plots/counts_stats/meanMethyl_vs_meanAccess.pdf",
        d="plots/counts_stats/meanRatevsCoverage.pdf",
        e="plots/counts_stats/meanRateLinePlot.pdf",
        f="plots/counts_stats/coverageLinePlot.pdf",
        g="plots/counts_stats/meanRateBoxPlot.pdf",
        h="plots/counts_stats/coverageBoxPlot.pdf",
        i="plots/counts_stats/coverageDensityPlot.pdf",
        j="plots/counts_stats/qc_accessCoverageBar.pdf",
        k="plots/counts_stats/qc_methylCoverageBar.pdf"
     shell:
        "Rscript scripts/counts_stats.R --outdir=plots/counts_stats"

rule met_acc_QC:
     input: rules.counts_stats.output.a
     output:
        "plots/met_acc_QC/qc_accessCoverageBar.pdf",
        "plots/met_acc_QC/qc_methylCoverageBar.pdf",
        "plots/met_acc_QC/qc_methylMeanCoverageBar.pdf",
        "plots/met_acc_QC/qc_accessMeanCoverageBar.pdf",
        "plots/met_acc_QC/qc_meanRateBoxPlot.pdf",
        "plots/met_acc_QC/qc_coverageBoxPlot.pdf",
        "tables/sample_stats_qcPass.txt",
        "tables/sample_read_report_qcPASS.txt",
     shell:
        "Rscript scripts/met_acc_QC.R --outdir=plots/met_acc_QC"
