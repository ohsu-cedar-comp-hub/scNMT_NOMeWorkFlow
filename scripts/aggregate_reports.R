library(data.table)
library(scales)

temp <- list.files ("bismarkSE", full.names = TRUE, pattern = "final.txt$")
out <- rbindlist(lapply(temp, fread), fill = TRUE)

dedup <- list.files("bismarkSE/dedup", full.names = TRUE, pattern = "final.txt$")
dedup_out <- rbindlist(lapply(dedup, fread), fill = TRUE)

colnames(out) <- c("Sample", "TotalReads", "MappedReads", "Total_Cs", "CpG_Methyl", "CHG_Methyl", "CHH_Methyl", "CpG_Unmethyl", "CHG_Unmethyl", "CHH_Unmethyl")
colnames(dedup_out) <- c("Sample", "Dedup_Num", "Dedup_Percent")

test_out <- aggregate(.~Sample, data=out, FUN=sum)

report_tab <- merge(test_out, dedup_out, by = "Sample")

report_tab$PercentMapped <- percent(report_tab$MappedReads/report_tab$TotalReads)
report_tab$CpG <- percent(report_tab$CpG_Methyl/(report_tab$CpG_Methyl + report_tab$CpG_Unmethyl))
report_tab$CHG <- percent(report_tab$CHG_Methyl/(report_tab$CHG_Methyl + report_tab$CHG_Unmethyl))
report_tab$CHH <- percent(report_tab$CHH_Methyl/(report_tab$CHH_Methyl + report_tab$CHH_Unmethyl))

write.table(report_tab, "tables/bismarkSE_mapping_report.txt", sep = "\t", row.names = F, quote = F)
