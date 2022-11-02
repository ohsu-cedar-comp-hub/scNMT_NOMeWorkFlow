echo -e "Sample TotalReads MappedReads PercentMapped Total_Cs CpG CHG CHH Dedup_Num Dedup_Percent" | tr ' ' '\t'  > tables/bismarkSE_mapping_report.txt
grep -v Project data/metadata.txt  | cut -f 2 | while read Name; do
    for x in 1 2; do
        Reads=`grep "^Sequences analysed in total" bismarkSE/${Name}_R${x}.${Name}_R${x}_val_${x}_bismark_bt2_SE_report.txt | awk '{print $5}'`
        Mapped=`grep "^Number of alignments with" bismarkSE/${Name}_R${x}.${Name}_R${x}_val_${x}_bismark_bt2_SE_report.txt | awk '{print $13}'`
        Percent=`grep "^Mapping efficiency" bismarkSE/${Name}_R${x}.${Name}_R${x}_val_${x}_bismark_bt2_SE_report.txt | awk '{print $3}'`
        Cs=`grep "^Total number of C's analysed:" bismarkSE/${Name}_R${x}.${Name}_R${x}_val_${x}_bismark_bt2_SE_report.txt | awk '{print $6}'`
        CpG=`grep "^C methylated in CpG context:" bismarkSE/${Name}_R${x}.${Name}_R${x}_val_${x}_bismark_bt2_SE_report.txt | awk '{print $6}'`
        CHG=`grep "^C methylated in CHG context:" bismarkSE/${Name}_R${x}.${Name}_R${x}_val_${x}_bismark_bt2_SE_report.txt | awk '{print $6}'`
        CHH=`grep "^C methylated in CHH context:" bismarkSE/${Name}_R${x}.${Name}_R${x}_val_${x}_bismark_bt2_SE_report.txt | awk '{print $6}'`
        dedup_num=`grep "^Total number duplicated alignments removed" bismarkSE/dedup/${Name}_R${x}.${Name}_R${x}_val_${x}_bismark_bt2.deduplication_report.txt | awk '{print $6}'`
        dedup_per=`grep "^Total number duplicated alignments removed" bismarkSE/dedup/${Name}_R${x}.${Name}_R${x}_val_${x}_bismark_bt2.deduplication_report.txt | awk '{print $7}' | tr -d '()'`
        echo ${Name}_R${x} $Reads $Mapped $Percent $Cs $CpG $CHG $CHH $dedup_num $dedup_per | tr ' ' '\t' >> tables/bismarkSE_mapping_report.txt
	echo ${Name}
    done;done
