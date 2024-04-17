while getopts a:b:c:d:o:p: flag
do
    case "${flag}" in
        a) input1=${OPTARG};;
        b) input2=${OPTARG};;
        c) input3=${OPTARG};;
        d) input4=${OPTARG};;
        o) output=${OPTARG};;
        p) prefix=${OPTARG};;
    esac
done

#echo -e "Sample TotalReads MappedReads Total_Cs CpG_Num CHG_Num CHH_Num Dedup_Num Dedup_Percent" | tr ' ' '\t'  > $output
for FILE in $input1 $input2 $input3 $input4; do
        Reads=`grep "^Sequences analysed in total" ${FILE} | awk '{print $5}'`
        Mapped=`grep "^Number of alignments with" ${FILE} | awk '{print $13}'`
        Cs=`grep "^Total number of C's analysed:" ${FILE} | awk '{print $6}'`
        CpG_methyl=`grep "^Total methylated C's in CpG context:" ${FILE} | awk '{print $7}'`
        CHG_methyl=`grep "^Total methylated C's in CHG context:" ${FILE} | awk '{print $7}'`
        CHH_methyl=`grep "^Total methylated C's in CHH context:" ${FILE} | awk '{print $7}'`
        CpG_unmethyl=`grep "^Total unmethylated C's in CpG context:" ${FILE} | awk '{print $7}'`
        CHG_unmethyl=`grep "^Total unmethylated C's in CHG context:" ${FILE} | awk '{print $7}'`
        CHH_unmethyl=`grep "^Total unmethylated C's in CHH context:" ${FILE} | awk '{print $7}'`
        echo ${prefix} $Reads $Mapped $Cs $CpG_methyl $CHG_methyl $CHH_methyl $CpG_unmethyl $CHG_unmethyl $CHH_unmethyl | tr ' ' '\t' >> ${output}
	echo ${prefix}
done
