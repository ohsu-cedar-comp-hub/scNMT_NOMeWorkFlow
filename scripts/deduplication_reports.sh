while getopts a:o:p: flag
do
    case "${flag}" in
        a) input1=${OPTARG};;
        o) output=${OPTARG};;
        p) prefix=${OPTARG};;
    esac
done

dedup_num=`grep "^Total number duplicated alignments removed" $input1 | awk '{print $6}'`
dedup_per=`grep "^Total number duplicated alignments removed" $input1 | awk '{print $7}' | tr -d '()'`
echo $prefix $dedup_num $dedup_per | tr ' ' '\t' >> $output
echo $prefix

