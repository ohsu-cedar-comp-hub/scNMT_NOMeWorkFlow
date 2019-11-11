#!/bin/bash
set -o errexit
set -o pipefail

################################################################################
# Set default values
################################################################################
Cores=3
tempDir="/tmp/"

################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
    echo >&2 "
Author: Ashley R Woodfin
$(basename $0) - Filter bismark CG methylation and GC accessibility GpC.cov.gz and CpG.cov.gz files
- <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
- Remove sites that were not covered by any read 
- Sometimes Bismark spits out a file were the non-observed sites have been removed, but sometimes it does not.
- Output is <chromosome> <start position> <count methylated> <count unmethylated>
                
USAGE: $(basename $0) -i <Input file> -o <Output directory> [OPTIONS]
 -i     Input file ( {CpG,GpC}_report.txt.gz format )      [ required ]
 -o     Output file path (must include context, and .gz)   [ required ]
 -c     Total cores (including 3 already used)             [ default: $Cores ]
 -t     Directory to which temporary files will be written [ default: $tempDir ]
NOTES:
The program will use 3 cores, but you can add additional for the sorting.
If 6 provided, 4 used in sorting step. 
"
    exit 1
fi

################################################################################
# Parse input and check for errors
################################################################################

while getopts "i:o:c:t:" o
do
    case "$o" in
        i) Infile="$OPTARG";;
        o) Outfile="$OPTARG";;
	c) Cores="$OPTARG";;
	t) tempDir="$OPTARG";;
       \?) exit 1;;
    esac
done

if [ -z "$Infile"  -o -z "$Outfile" ]; then
    echo >&2 "ERROR: -i -o required!"; exit 1
fi

################################################################################
# 1. Create outdir if needed
# 2. Filter coverage reports
################################################################################
Outdir=`dirname $Outfile`
if [ ! -d "$Outdir" ]; then
    mkdir -p $Outdir
fi

TMP=$(mktemp -d -p "${tempDir}")

#Sample=`basename $Infile | sed 's/.NOMe./_/g' | tr -d '.cov.gz' | cut -d "_" -f 1-2`
#Context=`basename $Infile | sed 's/.NOMe./_/g' | tr -d '.cov.gz' | cut -d "_" -f 3`
echo "Processing ${Infile}"
zmore $Infile | awk  -vOFS='\t' '{ if ($5>0 || $6>0)  {print $1,$2,$5,$6}}' | sort -k 1,1 -k2,2n -S 1G --parallel $(echo $Cores | awk '{print $1-3}') | gzip > ${Outfile}

##############################################################################################
# exit
##############################################################################################

rm -rf $TMP

exit 0
