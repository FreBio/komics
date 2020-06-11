#!/bin/bash

if [ $# -eq 0 ]; then
    echo ""
    echo "FVDB"
    echo ""
    echo "Calculate read depth statistics from SAM/BAM files"
    echo ""
    echo "Usage: mapping_stats.sh inputfile outputfile"
    echo "Requires: SAMTOOLS"
    echo ""
else

samfile=$1
textfile=$2
genomehapdepth=$(($3/2))

samtools depth -a $samfile > $samfile\.depth
MinicircleContigs=$(awk '{print $1}' $samfile\.depth | sort | uniq)

echo "CONTIG" "MEDIAN.DEPTH" "MINICIRCLE.COPY.NUMBER">> $textfile
for mini in MinicircleContigs
	depth=$(grep "^$mini" $samfile\.depth | awk ' { a[i++]=$3; } END { print a[int(i/2)]; }')
	copies=$(($depth/$genomehapdepth))
	do echo $mini $depth $copies>> $textfile
done

fi