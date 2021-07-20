#!/bin/bash

if [ $# -eq 0 ]; then
    echo ""
    echo "FVDB"
    echo ""
    echo "Calculate mapping and depth statistics from SAM/BAM files"
    echo ""
    echo "Usage: mapping_stats.sh samfile"
    echo ""
else

samfile=$1
out=$(echo $samfile | sed 's/\.sam//g')

Nreads=$(samtools view -c $samfile)
NMreads=$(samtools view -c -F 4 $samfile)
NMQ20reads=$(samtools view -c -F 4 -q 20 $samfile)

NCSBreads=$(samtools view $samfile | egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC')
NCSBMreads=$(samtools view -F 4 $samfile | egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC')
NCSBMQ20reads=$(samtools view -F 4 -q 20 $samfile | egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC')

echo "Number of reads:" $Nreads >> $out.mapping.stats.txt
echo "Number of mapped reads:" $NMreads >> $out.mapping.stats.txt
echo "Number of reads with mapping quality >= 20:" $NMQ20reads >> $out.mapping.stats.txt
echo "Number of CSB3-containing reads:" $NCSBreads >> $out.mapping.stats.txt
echo "Number of mapped CSB3-containing reads:" $NCSBMreads >> $out.mapping.stats.txt
echo "Number of CSB3-containing reads with mapping quality >= 20:" $NCSBMQ20reads >> $out.mapping.stats.txt


# Calculating read depths
samtools sort $samfile > $samfile.sorted.sam
samtools depth -a $samfile.sorted.sam > $samfile\.depth
MinicircleContigs=$(awk '{print $1}' $samfile\.depth | sort | uniq)

calcs () {
sort -n | awk '
  BEGIN {
    c = 0;
    sum = 0;
  }
  $1 ~ /^(\-)?[0-9]*(\.[0-9]*)?$/ {
    a[c++] = $1;
    sum += $1;
  }
  END {
    ave = sum / c;
    if( (c % 2) == 1 ) {
      median = a[ int(c/2) ];
    } else {
      median = ( a[c/2] + a[c/2-1] ) / 2;
    }
    OFS="\t";
    print ave, median, a[0], a[c-1];
  }
'
}

echo "CONTIG" "AVERAGE.DEPTH" "MEDIAN.DEPTH" "MIN.DEPTH" "MAX.DEPTH ">> $out.depth.stats.txt

for mini in $MinicircleContigs
	do echo $mini $(grep "$mini" $samfile\.depth | awk '{print $3}' | calcs)  >> $out.depth.stats.txt
done

fi

