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
## added this:
NMPreads=$(samtools view $samfile | awk '$2==81 || $2==83 || $2==161 || $2==163 || $2 == 97 || $2 == 99 || $2 == 145 || $2 == 147 { tot ++ $2 } END { print +tot}')
NMQ20reads=$(samtools view -c -F 4 -q 20 $samfile)

NCSBreads=$(samtools view $samfile | egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC')
NCSBMreads=$(samtools view -F 4 $samfile | egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC')

## changed high quality to perfect alignments:
NCSBperfectalignments=$(samtools view $samfile | egrep 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC' | awk '{print $6}' | grep -cvE '[S/D/I/*]')

echo "N_R:" $Nreads > $out.mapping.stats.txt # changed '>' to '>>'
echo "N_MR:" $NMreads >> $out.mapping.stats.txt
echo "N_MR_PP:" $NMPreads >> $out.mapping.stats.txt
echo "N_MR_HQ:" $NMQ20reads >> $out.mapping.stats.txt
echo "N_CSB3-R:" $NCSBreads >> $out.mapping.stats.txt
echo "N_CSB3-R_M:" $NCSBMreads >> $out.mapping.stats.txt
echo "N_CSB3-R_PA:" $NCSBperfectalignments >> $out.mapping.stats.txt

# Calculating read depths
samtools sort $samfile > $out.sorted.sam
samtools depth -aa -d 0 $out.sorted.sam > $out\.depth ## changed 'a' to 'aa', changed '>' to '>>'
MinicircleContigs=$(awk '{print $1}' $out\.depth | sort | uniq)

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

echo "CONTIG" "AVERAGE.DEPTH" "MEDIAN.DEPTH" "MIN.DEPTH" "MAX.DEPTH " > $out.depth.stats.txt # changed '>' to '>>'

for mini in $MinicircleContigs
	do echo $mini $(grep "$mini" $out\.depth | awk '{print $3}' | calcs)  >> $out.depth.stats.txt
done

fi

