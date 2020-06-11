#!/bin/bash

if [ $# -eq 0 ]; then
    echo ""
    echo "FVDB"
    echo ""
    echo "Calculate mapping statistics from SAM/BAM files"
    echo ""
    echo "Usage: mapping_stats.sh inputfile"
    echo ""
else

samfile=$1
textfile=$2

Nreads=$(samtools view -c $samfile)
NMreads=$(samtools view -c -F 4 $samfile)
NMQ20reads=$(samtools view -c -F 4 -q 20 $samfile)
NMPPreads=$(samtools view -c -f 3 $samfile)

NCSBreads=$(samtools view $samfile | egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC')
NCSBMreads=$(samtools view -F 4 $samfile | egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC')
NCSBMQ20reads=$(samtools view -F 4 -q 20 $samfile | egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC')
NCSBPPreads=$(samtools view -f 3 $samfile | egrep -c 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC')

echo "Number of reads:" $Nreads >> $textfile
echo "Number of mapped reads:" $NMreads >> $textfile
echo "Number of reads with mapping quality >= 20:" $NMQ20reads >> $textfile
echo "Number of properly paired reads:" $NMPPreads >> $textfile
echo "Number of CSB3-containing reads:" $NCSBreads >> $textfile
echo "Number of mapped CSB3-containing reads:" $NCSBMreads >> $textfile
echo "Number of CSB3-containing reads with mapping quality >= 20:" $NCSBMQ20reads >> $textfile
echo "Number of properly paired CSB3-containing reads:" $NCSBPPreads >> $textfile

fi