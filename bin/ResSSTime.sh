#!/bin/bash

for f in `ls M*/Rep*/ss/*_ss.txt`;
do
echo $f;mut=`echo $f|cut -d '/' -f1`;
l=`grep -v "^Average" $f|sed 's/ *//g;'|cut -d ':' -f2`;
for ll in `echo $l`;
do tempfile=`basename $f .txt`.txt;
echo "$ll,$mut" >> $tempfile;
done;
paste $tempfile temporary > $tempfile.pasted;rm $tempfile;done

for f in `ls M*/Rep*/ss/*_ss.txt`;do echo $f;done > listorder.txt

for f in `cat listorder.txt`;
do
 cat `basename $f`.pasted >> 2lrw9e_ss_avg_per_resi.csv;
done

perl -pi -e 's/[^\S\n]/,/g' 2lrw9e_ss_avg_per_resi.csv;
