#!/bin/bash


#for f in `ls M*/Rep*/ss/*.txt`;
for f in `ls M*/ss/*.txt`;
do
 mut=`echo $f|cut -d '/' -f1`;
 echo -n $mut",";grep "^Average" $f|awk  -F "=" '{ print $2}';
done

#ab samajh main asani se ajayega
#:)
