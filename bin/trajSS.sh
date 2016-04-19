#!/bin/bash

name=$1;
# Make lists of 1000 frames (pdbsecstr output .ss files)
ls -tr -1 *.ss | split -l 1000 -d - lists
# Merging all .ss files column-wise in each list (list00..10)
for list in lists*;
do
    paste -d "" $(cat $list) > merge${list##lists};
done

# Column-wise merging of each list (merge00..10)
paste -d "" merge* >$name"_temp.dat"
rm lists*
rm merge*