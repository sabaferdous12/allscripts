#!/bin/sh
### Usage: 
## This script runs in the directory containing the mutants and
## creates a directory ss and creates a file (.dat) for SS assignments
## for each frame of the trajectory by calling the script <trajSSHandling.pl>

### Then it calls <getTimepercentageofSS.pl> to calculate the percentage of time 
### a peptide stayed in its initial conformation
echo -n "enter simulation time in ps:"
read timeps
echo -n "enter start position:"
read start
echo -n "enter lenght of peptide:"
read length



for f in `ls -d M*`;
do
	cd $f
	mkdir -p ss
	cp input_1000.tpr ss/
	cp *.xtc ss/
	cd ss
	traj=`ls *.xtc`
	echo "running trajSSHandling on $m"
	~/allscript/bin/trajSSHandling.pl -s input_1000.tpr -t $traj <<< $'1\n' 
	traj_ss_dat=`ls *.dat`
	echo "running getTimepercentageofSS on $traj_ss_dat"
	~/allscript/bin/getTimePercentageOfSS.pl -n 10 -ls $timeps -s $start -lp $length -f $traj_ss_dat > `basename $traj .xtc`.txt
	cd ../..
done
	
