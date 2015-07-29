#!/bin/bash


#############################
# Description: Performs  processing on directories 
#              Creates Combined directory for antibody and complexes generated
#              by program get_antibody_complex.pl. It creates directory for 
#             redundant file. Then it runs program get-redund.pl for redundancy#              test and stores results in directory Redundant_files. Finally 
#              moves all the directories in to the main directory i.e Data 
#############################
function process
{
nsch=$1 # parameter

#dir=Combined_$nsch
Redundant="Redundant_files"
mkdir -p Combined_$nsch
mkdir -p $Redundant
cp ./Antibody_$nsch/* ./Combined_$nsch
cp ./Complex_$nsch/* ./Combined_$nsch
cd Antibody_$nsch
perl ~/scripts/bin/getRedundantAntibodyClusters.pl;
mv *.txt ../$Redundant
cd ..
cd Complex_$nsch
perl ~/scripts/bin/getRedundantAntibodyClusters.pl; 
mv *.txt ../$Redundant
cd ..
cd Combined_$nsch
perl ~/scripts/bin/getRedundantAntibodyClusters.pl;
mv *.txt ../$Redundant
cd ..
data="Data"
mkdir -p $data
mv Antibody_$nsch ./$data
mv Complex_$nsch ./$data
mv Combined_$nsch ./$data
} # Function ends
##############################

function Split 
{
pattern=$1;
String=`grep $pattern masterlog-a.log` # e.g: returns complexes=400
parameter=$(echo $String | cut -f1 -d=) # splits the string and returns number
count=$(echo $String | cut -f2 -d=) # after =, the delimiter
}
##############################

# The main program begins here

# Running get_antibody_complex program for 3 numbering schemes
echo "get_antibody_complex program is running for Kabat numbering";
scheme="Kabat"
perl ~/scripts/bin/get_antibody_complex.pl -k SACS_updated.txt
process $scheme

echo "get_antibody_complex program is running for Chothia numbering";
scheme="Chothia"
perl ~/scripts/bin/get_antibody_complex.pl -c SACS_updated.txt 
process $scheme

echo "get_antibody_complex program is running for Martin numbering";
scheme="Martin"
perl ~/scripts/bin/get_antibody_complex.pl -a SACS_updated.txt
process $scheme


log="Logs"
mkdir -p $log
mv *.log ./$log 
mv -f $Redundant ./$data

cd $data/$Redundant
sort Redundant_Complex_Martin.txt -o Redundant_Complex_Martin.txt
sort Redundant_Antibody_Martin.txt -o Redundant_Antibody_Martin.txt
sort Redundant_Combined_Martin.txt -o Redundant_Combined_Martin.txt
sort Redundant_Complex_Kabat.txt -o Redundant_Complex_Kabat.txt
sort Redundant_Antibody_Kabat.txt -o Redundant_Antibody_Kabat.txt
sort Redundant_Combined_Kabat.txt -o Redundant_Combined_Kabat.txt
sort Redundant_Complex_Chothia.txt -o Redundant_Complex_Chothia.txt
sort Redundant_Antibody_Chothia.txt -o Redundant_Antibody_Chothia.txt
sort Redundant_Combined_Chothia.txt -o Redundant_Combined_Chothia.txt
cd ..

# To generate non redundant data set (Only for Martin numbering)
mkdir -p NR_Complex_Martin 
cd Complex_Martin  # runs perl program for non-redundant data   
perl ~/scripts/bin/copy_non_redundant.pl ../$Redundant/Redundant_Complex_Martin.txt NR_Complex_Martin
processed_complexes=`ls|wc -l` # counts number of PDB complexes formed
cd ..

mkdir -p NR_Antibody_Martin
cd Antibody_Martin 
perl ~/scripts/bin/copy_non_redundant.pl ../$Redundant/Redundant_Antibody_Martin.txt NR_Antibody_Martin
processed_antibody=`ls|wc -l` 
cd ..

mkdir -p NR_Combined_Martin
cd Combined_Martin
perl ~/scripts/bin/copy_non_redundant.pl ../$Redundant/Redundant_Combined_Martin.txt NR_Combined_Martin
processed_combined=`ls|wc -l` 
cd ..
 
#rm *.tar.bz2 # checks if there is any compressed folder already. remove them
for directory in  *
do
    dirname=$(basename "$directory")
    tar -jcvf $dirname.tar.bz2 $dirname
done

cd NR_Complex_Martin # counts number of non redundant complexes
NR_Complex=`ls|wc -l`
cd ..

cd NR_Antibody_Martin
NR_Antibody=`ls|wc -l`
cd ..

cd NR_Combined_Martin
NR_Combined=`ls|wc -l`
cd ../..
mv *.list ./$log
# The following bit of code finds stats of processed PDBs from masterlog.log 

cd Logs

pattern="Complexes="
Split $pattern
complex=$count

pattern="Free_Antibody="
Split $pattern
antibody=$count

pattern="Complete_Dataset="
Split $pattern
completeDataset=$count

pattern="Light-Antigen="
Split $pattern
lightAntigen=$count

pattern="Heavy-Antigen="
Split $pattern
heavyAntigen=$count

pattern="Bens-Jones="
Split $pattern
bensJones=$count

pattern="Camelids="
Split $pattern
camelids=$count

pattern="Fc="
Split $pattern
fc=$count

pattern="Failed="
Split $pattern
kabatFailed=$count

pattern="Superseded="
Split $pattern
superseded=$count
cp Kabat_Failed.list ../Data
#webData="/acrm/www/html/abs/abdb/"

cd ..

bash ~/scripts/bin/statsProcessed.sh $complex $processed_complexes $NR_Complex $antibody $processed_antibody $NR_Antibody $completeDataset $processed_combined $NR_Combined >stats_processed.tt
bash ~/scripts/bin/statsUnprocessed.sh $bensJones $lightAntigen $camelids $heavyAntigen $fc $kabatFailed $superseded >stats_unprocessed.tt

#cp stats_processed.tt $webData
#cp stats_unprocessed.tt $webData

cd Data

mkdir -p  NR_Martin_merged
cd NR_Complex_Martin
cp *.pdb ../NR_Martin_merged
cd ..

cd NR_Antibody_Martin
cp *.pdb ../NR_Martin_merged
cd ..

## This short perl script reads the 2 directories and then find the difference between the two. 
perl ~/scripts/bin/findDifference.pl ./NR_Combined_Martin/ ./NR_Martin_merged/
comm -13 <(sort Combined.txt) <(sort Merged.txt) >difference.txt
perl ~/scripts/bin/getFormattedList.pl

cd ..
exit;
# This shell script does not update the website...

 
