 #!/bin/bash
#*************************************************************************
#
#   Program:    UpdateAbDb
#   File:       UpdateAbDb.sh
#   
#   Version:    V1.1
#   Date:       22.04.14
#   Function:   This is an automatic pipeline to process the Data for AbDb and 
#               it also updates the Data directory at /acrm/www/html/abs/abdb/ 
#   Usage:      ./UpdateAbDb.sh 
#               
#
#   Copyright:  (c) UCL, Saba Ferdous, 2014
#   Author:     Miss Saba Ferdous
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College
#               Gower Street
#               London
#               WC1E 6BT
#   EMail:      saba@bioinf.org.uk
#               
#*************************************************************************


function process
{
nsch=$1 # parameter

#dir=Combined_$nsch
Redundant="Redundant_files"
mkdir -p Combined_$nsch
mkdir -p $Redundant
# copy all the complexes in Combined directory
cp ./Antibody_$nsch/* ./Combined_$nsch 
cp ./Complex_$nsch/* ./Combined_$nsch
cp ./Hapten_$nsch/* ./Combined_$nsch

cd Antibody_$nsch
perl ~/scripts/bin/getRedundantAntibodyClusters.pl;
mv *.txt ../$Redundant
cd ..

cd Hapten_$nsch
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
mv Hapten_$nsch ./$data
} # Function ends
##############################

function Split 
{
pattern=$1;
String=`grep $pattern masterlog-a.log` # e.g: returns complexes=400
parameter=$(echo $String | cut -f1 -d=) # splits the string and returns number
count=$(echo $String | cut -f2 -d=) # after =, the delimiter
}

# ===========================
# The main program begins here
# ===========================
# Running processAntibodyPDBs program for 3 numbering schemes
echo "processAntibodyPDBs program is running for Kabat numbering";
scheme="Kabat"
perl ~/scripts/bin/processAntibodyPDBs.pl -k $1
process $scheme

echo "processAntibodyPDBs program is running for Chothia numbering";
scheme="Chothia"
perl ~/scripts/bin/processAntibodyPDBs.pl -c $1
process $scheme

echo "processAntibodyPDBs.pl program is running for Martin numbering";
scheme="Martin"
perl ~/scripts/bin/processAntibodyPDBs.pl -a $1
process $scheme


log="Logs"
mkdir -p $log
mv *.log ./$log 
mv *.list ./$log
mv -f $Redundant ./$data

# Sorting the redundant clusters files
cd $data/$Redundant
sort Redundant_Complex_Martin.txt -o Redundant_Complex_Martin.txt
sort Redundant_Antibody_Martin.txt -o Redundant_Antibody_Martin.txt
sort Redundant_Hapten_Martin.txt -o Redundant_Hapten_Martin.txt
sort Redundant_Combined_Martin.txt -o Redundant_Combined_Martin.txt
sort Redundant_Complex_Kabat.txt -o Redundant_Complex_Kabat.txt
sort Redundant_Antibody_Kabat.txt -o Redundant_Antibody_Kabat.txt
sort Redundant_Hapten_Kabat.txt -o Redundant_Hapten_Kabat.txt
sort Redundant_Combined_Kabat.txt -o Redundant_Combined_Kabat.txt
sort Redundant_Complex_Chothia.txt -o Redundant_Complex_Chothia.txt
sort Redundant_Hapten_Chothia.txt -o Redundant_Hapten_Chothia.txt
sort Redundant_Antibody_Chothia.txt -o Redundant_Antibody_Chothia.txt
sort Redundant_Combined_Chothia.txt -o Redundant_Combined_Chothia.txt

cd ..

# To generate non redundant data set (Only for Martin numbering)
mkdir -p NR_Complex_Martin 
cd Complex_Martin  # runs perl program for non-redundant data   
perl ~/scripts/bin/constructNonRedundantData.pl ../$Redundant/Redundant_Complex_Martin.txt NR_Complex_Martin
processed_complexes=`ls|wc -l` # counts number of PDB complexes formed
cd ..

mkdir -p NR_Antibody_Martin
cd Antibody_Martin 
perl ~/scripts/bin/constructNonRedundantData.pl ../$Redundant/Redundant_Antibody_Martin.txt NR_Antibody_Martin
processed_antibody=`ls|wc -l` 
cd ..

mkdir -p NR_Hapten_Martin
cd Hapten_Martin
perl ~/scripts/bin/constructNonRedundantData.pl ../$Redundant/Redundant_Hapten_Martin.txt NR_Hapten_Martin
processed_hapten=`ls|wc -l`
cd ..

mkdir -p NR_Combined_Martin
cd Combined_Martin
perl ~/scripts/bin/constructNonRedundantData.pl ../$Redundant/Redundant_Combined_Martin.txt NR_Combined_Martin
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

cd NR_Hapten_Martin
NR_Hapten=`ls|wc -l`
cd ..

cd NR_Combined_Martin
NR_Combined=`ls|wc -l`
ls | grep "_" | cut -f1 -d_ > ../Combined.txt
cd ../..


# The following bit of code finds stats of processed PDBs from masterlog.log 
cd $log
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

pattern="Hapten="
Split $pattern
hapten=$count

cp Kabat_Failed.list ../$data
cp ../header.dat ../$data
cd ../$data
# To merge 2 consective lines into one
awk 'NR%2{printf $0" ";next;}1' header.dat >headerProcessed.dat
cd ..

bash ~/scripts/bin/statsProcessed.sh $complex $processed_complexes $NR_Complex $hapten $processed_hapten $NR_Hapten $antibody $processed_antibody $NR_Antibody $completeDataset $processed_combined $NR_Combined >stats_processed.tt
bash ~/scripts/bin/statsUnprocessed.sh $bensJones $lightAntigen $camelids $heavyAntigen $fc $kabatFailed $superseded >stats_unprocessed.tt

cd $data

mkdir -p  NR_Martin_merged
cd NR_Complex_Martin
cp *.pdb ../NR_Martin_merged
cd ..

cd NR_Antibody_Martin
cp *.pdb ../NR_Martin_merged
cd ..

cd NR_Hapten_Martin
cp *.pdb ../NR_Martin_merged
cd ..

# To copy the PDB codes from the directory
cd NR_Martin_merged
ls | grep "_" | cut -f1 -d_ > ../Merged.txt
cd .. 

# This shell command finds diffrence between Combined.txt and Merged.txt
comm -13 <(sort Combined.txt) <(sort Merged.txt) >difference.txt
# This script produces list of antibody PDBs that are present in both forms 
# i.e Free or Complexed
perl ~/scripts/bin/FreeComplexedAntibody.pl

cd ..

# To update the Wed Data
#=======================

webData="/acrm/www/html/abs/abdb/"
cp stats_processed.tt $webData
cp stats_unprocessed.tt $webData

cdate=$(date +"%Y-%m-%d") # Date of the Backup
# If Data folder already exists on $webData then rename with that date
if [ -d "$webData/Data" ]; then 
	mv "$webData/Data" $webData/"Data_Old_".$datestring
fi

 	cp -r ./Data $webData

cd $webData
make



