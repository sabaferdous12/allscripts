#!/bin/bash                                                                    
#*************************************************************************      
#                                                                               
#   Program:    processAntibodyMartin                                                     
#   File:       processAntibodyMartin.sh                                                  
#                                                                               
#   Version:    V1.1                                                            
#   Date:       22.04.14                                                        
#   Function:   This is an automatic pipeline to process the Data for AbDb and 
#               It generates data just for Martin numbering scheme 
#   Usage:      ./processAntibodyMartin.sh                                                 
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
 # parameter
free=$1
proAntigen=$2
npAntigen=$3
Combined=$4
#dir=Combined_$nsch
Redundant="Redundant_files"
mkdir -p $Combined
mkdir -p $Redundant
# copy all the complexes in Combined directory
cp ./$free/* ./$Combined 
cp ./$proAntigen/* ./$Combined
cp ./$npAntigen/* ./$Combined

cd $free
perl ~/scripts/bin/getRedundantAntibodyClusters.pl;
mv *.txt ../$Redundant
cd ..

cd $npAntigen
perl ~/scripts/bin/getRedundantAntibodyClusters.pl; 
mv *.txt ../$Redundant
cd ..

cd $proAntigen
perl ~/scripts/bin/getRedundantAntibodyClusters.pl;
mv *.txt ../$Redundant
cd ..

cd $Combined
perl ~/scripts/bin/getRedundantAntibodyClusters.pl;
mv *.txt ../$Redundant
cd ..
makeNRdata

data="Data"
mkdir -p $data
mv $free ./$data
mv $proAntigen ./$data
mv $Combined ./$data
mv $npAntigen ./$data
mv "NR_"* ./$data

#mv -f $Redundant ./$data

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

function makeNRdata
{
mkdir -p "NR_"$proAntigen
cd $proAntigen  # runs perl program for non-redundant data   
perl ~/scripts/bin/constructNonRedundantData.pl ../$Redundant/"Redundant_"$proAntigen".txt" "NR_"$proAntigen
processed_complexes=`ls|wc -l` # counts number of PDB complexes formed
cd ..

mkdir -p "NR_"$free
cd $free
perl ~/scripts/bin/constructNonRedundantData.pl ../$Redundant/"Redundant_"$free".txt" "NR_"$free
processed_antibody=`ls|wc -l` 
cd ..

mkdir -p "NR_"$npAntigen
cd $npAntigen
perl ~/scripts/bin/constructNonRedundantData.pl ../$Redundant/"Redundant_"$npAntigen".txt" "NR_"$npAntigen
processed_hapten=`ls|wc -l`
cd ..

mkdir -p "NR_"$Combined
cd $Combined
perl ~/scripts/bin/constructNonRedundantData.pl ../$Redundant/"Redundant_"$Combined.txt "NR_"$Combined
processed_combined=`ls|wc -l`
cd ..

#for directory in  *
#do
 #   dirname=$(basename "$directory")
  #  tar -jcvf $dirname.tar.bz2 $dirname
#done

#cd "NR_"$proAntigen # counts number of non redundant complexes
#NR_Complex=`ls|wc -l`
#cd ..

#cd "NR_"$free
#NR_Antibody=`ls|wc -l`
#cd ..

#cd "NR_"$npAntigen
#NR_Hapten=`ls|wc -l`
#cd ..

#cd "NR_"$Combined
#NR_Combined=`ls|wc -l`
#ls | grep "_" | cut -f1 -d_ > ../Combined.txt
#cd ../..
}





# The main program begins here

# Running get_antibody_complex program for 3 numbering schemes
echo "get_antibody_complex program is running for Martin numbering";
scheme="Martin"
free="FreeAntibody_"$scheme
proAntigen="AntibodyAntigen_"$scheme
npAntigen="AntibodyHapten_"$scheme
Combined="CombinedAb_"$scheme
perl ~/scripts/bin/processAntibodyPDBs.pl -a $1
process $free $proAntigen $npAntigen $Combined

free="LightChain_"$scheme
proAntigen="LightAntigen_"$scheme
npAntigen="LightHapten_"$scheme
Combined="CombinedLg_"$scheme
process $free $proAntigen $npAntigen $Combined

free="HeavyChain_"$scheme
proAntigen="HeavyAntigen_"$scheme
npAntigen="HeavyHapten_"$scheme
Combined="CombinedHv_"$scheme
process $free $proAntigen $npAntigen $Combined

log="Logs"
mkdir -p $log
mv *.log ./$log 
mv *.list ./$log


mv -f $Redundant ./$data
cd $data/$Redundant

# To sort all the files in directory
shopt -s nullglob
filearray=( * ) # Reading directory files into an array
 for i in "${filearray[@]}"
 do
     sort $i -o $i
 done

cd ..


# To generate non redundant data set (Only for Martin numbering) 
#mkdir -p "NR_"$proAntigen
#cd $proAntigen  # runs perl program for non-redundant data   
#perl ~/scripts/bin/constructNonRedundantData.pl ../$Redundant/"Redundant_"$proAntigen".txt" "NR_"$proAntigen
#processed_complexes=`ls|wc -l` # counts number of PDB complexes formed
#cd ..

#mkdir -p "NR_"$free
#cd $free
#perl ~/scripts/bin/constructNonRedundantData.pl ../$Redundant/"Redundant_"$free".txt" "NR_"$free
#processed_antibody=`ls|wc -l` 
#cd ..

#mkdir -p "NR_"$npAntigen
#cd $npAntigen
#perl ~/scripts/bin/constructNonRedundantData.pl ../$Redundant/"Redundant_"$npAntigen".txt" $npAntigen
#processed_hapten=`ls|wc -l`
#cd ..

#mkdir -p "NR_"$Combined
#cd $Combined
#perl ~/scripts/bin/constructNonRedundantData.pl ../$Redundant/"Redundant_"$Combined.txt "NR_"$Combined
#processed_combined=`ls|wc -l` 
#cd ..
 
#rm *.tar.bz2 # checks if there is any compressed folder already. remove them
#for directory in  *
#do
 #   dirname=$(basename "$directory")
  #  tar -jcvf $dirname.tar.bz2 $dirname
#done

#cd NR_Complex_Martin # counts number of non redundant complexes
#NR_Complex=`ls|wc -l`
#cd ..

#cd NR_Antibody_Martin
#NR_Antibody=`ls|wc -l`
#cd ..

#cd NR_Hapten_Martin
#NR_Hapten=`ls|wc -l`
#cd ..

#cd NR_Combined_Martin
#NR_Combined=`ls|wc -l`
#ls | grep "_" | cut -f1 -d_ > ../Combined.txt
#cd ../..
exit;

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

#cd ..
bash ~/scripts/bin/statsProcessed.sh $complex $processed_complexes $NR_Complex $hapten $processed_hapten $NR_Hapten $antibody $processed_antibody $NR_Antibody $completeDataset $processed_combined $NR_Combined >stats_processed.tt
bash ~/scripts/bin/statsUnprocessed.sh $bensJones $lightAntigen $camelids $heavyAntigen $fc $kabatFailed $superseded >stats_unprocessed.tt

cd $data
mkdir -p NR_Martin_merged
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
exit;

# This shell script does not update the website...
