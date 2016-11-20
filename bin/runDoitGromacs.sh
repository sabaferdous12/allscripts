#!/bin/sh
echo -n "Enter Name:"
read nameM
for f in `ls -d Rep*`;
do
    cd $f
#    IFS='_' read -a fn <<<"$f"
#    name=${fn[0]}"_"${fn[1]}
    name=$nameM"_"$f
    ~/git/doitGROMACS_v2.1/doitGROMACS.sh -g -b acrm -n $name -t 500 -s input_500.tpr -f $name"_500.xtc" <<< 'h20'
    ~/git/doitGROMACS_v2.1/doitGROMACS.sh -g -b acrm -n $name -t 500 -s input_500.tpr -f $name"_500.xtc" <<< 'cond'
   ~/git/doitGROMACS_v2.1/doitGROMACS.sh -g -b acrm -n $name -t 500 -s input_500.tpr -f $name"_500.xtc" <<< 'rmsdfg'
   ~/git/doitGROMACS_v2.1/doitGROMACS.sh -g -b acrm -n $name -t 500 -s input_500.tpr -f $name"_500.xtc" <<< 'ggplot'

   mkdir -p ss
   cp input_500.tpr ss/
   cp *.xtc ss/
   cd ss
   traj=`ls *.xtc`
   echo "running trajSSHandling on $m"
   ~/allscript/bin/trajSSHandling.pl -s input_500.tpr -t $traj <<< $'1\n'

   cd ..   
   cd ..
done