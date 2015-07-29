#!/bin/bash

# Make html page

#echo $1 $2 $3 $4 $5 $6 $7 $8 $9 ' -> echo $1 $2 $3 $4 $5 $6 $7 $8 $9'

args=("$@")


cat << __EOF__
<div class='abdb_stats'>
<h2>Database Statistics</h2>

<table class="mytable">
<tr>
<th>Antibody Structures</th>
<th>Total Processed PDBs </th>
<th>Total Resultant Complexes </th>
<th>Non-Redundant Antibodies</th>
</tr>


<tr>
<td>Complexes</td>
<td>${args[0]}</td>
<td>${args[1]}</td>
<td>${args[2]}</td>
</tr>

<tr>
<td>Free Antibody</td>
<td>${args[3]}</td>
<td>${args[4]}</td>
<td>${args[5]}</td>
</tr>

<tr>
<td>Complete Dataset</td>
<td>${args[6]}</td>
<td>${args[7]}</td>
<td>${args[8]}</td>
</tr>

</table>

</div>

__EOF__
