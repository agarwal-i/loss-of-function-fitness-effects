#!/bin/bash

nsim=$1
mutationrate=$2
ufactor=$3
u_sd=$4
dom=$5
arrayid=$6
dirlabel=$7
type=$8 #fixed_sel prior neutral

mkdir -p out/${type}_${dirlabel}

Rscript get_${type}_params.R $mutationrate $dom $nsim $arrayid $ufactor $u_sd
wait
while read i ; do ./simulator_mod_sd $i ; done < $type.params.$arrayid.txt > out/${type}_${dirlabel}/$arrayid.txt
