#!/bin/bash
#final stage resp fitting
ml purge
ml amber/a22t23
[[ $# < 2 ]] && echo $0 prefix rn && exit
prefix=$1
esp=$2
pdb=$3
#sed -i "/CONECT/d" $prefix.pdb
antechamber -fi gesp -fo ac -i $esp -o ${prefix}_resp.ac -c resp -ge ${prefix}.esp -pf y
antechamber -fi ac -i ${prefix}_resp.ac -c wc -cf ${prefix}.crg -pf y
antechamber -fi pdb -i $pdb -c rc -cf ${prefix}.crg -fo ac -o ${prefix}_resp_pdb.ac -pf y
#atomtype -i ${prefix}_resp_pdb.ac -o ${prefix}_resp_pdb_gaff.ac -p gaff
prepgen -i ${prefix}_resp_pdb.ac -o ${prefix}.prepc -f car -rn LIG
parmchk2 -i ${prefix}.prepc -o ${prefix}.frcmod -f prepc

