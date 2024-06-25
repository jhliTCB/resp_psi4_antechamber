#!/bin/bash
sb() {
cat <<_EOF
#!/bin/bash
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem 150GB
#SBATCH -t 32:00:00
#SABTCH -w comput$2
#SBATCH -p cpu_part
#SBATCH -J $1

script=/public/home/lijunhao/scripts/resp_psi4_antechamber
inpPDB=$1
scratch=$HOME/scratch
echo \$inpPDB
check() { [[ \$1 != 0 ]] && echo "\$2 not correctly finished!" && exit 1 || echo "\$2 Done!" ; }
chg=1
\$script/psi4_esp.py $3 -ipdb \${inpPDB}.pdb -cm "\$chg,1" -scr \$scratch -c \$SLURM_CPUS_PER_TASK -m \$[SLURM_MEM_PER_NODE/1024]
check \$? psi4_esp.py
#\$script/fit_cov_resp.py \${inpPDB}_opt.pdb results_\${inpPDB}.esp
#check \$? fit_cov_resp.py
\$script/gesp2pdb.py results_\${inpPDB}.esp \${inpPDB}.pdb
check \$? gesp2pdb.py
_EOF
}
i=0
cs=(5 5 8 8)
for f in cmpd45.pdb  cmpd56.pdb  cmpd65.pdb  MT17141.pdb; do
  mkdir -p noopt_${f%.*} && cd noopt_${f%.*}
  ln -s ../$f .
  sb ${cs[$i]} noopt_${f%.*} "-no_opt" > ttt.sh
  [[ -z $1 ]] && sbatch ttt.sh
  ((i+=1))
  cd ..
done
  
