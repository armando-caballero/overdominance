#!/bin/bash
#$ -cwd

rm script_totalvclass_SD.sh.*

if [ $# -ne 5 ]  
then
	echo "Usage: $0 <PS> <LOD> <sod1> <sod2> <n> " 
	exit 1
fi

#Set arguments
PS=$1
LOD=$2
sod1=$3
sod2=$4
n=$5

WDIR=$PWD

if [ -d "RESULTS$n" ]
then
rm -r RESULTS$n
fi
mkdir -p $WDIR/RESULTS$n

mkdir -p /state/partition1/armandoPURGA$n/$SLURM_JOBID/

cp naturalvclass /state/partition1/armandoPURGA$n/$SLURM_JOBID/
cp purgingvclass /state/partition1/armandoPURGA$n/$SLURM_JOBID/

touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`

cd /state/partition1/armandoPURGA$n/$SLURM_JOBID

num=$RANDOM
echo "$num" > seedfile

time ./naturalvclass>>out<<@
0
-99
1000	N
$PS	PS(99=random)
99	Lenght genome (99=free)
200	NCRO (max 2000)(Neu=Ncro)
30	NLOCI (2-30)
0.2	Lambda_a
0.015	Lambda_L
$LOD	Lambda_OD
$sod1	sod1
$sod2	sod2
0.0	absolute effect of lethal (QT): normal (aL,aL)
0	random proportion(0) or large mutants (1)
1.0	Psi
0.33	beta_s
0.33	beta_a
0.2	ave |s|
0.2	ave |a|
0.0	PP_s
0.0	PP_a
2	dom model (0=cnt; 1:Deng, 2:CK94 gamma)
0	h_s (mod 0), k_s (mod 1)
0.283	ave h_s (mod 2)
0	h_a (mod 0), k_a (mod 1)
0.283	ave h_a (mod 2)
99	rho (99:a=s)
0	Vs
1	multi(1), add(2)
10000	generations
2000	gen/block
@

time ./purgingvclass>>out<<@
0
-99
50	N
0.0	A
$PS	PS
10	J (Family size)
99	L (99: free recombination)
200	NCRO (max 500)(Neu=Ncro)
30	NLOCI (2-30)
0.2	Lambda_s
0.015	Lambda_L
$LOD	Lambda_OD
$sod1	sod1
$sod2	sod2
0.33	beta_s
0.2	ave |s|
0.283	ave h_s
250	generations
10	lines
20	replicates
@

cp -r /state/partition1/armandoPURGA$n/$SLURM_JOBID/genfile $WDIR/RESULTS$n/
cp -r /state/partition1/armandoPURGA$n/$SLURM_JOBID/genfile.dat $WDIR/RESULTS$n/
cp -r /state/partition1/armandoPURGA$n/$SLURM_JOBID/dfilename $WDIR/RESULTS$n/
cp -r /state/partition1/armandoPURGA$n/$SLURM_JOBID/popfile $WDIR/RESULTS$n/
cp -r /state/partition1/armandoPURGA$n/$SLURM_JOBID/datafile $WDIR/RESULTS$n/

#Una vez copiado todo lo necesario hemos de limpiar el scratch eliminando el subdirectorio creado por el trabajo:
rm -r /state/partition1/armandoPURGA$n/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*
