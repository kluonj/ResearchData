#!/bin/bash

#TEMPLIST="4000 4500 5000 5500 6000 6500 7000 7500 8000"
TEMPLIST="5000 5500 6000 6500 7000 7500 8000"
#RHOLIST="0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4"
#RHOLIST="0.4  0.5  0.6  0.7 0.8 0.9 1.0" 
RHOLIST="1.1 1.2 1.3 1.4"

SUBDIR=$(pwd)

for TEMP in $TEMPLIST
do
  for RHO in $RHOLIST
  do
    WORKDIR=$SUBDIR/RUN/$RHO-$TEMP
    mkdir -p $WORKDIR

    cd  $WORKDIR

    echo "Inside $WORKDIR"
    cat > job.sh  <<EOF
#!/bin/bash -l
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --partition amd_a8_384

source ~/.bashrc
conda activate ~/soft/deepmd-kit/


ulimit -s unlimited
export OMP_NUM_THREADS=1


mpirun -n \$SLURM_NTASKS lmp -in in.lammps -var rho ${RHO} \
	-var natoms 4096 -var coa 1.0 -var nsteps 100000 \
	-var highT ${TEMP} -var lowT ${TEMP}
EOF

	cp $SUBDIR/in.lammps .
  
	#ln -snf ~/kluo/liq-vap/run/long_train/graph.pb pair.pb
	ln -snf $SUBDIR/graph-compress.pb pair.pb

	sbatch job.sh 

	echo "Slurm job submitted!"
done
done
