#!/bin/bash

TEMPLIST=$(seq 3500 250 6000)
RHOLIST="0.6"

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
#SBATCH --ntasks-per-node 8
#SBATCH --partition amd_a8_384

source ~/.bashrc
conda activate ~/soft/deepmd-kit/


ulimit -s unlimited
export OMP_NUM_THREADS=1


mpirun -n \$SLURM_NTASKS lmp -in in.lammps -var rho ${RHO} \
	-var natoms 4000 -var coa 5.0 -var nsteps 500000 \
	-var highT 12000 -var lowT ${TEMP}

conda activate ase

python make_gif.py
python plot_z_density.py
EOF

	cp $SUBDIR/make_gif.py .
	cp $SUBDIR/plot_z_density.py .
	cp $SUBDIR/in.lammps .
	ln -snf $SUBDIR/graph.pb pair.pb
  
	#ln -snf ~/kluo/liq-vap/run/long_train/graph.pb pair.pb
	#ln -snf ~/kluo/liq-vap/run/rcut6/graph-compress.pb pair.pb
	#ln -snf ~/kluo/liq-vap-production/dpgen/long_train/rcut6-high-learning/graph-compress.pb pair.pb

	sbatch job.sh 
	echo "Slurm job submitted!"
done
done
