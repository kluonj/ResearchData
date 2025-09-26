#!/bin/bash

TEMPLIST=$(seq 3500 250 6000)
#TEMPLIST=5750
RHOLIST="0.6"

SUBDIR=$(pwd)

for TEMP in $TEMPLIST
do
  for RHO in $RHOLIST
  do
    WORKDIR=$SUBDIR/PROD/$RHO-$TEMP
    WORKDIROLD=$SUBDIR/RUN/$RHO-$TEMP
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


mpirun -n \$SLURM_NTASKS lmp -in re.lammps -var rho ${RHO} \
	-var natoms 4000 -var coa 5.0 -var nsteps 50000 \
	-var highT ${TEMP} -var lowT ${TEMP}

conda activate ase
python make_gif.py
python plot_z_density.py
python gmm_analysis_filter.py
EOF

	cp $SUBDIR/make_gif.py .
	cp $SUBDIR/plot_z_density.py .
	cp $SUBDIR/gmm_analysis_filter.py .
	cp $SUBDIR/re.lammps .

	ln -snf $SUBDIR/graph.pb pair.pb
	#ln -snf /public1/home/a8s000527/kluo/liq-vap-production/dpgen/long_train/rcut6/graph-compress.pb pair.pb
  
	cp $WORKDIROLD/*.data .
	#ln -snf ~/kluo/liq-vap/run/long_train/graph.pb pair.pb
	#ln -snf ~/kluo/liq-vap/run/rcut6/graph-compress.pb pair.pb

	sbatch job.sh 

	echo "Slurm job submitted!"
done
done
