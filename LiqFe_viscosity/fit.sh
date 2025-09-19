#!/bin/bash


#sys_size="3 4 6 8 10"
#sys_size="3 4 6"
#sys_size="8"

RHOS="10.7 11.0 11.3 11.5 12.0 12.5 13.0 13.3"
TEMPS="3000 3500 4300 4800 5000 5500 6000 6500 7000"

PREFIX="WORKDIR"
#PREFIX="WORKDIR-8x8x8"
SUBDIR=$(pwd)

#for size in $sys_size
for rho in $RHOS
do
for temp in $TEMPS
do
   WORKDIR=$SUBDIR/${PREFIX}/T${temp}-rho${rho}
   #mkdir -p ${WORKDIR}
   cd ${WORKDIR}
   pwd
   #size=8
   #cp ${SUBDIR}/create_volume/rho${rho}-${size}${size}${size}.lmp pos.lmp

   #cp ${SUBDIR}/in.lmp in.lammps.orig
   #cp ${SUBDIR}/in.lmp in.lammps
   #sed 's/T equal 6000.0/T equal '"${temp}/1" in.lammps.orig > in.lammps
	

   #cp ${SUBDIR}/lmp.sh .
   #cp ${SUBDIR}/graph.pb .
   
   #sbatch lmp.sh
   #cat > ther.sh <<EOF
#~/work/pkg/gk/build/bin/gk --task ther --input current.dat -d 400000 -w 1000 -s 500 --tc 300 -v 0 -c 8 
#EOF
#   cat > visc.sh <<EOF
#~/work/pkg/gk/build/bin/gk --task visc --input stress.dat -d 40000 -w 1000 -s 800 --tc 150 -v 0
#EOF
   #bash ther.sh
 #  bash visc.sh
 gkplot.sh
 #gnuplot plt.gpl
   
   cat > fit.sh <<EOF
#!/bin/bash

visfile="eta.dat"

cat > fit_viscosity.gnu << END
set term postscript eps color enhanced
set output "fit_viscosity.eps"
f(x)=A*a*b1*(1-exp(-x/b1)) +  A*(1-a)*b2*(1-exp(-x/b2)) 
A=0.2; a=0.5; b1=0.05; b2=0.2;
fit f(x) "\$visfile" u 1:5 via A,a,b1,b2

# original one exponential
#f(x)=a*(1-exp(-x/b))
#a=0.07; b=10;
#fit f(x) "\$visfile" u 1:5 via a,b


set xlabel "time (ps)"
set ylabel "{/Symbol h} (Pa s)"

#set label 1 sprintf("a = %3.4f",a) at graph 0.5, graph 0.6 font ",18"
#set label 2 sprintf("b = %3.4f",b) at graph 0.5, graph 0.5 font ",18"

set label 1 sprintf("A=%3.4f, a=%3.4f",A, a) at graph 0.5, graph 0.6 font ",18"
set label 2 sprintf("b1=%3.4f, b2=%3.4f",b1,b2) at graph 0.5, graph 0.5 font ",18"
set label 3 sprintf("eta=%3.4f",A*a*b1+A*(1-a)*b2) at graph 0.5, graph 0.4 font ",18"

p "\$visfile" u 1:5 ti 'data' w l, f(x) ti "fit"
END
   gnuplot fit_viscosity.gnu
EOF
   rm fit.log
   bash fit.sh
   mv fit.log fit_visc.log


#############################################################
# diffusion coefficient
   cat > fit_msd.sh <<EOF
#!/bin/bash

visfile="msd.dat"

cat > fit_diff.gnu << END
set term postscript eps color enhanced
set output "fit_diff.eps"
#f(x)=a*(1-exp(-x/b))
f(x)=1e-4*6.0*D*x
D=0.1;
fit f(x) "\$visfile" u 1:5 via D

set xlabel "time (ps)"
set ylabel "D (10^{-9} m^2/s)"

set label 1 sprintf("D = %3.4f",D) at graph 0.5, graph 0.6 font ",18"
#set label 2 sprintf("b = %3.4f",b) at graph 0.5, graph 0.5 font ",18"

p "\$visfile" u 1:5 ti 'data' w l, f(x) ti "fit"
END
   gnuplot fit_diff.gnu
EOF
   rm fit.log
   bash fit_msd.sh
   cp fit.log fit_diff.log

   cd ${SUBDIR}
done
done


