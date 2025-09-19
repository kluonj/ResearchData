#!/bin/bash

visfile="eta.dat"

cat > fit_viscosity.gnu << END
set term postscript eps color enhanced
set output "fit_viscosity.eps"

f(x)=A*a*b1*(1-exp(-x/b1)) +  A*(1-a)*b2*(1-exp(-x/b2)) 
A=0.2; a=0.80; b1=0.01; b2=0.02;
fit f(x) "$visfile" u 1:5 via A,a,b1,b2

#f(x)=a*(1-exp(-x/b))
#a=0.01; b=1;
#fit f(x) "$visfile" u 1:5 via a,b

set xlabel "time (ps)"
set ylabel "{/Symbol h} (Pa s)"

#set label 1 sprintf("a = %3.4f",a) at graph 0.5, graph 0.6 font ",18"
#set label 2 sprintf("b = %3.4f",b) at graph 0.5, graph 0.5 font ",18"

set label 1 sprintf("A=%3.4f, a=%3.4f",A, a) at graph 0.5, graph 0.6 font ",18"
set label 2 sprintf("b1=%3.4f, b2=%3.4f",b1,b2) at graph 0.5, graph 0.5 font ",18"
set label 3 sprintf("eta=%3.4f",A*a*b1+A*(1-a)*b2) at graph 0.5, graph 0.4 font ",18"


p "$visfile" u 1:5 ti 'data' w l, f(x) ti "fit"
END
   gnuplot fit_viscosity.gnu
