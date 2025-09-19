set term postscript eps color enhanced
set output "fit_diff.eps"
#f(x)=a*(1-exp(-x/b))
f(x)=1e-4*6.0*D*x
D=0.1;
fit f(x) "msd.dat" u 1:5 via D

set xlabel "time (ps)"
set ylabel "D (10^{-9} m^2/s)"

set label 1 sprintf("D = %3.4f",D) at graph 0.5, graph 0.6 font ",18"
#set label 2 sprintf("b = %3.4f",b) at graph 0.5, graph 0.5 font ",18"

p "msd.dat" u 1:5 ti 'data' w l, f(x) ti "fit"
