set term postscript eps color enhanced
set output "fit_press.eps"
f(x)=a
fit f(x) "press.dat" u 1:2 via a

set xlabel "time step"
set ylabel "Pressure (bar)"


set label 1 sprintf("P=%6.1f",a/1e4) at graph 0.5, graph 0.6 font ",18"
#set label 1 sprintf("A=%3.4f, a=%3.4f",A, a) at graph 0.5, graph 0.6 font ",18"
#set label 2 sprintf("b1=%3.4f, b2=%3.4f",b1,b2) at graph 0.5, graph 0.5 font ",18"
#set label 3 sprintf("eta=%3.4f",A*a*b1+A*(1-a)*b2) at graph 0.5, graph 0.4 font ",18"

p "press.dat" u 1:2 ti 'data' w l, f(x) ti "fit"
