#set title "Cubic Hermite basis functions"
set nokey
set xlabel "$\xi$"
set label 1 "$\chbfn{1}{0}{\xi}$" at 0.30, 0.648,0 centre
set label 2 "$\chbfn{1}{1}{\xi}$" at 0.50, 0.200,0 centre
set label 3 "$\chbfn{2}{0}{\xi}$" at 0.70, 0.648,0 centre
set label 4 "$\chbfn{2}{1}{\xi}$" at 0.65,-0.080,0 centre
set label 5 "slope=1" at 0.08,0.30 left
set label 6 "slope=1" at 0.92,0.30 right
set arrow 1 from 0.15,0.24 to 0.15, 0.15
set arrow 2 from 0.85,0.24 to 0.85,-0.15
set arrow 3 from 0.00, 0.00 to 0.18,0.18 nohead
set arrow 4 from 0.82,-0.18 to 1.00,0.00 nohead
#set xtics  0.00,0.25,1
#set ytics -0.25,0.25,1
psi10(x)=1-3*x*x+2*x*x*x
psi20(x)=x*x*(3-2*x)
psi11(x)=x*(x-1)*(x-1)
psi21(x)=x*x*(x-1)
plot[0:1] psi10(x),psi20(x),psi11(x),psi21(x)
