set terminal postscript eps enhanced dashed size 8,9 font "cmr10,20"
set output "quantum_corr.eps"
set multiplot layout 3,2
set tics nomirror
set yrange [-0.1:1.1]
set xlabel "t"
set ylabel "{/Symbol \341}{/:Bold I}_1{/:Bold I}_2{/Symbol \361}"
set key font ",15"
set title "{/Symbol n} = 10^{-2}"
plot "../results/quantum_0.dat" u 1:3 w l lc 0 lw 1.5 notitle, "" u 1:4 every 100 pt 2 lc 0 notitle, 0 dt 2 lc 0 notitle
set title "{/Symbol n} = 10^{-1}"
plot "../results/quantum_1.dat" u 1:3 w l lc 0 lw 1.5 notitle, "" u 1:4 every 100 pt 2 lc 0 notitle, 0 dt 2 lc 0 notitle
set title "{/Symbol n} = 0.5"
plot "../results/quantum_2.dat" u 1:3 w l lc 0 lw 1.5 notitle, "" u 1:4 every 100 pt 2 lc 0 notitle, 0 dt 2 lc 0 notitle
set title "{/Symbol n} = 1"
plot "../results/quantum_3.dat" u 1:3 w l lc 0 lw 1.5 notitle, "" u 1:4 every 100 pt 2 lc 0 notitle, 0 dt 2 lc 0 notitle
set title "{/Symbol n} = 10^{1}"
plot "../results/quantum_4.dat" u 1:3 w l lc 0 lw 1.5 notitle, "" u 1:4 every 100 pt 2 lc 0 notitle, 0 dt 2 lc 0 notitle
set title "{/Symbol n} = 10^{2}"
plot "../results/quantum_5.dat" u 1:3 w l lc 0 lw 1.5 notitle, "" u 1:4 every 100 pt 2 lc 0 notitle, 0 dt 2 lc 0 notitle
unset multiplot
set output "quantum_entropy.eps"
set multiplot layout 3,2
set tics nomirror
set yrange [0:1.6]
set ylabel "S"
set key font ",15"
set title "{/Symbol n} = 10^{-2}"
plot "../results/quantum_0.dat" u 1:5 w l lc 0 lw 1.5 notitle, "" u 1:6 every 100 pt 2 lc 0 notitle, log(4) dt 2 lc 0 notitle
set title "{/Symbol n} = 10^{-1}"
plot "../results/quantum_1.dat" u 1:5 w l lc 0 lw 1.5 notitle, "" u 1:6 every 100 pt 2 lc 0 notitle, log(4) dt 2 lc 0 notitle
set title "{/Symbol n} = 0.5"
plot "../results/quantum_2.dat" u 1:5 w l lc 0 lw 1.5 notitle, "" u 1:6 every 100 pt 2 lc 0 notitle, log(4) dt 2 lc 0 notitle
set title "{/Symbol n} = 1"
plot "../results/quantum_3.dat" u 1:5 w l lc 0 lw 1.5 notitle, "" u 1:6 every 100 pt 2 lc 0 notitle, log(4) dt 2 lc 0 notitle
set title "{/Symbol n} = 10^{1}"
plot "../results/quantum_4.dat" u 1:5 w l lc 0 lw 1.5 notitle, "" u 1:6 every 100 pt 2 lc 0 notitle, log(4) dt 2 lc 0 notitle
set title "{/Symbol n} = 10^{2}"
plot "../results/quantum_5.dat" u 1:5 w l lc 0 lw 1.5 notitle, "" u 1:6 every 100 pt 2 lc 0 notitle, log(4) dt 2 lc 0 notitle
unset multiplot
set output