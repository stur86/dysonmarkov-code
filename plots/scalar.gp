set terminal postscript eps enhanced size 8,9 font "cmr10,20"
set output "scalar.eps"
set multiplot layout 3,2
set tics nomirror
set yrange [-1:1]
set xlabel "t"
set ylabel "S"
set title "{/Symbol n} = 10^{-2}"
plot "../results/scalar_0.dat" u 1:4 w l lc 0 lw 1.5 notitle, "" u 1:3 every 100 pt 2 lc 0 notitle
set title "{/Symbol n} = 10^{-1}"
plot "../results/scalar_1.dat" u 1:4 w l lc 0 lw 1.5 notitle, "" u 1:3 every 100 pt 2 lc 0 notitle
set title "{/Symbol n} = 1"
plot "../results/scalar_2.dat" u 1:4 w l lc 0 lw 1.5 notitle, "" u 1:3 every 100 pt 2 lc 0 notitle
set title "{/Symbol n} = 5"
plot "../results/scalar_3.dat" u 1:4 w l lc 0 lw 1.5 notitle, "" u 1:3 every 100 pt 2 lc 0 notitle
set title "{/Symbol n} = 10^{1}"
plot "../results/scalar_4.dat" u 1:4 w l lc 0 lw 1.5 notitle, "" u 1:3 every 100 pt 2 lc 0 notitle
set title "{/Symbol n} = 10^{2}"
plot "../results/scalar_5.dat" u 1:4 w l lc 0 lw 1.5 notitle, "" u 1:3 every 100 pt 2 lc 0 notitle
unset multiplot
set output