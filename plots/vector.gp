set terminal postscript eps enhanced size 8,9 font "cmr10,20"
set output "vector.eps"
set multiplot layout 3,2
set tics nomirror
set yrange [-0.45:1.1]
set xlabel "t"
set ylabel "{/Symbol \341}v^{(i)}{/Symbol \361}"
set key font ",15"
set title "{/Symbol n} = 10^{-2}"
plot "../results/vector_0.dat" u 1:6 w l lc 0 lw 1.5 notitle, "" u 1:3 every 200 pt 4 lc 0 ti "x", "" u 1:7 w l lc 0 lw 1.5 notitle, "" u 1:4 every 200 pt 6 lc 0 ti "y", "" u 1:8 w l lc 0 lw 1.5 notitle, "" u 1:5 every 200 pt 8 lc 0 ti "z"
set title "{/Symbol n} = 10^{-1}"
plot "../results/vector_1.dat" u 1:6 w l lc 0 lw 1.5 notitle, "" u 1:3 every 200 pt 4 lc 0 notitle, "" u 1:7 w l lc 0 lw 1.5 notitle, "" u 1:4 every 200 pt 6 lc 0 notitle, "" u 1:8 w l lc 0 lw 1.5 notitle, "" u 1:5 every 200 pt 8 lc 0 notitle
set title "{/Symbol n} = 0.5"
plot "../results/vector_2.dat" u 1:6 w l lc 0 lw 1.5 notitle, "" u 1:3 every 200 pt 4 lc 0 notitle, "" u 1:7 w l lc 0 lw 1.5 notitle, "" u 1:4 every 200 pt 6 lc 0 notitle, "" u 1:8 w l lc 0 lw 1.5 notitle, "" u 1:5 every 200 pt 8 lc 0 notitle
set title "{/Symbol n} = 1"
plot "../results/vector_3.dat" u 1:6 w l lc 0 lw 1.5 notitle, "" u 1:3 every 200 pt 4 lc 0 notitle, "" u 1:7 w l lc 0 lw 1.5 notitle, "" u 1:4 every 200 pt 6 lc 0 notitle, "" u 1:8 w l lc 0 lw 1.5 notitle, "" u 1:5 every 200 pt 8 lc 0 notitle
set title "{/Symbol n} = 10^{1}"
plot "../results/vector_4.dat" u 1:6 w l lc 0 lw 1.5 notitle, "" u 1:3 every 200 pt 4 lc 0 notitle, "" u 1:7 w l lc 0 lw 1.5 notitle, "" u 1:4 every 200 pt 6 lc 0 notitle, "" u 1:8 w l lc 0 lw 1.5 notitle, "" u 1:5 every 200 pt 8 lc 0 notitle
set title "{/Symbol n} = 10^{2}"
plot "../results/vector_5.dat" u 1:6 w l lc 0 lw 1.5 notitle, "" u 1:3 every 200 pt 4 lc 0 notitle, "" u 1:7 w l lc 0 lw 1.5 notitle, "" u 1:4 every 200 pt 6 lc 0 notitle, "" u 1:8 w l lc 0 lw 1.5 notitle, "" u 1:5 every 200 pt 8 lc 0 notitle
unset multiplot
set output