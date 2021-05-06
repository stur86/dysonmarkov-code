set terminal postscript eps enhanced dashed size 6,4.5 font "cmr10,20"
set output "benchmark.eps"
set tics nomirror
set key font ",15"
set log
set xlabel "Run time (ms)"
set ylabel "{/Symbol \341}{/TimesNewRoman \174}S_{dy}(t)-S_{mc}(t){/TimesNewRoman \174}/{/Symbol s}_{mc}(t){/Symbol \361}"
set pointsize 2
set yrange [0.02:1]
set xrange [2:300]
plot "../results/benchmark_0.dat" u ($1*1e3):2 w lp pt 2 lc 0 ti "{/Symbol n} = 10^{-2}", "../results/benchmark_1.dat" u ($1*1e3):2 w lp pt 4 lc 0 ti "{/Symbol n} = 10^{-1}", "../results/benchmark_2.dat" u ($1*1e3):2 w lp pt 4 lc 0 ti "{/Symbol n} = 1", "../results/benchmark_3.dat" u ($1*1e3):2 w lp pt 6 lc 0 ti "{/Symbol n} = 5", "../results/benchmark_4.dat" u ($1*1e3):2 w lp pt 8 lc 0 ti "{/Symbol n} = 10^{1}", "../results/benchmark_5.dat" u ($1*1e3):2 w lp pt 10 lc 0 ti "{/Symbol n} = 10^{2}"
set output