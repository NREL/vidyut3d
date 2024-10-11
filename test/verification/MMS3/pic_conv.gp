set terminal pdf
set output "conv.pdf"
set termoption font "Helvetica,18"
set key spacing 1.1
#set lmargin 12.5
#set rmargin 3
#set bmargin 3
set xlabel "time step"
set ylabel "Error L2 norm"
set key bottom right
#set xtics format "%5.2e"
set ytics format "%5.2e"
set log
#set yrange [1e-7:10]
#set xrange [0.005:0.3]
plot 'err1' u 1:2 w lp  lw 2 ps 1 pt 5 lc 6 title "nsteps=1",\
    'err2' u 1:2 w lp  lw 3 pt 6 lc 7 title "nsteps=2",\
    'err4' u 1:2 w lp  lw 3 pt 7 lc 5 title "nsteps=4",\
    'err1' u 1:($1*5e5) w l  lw 3 lc 6 title "order=1",\
    'err1' u 1:($1*$1*5e5) w l  lw 3 lc 7 title "order=2"
