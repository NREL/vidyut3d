set terminal pdf
set output "conv_phi.pdf"
set termoption font "Helvetica,18"
set key spacing 1.1
#set lmargin 12.5
#set rmargin 3
#set bmargin 3
set xlabel "Grid size"
set ylabel "Error L2 norm"
set key bottom right
#set xtics format "%5.2e"
set ytics format "%5.2e"
set log
#set yrange [1e-7:10]
set xrange [0.005:0.3]
plot 'err1' u 1:3 w lp  lw 2 ps 1 pt 5 lc 6 title "hyporder=1",\
    'err2' u 1:3 w lp  lw 3 pt 6 lc 7 title "hyporder=2",\
    'err5' u 1:3 w lp  lw 3 pt 7 lc 5 title "hyporder=5",\
    'err1' u 1:($1/50) w l  lw 3 lc 6 title "",\
    'err1' u 1:($1*$1/100) w l  lw 3 lc 7 title ""
