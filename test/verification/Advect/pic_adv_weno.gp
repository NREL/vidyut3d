set terminal png
set output "conv_weno.png"
set termoption font "Helvetica,15"
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
set ytics 1e4
set yrange [1e-12:1]
set xrange [0.005:0.3]
plot 'errconv5' u 1:2 w lp  lw 4 ps 1 pt 5 lc 6 title "WENO-Z",\
    'errconv5' u 1:($1**5) w l  lw 4 title "5th order"
