set terminal pngcairo dashed size 800,500
#set terminal pdf
set output "compare_eden.png"
set termoption font "Helvetica,12"
set key spacing 1.2
set lmargin 9
set rmargin 5
set bmargin 3.5
set xlabel "non-dim dist"
set ylabel "Electron Density (#/m3)"
set yrange [1e13:3.5e15]
plot 'linedata_x0032.dat' u ($1/0.067):5 w l lw 3 title "Vidyut3d",\
 'somafoam_soln' u 1:($2*1e14) w p pt 4 ps 1 lw 3 title "Verma and Venkattraman, Comp. Phys. Comm., 2021",\
 'Turner_Momsoln' u 1:($2*1e15) w p pt 4 ps 1 lw 3 title "Turner et al., POP, 2013 (Moment)",\
 'Turner_PICsoln' u 1:2 w p pt 4 ps 1 lw 3 title "Turner et al., POP, 2013 (PIC)"
