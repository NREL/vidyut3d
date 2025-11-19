#main code==================================================
command1=sprintf("ls %s",ARG1)
filenames = system(command1)
command2=sprintf("ls -1 %s | wc -l",ARG1)
N_str=system(command2)
N=int(N_str)
colnum=int(ARG2)

set palette defined (1 "blue", N "red")
#other fancy palettes
#set palette rgbformulae 7,5,15
#set palette defined (0 "red", 1 "orange", 2 "yellow", 3 "green", 4 "blue", 5 "violet")

plot for [i=1:N] word(filenames,i) using 1:colnum:($1*0+i) with lines lw 1.5 lc palette notitle

#==================================================
#trying other fancy things, may be useful in future
#==================================================

#just learning to create arrays
#array pltid[N]
#do for [i=1:N] {
#    pltid[i] = i  
#}

#trying other fancy things
#doesnt work as it updates the plot to the latest file
#do for [i=1:words(filenames)] {
#        current_file = word(filenames, i)
#        plot current_file using 1:colnum:($1*0+i) with lines lc palette notitle
#    }

#another way with arrays
#plot for [i=1:N] word(filenames,i) using 1:colnum:($3*0+pltid[i]) with lines lc palette notitle

#another way
#plot for [i=0:N-1] sprintf("plt%5.5d.slice", i) using 1:colnum:($3*0+i) with lines lc palette notitle
#================================================
