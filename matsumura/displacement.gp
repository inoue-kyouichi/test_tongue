set terminal postscript eps enhanced color "Helvetica" 35 size 4in, 4in

set border 255 lw 2
set size 2.5,1.5
#set size square
set size ratio 0.5

set xlabel "{Muscle Fiber}" font "Helvetica,40"
set ylabel "{Displacement} [mm]" font "Helvetica,40" offset 0,0

#set style data lines

#set xrange [0:21]
#set yrange [0.0:0.15]

set bmargin 4
set tmargin 2
#set lmargin 4
#set rmargin 4

#set tics scale 3 ,2 font "helvetica-bold,35"
#set tics scale 3 ,1 font "Helvetica,40"
#set xtics  (0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
#set ytics  (0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15) 
#set mytics 0.5
#set logscale x
#set logscale y
#set format y "10^{%L}"

#set key left top
#set key at graph 1.35, graph 0.96
#set key right top
#set key at 1000.0 ,67.0
set key spacing 1.2
#set key box 

set style line 1 lt 1 lc rgb "#a9a9a9" lw 6 ps 4 pt 12 
set boxwidth 1.5 absolute
set style fill solid


set key at graph 1.05, graph 1.15
set output "displacement.eps"
p\
 "displacement.dat" u 1:2 w boxes  title "Number of bonding area : 8" ls 1


