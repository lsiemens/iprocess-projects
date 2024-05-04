# Plot data from spring simple simulation
set title "Simple Spring Simulation"
set xlabel "Time t"
set multiplot layout 1,2 columns
plot "spring_simple.dat" using 1:2 index 0 title "x" with lines, "" using 1:4 index 0 title "energy" with lines
plot "spring_simple.dat" using 1:2 index 1 title "x" with lines, "" using 1:4 index 1 title "energy" with lines
unset multiplot