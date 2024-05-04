# Plot data from Euler Spring Simulation
set xlabel "Time t"
set multiplot layout 1,2 columns title "Euler Spring Simulation"
set title "Short duration"
plot "euler_spring.dat" using 1:2 index 0 title "x" with lines, "" using 1:4 index 0 title "energy" with lines
set title "Long duration"
plot "euler_spring.dat" using 1:2 index 1 title "x" with lines, "" using 1:4 index 1 title "energy" with lines
unset multiplot
