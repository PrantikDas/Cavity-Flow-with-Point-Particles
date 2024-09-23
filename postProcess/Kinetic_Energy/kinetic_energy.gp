# plot_script.gp

# Set the title for the plot
set title "Kinetic Energy vs Simulation Time"

# Set labels for the axes
set xlabel "Simulation Time"
set ylabel "Kinetic Energy"

# Range 
set xrange [0:20]
set yrange [0:10]


# Set grid
set grid

# Plot the data from file
plot 'Kinetic_Energy.txt' using 1:2 with lines title 'Kinetic Energy',

# Set output to a PNG file (optional, comment this block if you want to view the plot directly)
set terminal pngcairo
set output 'kinetic_energy.png'
replot
unset output
