# plot_script.gp

# Set the title for the plot
set title "Residual vs Simulation Time"

# Set labels for the axes
set xlabel "Simulation Time"
set ylabel "Residual"

# Range 
set xrange [0:20]
# set yrange [-100:100]

set logscale y 10

# Set grid
set grid

# Plot the data from file
plot 'Residuals.txt' using 1:2 with lines title 'u', \
     'Residuals.txt' using 1:3 with lines title 'v', \
     'Residuals.txt' using 1:4 with lines title 'p'

# Set output to a PNG file (optional, comment this block if you want to view the plot directly)
set terminal pngcairo
set output 'Residual.png'
replot
unset output
