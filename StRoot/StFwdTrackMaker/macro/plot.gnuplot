set terminal png size 800,600
set output 'memory_plot.png'


set title "Memory Usage Over Time"
set xlabel "Time (s)"
set ylabel "Memory (MB)"
set grid

plot "last.dat" using 1:2 with lines title "Memory (MB)"