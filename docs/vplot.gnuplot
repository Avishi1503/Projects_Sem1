# Run from the Algorithm_Comparison folder after generating the long matrix.
# python src/vplot.py --input data/mapped.bed.gz --long --output results/vplot/vplot.long.tsv
set terminal pngcairo size 1400,850 enhanced font 'Arial,12'
set output 'results/vplot/vplot.gnuplot.png'
set title 'V-plot: fragment midpoint relative to region midpoint'
set xlabel 'Offset (bp)'
set ylabel 'Fragment length (bp)'
set cblabel 'Region / fragment pair count'
set xrange [-500:499]
set yrange [0:700]
set view map
set size ratio 0.6
unset key
set palette defined (0 '#fcfdbf', 1 '#fca50a', 2 '#dd513a', 3 '#932667', 4 '#420a68', 5 '#000004')
plot 'results/vplot/vplot.long.tsv' using 1:2:3 with image
