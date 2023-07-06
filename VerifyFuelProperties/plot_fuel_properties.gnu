#!/usr/bin/gnuplot
#
# Colored tics with the epslatex terminal
#
reset

# set path of config snippets
set loadpath '/home/ad738098/work/AtomizationSprays_2019/PLOTTING_CONFIGS'

# epslatex
set terminal epslatex size 15cm,7cm color colortext standalone header '\definecolor{t}{rgb}{0.0,0.0,0.0}'
#set terminal epslatex size 15cm,7cm color colortext header '\definecolor{t}{rgb}{0.0,0.0,0.0}'
set output 'HeatOfVaporization.tex'

# load plotting style
load 'mypalette.pal'

set xlabel 'Temperature / K'
set ylabel 'Heat of Vaporization / J/kg/K' offset 2
set xrange [300.0:700.0]
#set yrange [0.0:60.0]

set key right
#set key width -8
#set key outside above 

DBE = './dinbutylether/LiquidFuelProperties.dat'
OCT = './1octanol/LiquidFuelProperties.dat'

plot OCT u 1:7 t '1octanol' w l ls 201, \
     DBE u 1:7 t 'dinbutylether' w l ls 305, \
