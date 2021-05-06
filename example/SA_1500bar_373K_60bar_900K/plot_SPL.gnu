#!/usr/bin/gnuplot
#
# Colored tics with the epslatex terminal
#
reset

# set path of config snippets
set loadpath '/home/ad738098/work/AtomizationSprays_2019/PLOTTING_CONFIGS'

# epslatex
#set terminal epslatex size 15cm,8cm color colortext standalone header '\definecolor{t}{rgb}{0.0,0.0,0.0}'
set terminal epslatex size 15cm,8cm color colortext header '\definecolor{t}{rgb}{0.0,0.0,0.0}'
set output 'PL_sim_exp.tex'

# load plotting style
load 'mypalette.pal'

set xlabel 'Time / ms'
set ylabel 'Spray Penetration Length / mm' offset 2
set xrange [0.0:1.5]
set yrange [0.0:60.0]

#set key above right vertical maxrows 4
#set key left top
set key outside right

SSPL = 'spray_PL.out'

ELPL = 'Liq_SprayA_Evap.dat'
EVPL = 'Vap_SprayA_Evap.dat'


plot ELPL u 1:($2+$3):($2-$3) with filledcurves fs transparent solid 0.25 lc rgb color002 notitle,\
     EVPL u 1:($2+$3):($2-$3) with filledcurves fs transparent solid 0.25 lc rgb color001 notitle,\
     SSPL u 1:2 t 'Sim: LPL' w l ls 202, \
     ELPL u 1:2 every 3 t 'Exp: LPL' w p ls 210,\
     SSPL u 1:3 t 'Sim: VPL' w l ls 101, \
     EVPL u 1:2 every 3 t 'Exp: VPL' w p ls 120
