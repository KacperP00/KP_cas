set xrange [0:1]
set yrange [0:50]
plot 'Liq_SprayA_Evap.dat' u 1:2, 'Vap_SprayA_Evap.dat' u 1:2, 'PL_dsd_2d.out' u 1:2 w l , 'PL_dsd_2d.out' u 1:3 w l
pause 1
reread
