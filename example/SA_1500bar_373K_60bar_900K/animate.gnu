set term gif animate
set output 'output.gif'

# define fixed axis-ranges
#set xrange [0.0:60]
#set yrange [0.0:200.0]
#set zrange [-1:1]

# n=number of data files
n = system(sprintf('ls data.out_*|wc -l'))

do for [j=1:n-1] {
    	filename = sprintf("data.out_%06d",(j-1)*100)
    	t = system(sprintf("awk 'NR==%d{print $2}' spray_PL.out",(2+(j-1)*10)))
    	i = system(sprintf("awk 'NR==%d{print $1}' spray_PL.out",(2+(j-1)*10)))
	print i,", ",t
    	set title 'time '.t.' ms'
    	#plot filename u ($1*90e-3):($15) w l t 'b'
    	plot filename u 1:($4) w l t 'Y_v'
    	#plot filename u 1:9 w l t 'dm', filename u 1:10 w l t 'dvar'
    	#plot filename u 1:13 w l t 'Td', filename u 1:14 w l t 'Tg'
    	#plot filename u 1:19 w l t 'Zvar_g'
}
