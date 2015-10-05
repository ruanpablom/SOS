reset
set terminal pngcairo
#set autoscale
#set lmargin 8
#set bmargin 4
#set key 13000,340
set output "speedUp.png"
set grid
set title 'SOS Threads'
set xlabel 'Thread'
set ylabel 'Time Sequential/Time Thread(x)'
plot "aceleracao.txt" using ($1):($2) title "Speedup" w lines lw 4 lt 1 lc 2, \
	 "aceleracao.txt" using ($1):($2) title "Thread Speed" w points pt 7 ps 2, \
            "aceleracao.txt" using ($1):($3) title "Linear Speedup" w lines lw 4 lt 10 lc 5, \
	        "aceleracao.txt" using ($1):($3) title "Linear Thread Speed" w points pt 7 ps 2; 
