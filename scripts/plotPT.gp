tikzfile="PT.tex"
set term tikz standalone
set output tikzfile

# find . -regex ".*avspec\.dat\b" | sed -e 's/.*_\([^]]*\)_.*/\1/g'
# find . -regex ".*out\.avspec\.dat\b" | sed -e 's/.*_U\([^]]*\)out.*/\1/g'
files = system("echo $(ls -1 *avspec.dat)")
set xlabel "omega"
set ylabel "A(omega)"
#plot 'me_GImp_U0.out.avspec.dat' u 1:2 with lines title 
plot for [data in files] data u 1:2 with lines title data
#system(sprintf("sed -e 's/.*_U\([^]]*\)out.*/\1/g' %s", data))
#do for [file in files]{
#    plot file u 1:2 with lines title system(sprintf("sed -e 's/.*_U\([^]]*\)out.*/\1/g' %s", file)) 
#}
set output
