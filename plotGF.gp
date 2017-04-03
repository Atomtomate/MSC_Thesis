set terminal svg enhanced size 2000,1000
set object 1 rect from screen 0, 0, 0 to screen 1, 1, 0 behind
set object 1 rect fc  rgb "white"  fillstyle solid 1.0
set key autotitle columnhead
set samples 110




set output 'G0_Guess.svg'
set title 'G_{0,up}({/Symbol t})'
plot '0_G0_Guess_IT.out' using 1:3 ps 0.2  title 'Weiss function initial guess'

set output 'GImp.svg'
set title 'G_{Imp,up}(i {/Symbol w}_n)'
filename(n) = sprintf("%d_GImp_MF.out", n)
plot for [i=0:30] filename(i) using 1:3 ps 0.2  title sprintf("Iteration %d",i)

set output 'GImp_IT.svg'
set title 'G_{Imp,up}(i {/Symbol t}_n)'
filename(n) = sprintf("%d_GImp_IT.out", n)
plot for [i=0:30] filename(i) using 1:2 ps 0.2  title sprintf("Iteration %d",i)


set output 'G0u.svg'
set title 'G_{0,up}(\tau)'
filename(n) = sprintf("%d_G0_IT.out", n)
plot for [i=0:30] filename(i) using 1:2 ps 0.2 title sprintf("Iteration %d",i)

set output 'G0u_mf.svg'
set title 'G_{0,up}(tau)'
filename(n) = sprintf("%d_G0_MF.out", n)
plot for [i=0:20] filename(i) using 1:3 ps 0.2 title sprintf("Iteration %d",i)

set output 'G0d.svg'
set title 'G_{0,down}(\tau)'
filename(n) = sprintf("%d_G0_IT.out", n)
plot for [i=0:20] filename(i) using 1:3 ps 0.2 title sprintf("Iteration %d",i)

set output 'G0d_mf.svg'
set title 'G_{0,up}(tau)'
filename(n) = sprintf("%d_G0_MF.out", n)
plot for [i=0:20] filename(i) using 1:5 ps 0.2 title sprintf("Iteration %d",i)

set output 'SImp.svg'
set title 'Sigma_{Imp,up}(i {/Symbol w}_n)'
filename(n) = sprintf("%d_sImp_MF.out", n)
plot for [i=0:20] filename(i) using 1:3 ps 0.2  title sprintf("Iteration %d",i)

set output 'GLoc_Re.svg'
set title 'Re[G_{loc,up}(i {/Symbol w}_n)]'
filename(n) = sprintf("%d_GLoc_MF.out", n)
plot for [i=0:20] filename(i) using 1:2 with lines smooth csplines title sprintf("Iteration %d",i)

set output 'GLoc_Im.svg'
set title 'Im[G_{loc,up}(i {/Symbol w}_n)]'
filename(n) = sprintf("%d_GLoc_MF.out", n)
plot for [i=0:20] filename(i) using 1:3 ps 0.2 title sprintf("Iteration %d",i)

set output 'Delta.svg'
set title 'Delta_{up}(i {/Symbol w}_n)'
filename(n) = sprintf("%d_Delta_MF.out", n)
plot for [i=0:20] filename(i) using 1:3 ps 0.2  title sprintf("Iteration %d",i)

set output 'G0_guess_IT.svg'
set title 'G_{0,up}({/Symbol t}_n)'
plot "0_G0_Guess_IT.out" using 1:2 ps 0.2
set output 'G0_guess_MF.svg'
set title 'G_{0,up}(i {/Symbol w}_n)'
plot "0_G0_Guess_MF.out" using 1:3 ps 0.2
