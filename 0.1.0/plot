set term pos eps enh

set logscale xy
set xlabel 'Number of steps'
set ylabel 'Maximum error (m)'
set ytics format "%g"
set xrange [8000:700000]
set yrange [1e-5:2]

set out 'max-convergence.eps'
plot 'out2-rk-1-1-0' u 2:3 t'RK1-{/Symbol D}t' w linesp pt 11,\
	'out2-rk-1-1-1' u 2:3 t'RK1-k_t' w linesp pt 10,\
	'out2-rk-2-1-0' u 2:3 t'RK2-{/Symbol D}t' w linesp pt 9,\
	'out2-rk-2-1-1' u 2:3 t'RK2-k_t' w linesp pt 8,\
	'out2-rk-3-2-0' u 2:3 t'RK3-{/Symbol D}t' w linesp pt 7,\
	'out2-rk-3-2-1' u 2:3 t'RK3-k_t' w linesp pt 6,\
	'out2-rk-4-3-0' u 2:3 t'RK4-{/Symbol D}t' w linesp pt 5,\
	'out2-rk-4-3-1' u 2:3 t'RK4-k_t' w linesp pt 4,\

set ylabel 'RMSE (m)'

set out 'rmse-convergence.eps'
plot 'out2-rk-1-1-0' u 2:4 t'RK1-{/Symbol D}t' w linesp pt 11,\
	'out2-rk-1-1-1' u 2:4 t'RK1-k_t' w linesp pt 10,\
	'out2-rk-2-1-0' u 2:4 t'RK2-{/Symbol D}t' w linesp pt 9,\
	'out2-rk-2-1-1' u 2:4 t'RK2-k_t' w linesp pt 8,\
	'out2-rk-3-2-0' u 2:4 t'RK3-{/Symbol D}t' w linesp pt 7,\
	'out2-rk-3-2-1' u 2:4 t'RK3-k_t' w linesp pt 6,\
	'out2-rk-4-3-0' u 2:4 t'RK4-{/Symbol D}t' w linesp pt 5,\
	'out2-rk-4-3-1' u 2:4 t'RK4-k_t' w linesp pt 4,\

set key bottom
set xlabel 'k_t'
unset logscale x
set xrange [0:0.6]

set out 'rmse-kt.eps'
plot 'out2-rk-1-1-1' u 1:4 t'RK1-k_t' w linesp pt 10,\
	'out2-rk-2-1-1' u 1:4 t'RK2-k_t' w linesp pt 8,\
	'out2-rk-3-2-1' u 1:4 t'RK3-k_t' w linesp pt 6,\
	'out2-rk-4-3-1' u 1:4 t'RK4-k_t' w linesp pt 4,\

set ylabel 'Maximum error (m)'

set out 'max-kt.eps'
plot 'out2-rk-1-1-1' u 1:3 t'RK1-k_t' w linesp pt 10,\
	'out2-rk-2-1-1' u 1:3 t'RK2-k_t' w linesp pt 8,\
	'out2-rk-3-2-1' u 1:3 t'RK3-k_t' w linesp pt 6,\
	'out2-rk-4-3-1' u 1:3 t'RK4-k_t' w linesp pt 4,\
