for i in 1 2 3 4; do
	sed "s/rk = 1/rk = $i/g" case2 > tmp
	sed -i "s/convergence = 10/convergence = 2/g" tmp
	sed -i "s/size = 1/size = 0/g" tmp
	sed -i "s/kt = 0\.6/dt = 0.05/g" tmp
	./ballisticpgo tmp out2-rk-$i-1-0
done
for i in 2 3 4; do
	sed "s/rk = 1/rk = $i/g" case2 > tmp
	sed -i "s/convergence = 10/convergence = 2/g" tmp
	sed -i "s/land = 1/land = 2/g" tmp
	sed -i "s/size = 1/size = 0/g" tmp
	sed -i "s/kt = 0\.6/dt = 0.05/g" tmp
	./ballisticpgo tmp out2-rk-$i-2-0
done
for i in 3 4; do
	sed "s/rk = 1/rk = $i/g" case2 > tmp
	sed -i "s/convergence = 10/convergence = 2/g" tmp
	sed -i "s/land = 1/land = 3/g" tmp
	sed -i "s/size = 1/size = 0/g" tmp
	sed -i "s/kt = 0\.6/dt = 0.05/g" tmp
	./ballisticpgo tmp out2-rk-$i-3-0
done
sed "s/size = 1/size = 0/g" case2ms > tmp
sed -i "s/convergence = 9/convergence = 2/g" tmp
sed -i "s/kt = 0\.6/dt = 0.05/g" tmp
./ballisticpgo tmp out2-ms-2-1-1-0
sed -i "s/rk = 1/rk = 2/g" tmp
./ballisticpgo tmp out2-ms-2-2-1-0
sed -i "s/land = 1/land = 2/g" tmp
./ballisticpgo tmp out2-ms-2-2-2-0
sed -i "s/rk = 2/rk = 1/g" tmp
./ballisticpgo tmp out2-ms-2-1-2-0
sed -i "s/steps = 2/steps = 3/g" tmp
./ballisticpgo tmp out2-ms-3-1-2-0
sed -i "s/rk = 1/rk = 2/g" tmp
./ballisticpgo tmp out2-ms-3-2-2-0
sed -i "s/land = 2/land = 1/g" tmp
./ballisticpgo tmp out2-ms-3-2-1-0
sed -i "s/rk = 2/rk = 1/g" tmp
./ballisticpgo tmp out2-ms-3-1-1-0
for i in 1 2 3 4; do
	sed "s/rk = 1/rk = $i/g" case2 > tmp
	sed -i "s/convergence = 10/convergence = 2/g" tmp
	./ballisticpgo tmp out2-rk-$i-1-1
done
for i in 2 3 4; do
	sed "s/rk = 1/rk = $i/g" case2 > tmp
	sed -i "s/convergence = 10/convergence = 2/g" tmp
	sed -i "s/land = 1/land = 2/g" tmp
	./ballisticpgo tmp out2-rk-$i-2-1
done
for i in 3 4; do
	sed "s/rk = 1/rk = $i/g" case2 > tmp
	sed -i "s/convergence = 10/convergence = 2/g" tmp
	sed -i "s/land = 1/land = 3/g" tmp
	./ballisticpgo tmp out2-rk-$i-3-1
done
cp case2ms tmp
sed -i "s/convergence = 9/convergence = 2/g" tmp
./ballisticpgo tmp out2-ms-2-1-1-1
sed -i "s/rk = 1/rk = 2/g" tmp
./ballisticpgo tmp out2-ms-2-2-1-1
sed -i "s/land = 1/land = 2/g" tmp
./ballisticpgo tmp out2-ms-2-2-2-1
sed -i "s/rk = 2/rk = 1/g" tmp
./ballisticpgo tmp out2-ms-2-1-2-1
sed -i "s/steps = 2/steps = 3/g" tmp
./ballisticpgo tmp out2-ms-3-1-2-1
sed -i "s/rk = 1/rk = 2/g" tmp
./ballisticpgo tmp out2-ms-3-2-2-1
sed -i "s/land = 2/land = 1/g" tmp
./ballisticpgo tmp out2-ms-3-2-1-1
sed -i "s/rk = 2/rk = 1/g" tmp
./ballisticpgo tmp out2-ms-3-1-1-1
for i in 2 3; do
	sed "s/rk = 1/rk = $i/g" case2 > tmp
	sed -i "s/convergence = 10/convergence = 2/g" tmp
	sed -i "s/size = 1/size = 0/g" tmp
	sed -i "s/kt = 0\.6/dt = 0.05/g" tmp
	sed -i "s/error-dt = 0/error-dt = 1/g" tmp
	./ballisticpgo tmp out2-rk-$i-1-2
done
for i in 2 3; do
	sed "s/rk = 1/rk = $i/g" case2 > tmp
	sed -i "s/convergence = 10/convergence = 2/g" tmp
	sed -i "s/size = 1/size = 0/g" tmp
	sed -i "s/kt = 0\.6/dt = 0.05/g" tmp
	sed -i "s/land = 1/land = 2/g" tmp
	sed -i "s/error-dt = 0/error-dt = 1/g" tmp
	./ballisticpgo tmp out2-rk-$i-2-2
done
for i in 3; do
	sed "s/rk = 1/rk = $i/g" case2 > tmp
	sed -i "s/convergence = 10/convergence = 2/g" tmp
	sed -i "s/size = 1/size = 0/g" tmp
	sed -i "s/kt = 0\.6/dt = 0.05/g" tmp
	sed -i "s/land = 1/land = 3/g" tmp
	sed -i "s/error-dt = 0/error-dt = 1/g" tmp
	./ballisticpgo tmp out2-rk-$i-3-2
done
rm tmp
gnuplot plot