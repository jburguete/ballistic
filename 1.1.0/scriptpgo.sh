for i in 1 2 3 4; do
	sed "s/type=\"1\"/type=\"$i\"/g" ../tests/case2.xml > tmp.xml
	sed -i "s/convergence=\"10\"/convergence=\"2\"/g" tmp.xml
	sed -i "s/time-step=\"1\"/time-step=\"0\"/g" tmp.xml
	sed -i "s/kt=\"0.6\"/dt=\"0.05\"/g" tmp.xml
	./ballisticpgo tmp.xml out2-rk-$i-1-0
done
for i in 2 3 4; do
	sed "s/type=\"1\"/type=\"$i\"/g" ../tests/case2.xml > tmp.xml
	sed -i "s/convergence=\"10\"/convergence=\"2\"/g" tmp.xml
	sed -i "s/land=\"1\"/land=\"2\"/g" tmp.xml
	sed -i "s/time-step=\"1\"/time-step=\"0\"/g" tmp.xml
	sed -i "s/kt=\"0.6\"/dt=\"0.05\"/g" tmp.xml
	./ballisticpgo tmp.xml out2-rk-$i-2-0
done
for i in 3 4; do
	sed "s/type=\"1\"/type=\"$i\"/g" ../tests/case2.xml > tmp.xml
	sed -i "s/convergence=\"10\"/convergence=\"2\"/g" tmp.xml
	sed -i "s/land=\"1\"/land=\"3\"/g" tmp.xml
	sed -i "s/time-step=\"1\"/time-step=\"0\"/g" tmp.xml
	sed -i "s/kt=\"0.6\"/dt=\"0.05\"/g" tmp.xml
	./ballisticpgo tmp.xml out2-rk-$i-3-0
done
sed "s/time-step=\"1\"/time-step=\"0\"/g" ../tests/case2ms.xml > tmp.xml
sed -i "s/convergence=\"9\"/convergence=\"2\"/g" tmp.xml
sed -i "s/kt=\"0.6\"/dt=\"0.05\"/g" tmp.xml
./ballisticpgo tmp.xml out2-ms-2-1-1-0
sed -i "s/type=\"1\"/type=\"2\"/g" tmp.xml
./ballisticpgo tmp.xml out2-ms-2-2-1-0
sed -i "s/land=\"1\"/land=\"2\"/g" tmp.xml
./ballisticpgo tmp.xml out2-ms-2-2-2-0
sed -i "s/type=\"2\"/type=\"1\"/g" tmp.xml
./ballisticpgo tmp.xml out2-ms-2-1-2-0
sed -i "s/type=\"1\"/type=\"2\"/g" tmp.xml
./ballisticpgo tmp.xml out2-ms-3-1-2-0
sed -i "s/type=\"1\"/type=\"2\"/g" tmp.xml
./ballisticpgo tmp.xml out2-ms-3-2-2-0
sed -i "s/land=\"2\"/land=\"1\"/g" tmp.xml
./ballisticpgo tmp.xml out2-ms-3-2-1-0
sed -i "s/type=\"2\"/type=\"1\"/g" tmp.xml
./ballisticpgo tmp.xml out2-ms-3-1-1-0
for i in 1 2 3 4; do
	sed "s/type=\"1\"/type=\"$i\"/g" ../tests/case2.xml > tmp.xml
	sed -i "s/convergence=\"10\"/convergence=\"2\"/g" tmp.xml
	./ballisticpgo tmp.xml out2-rk-$i-1-1
done
for i in 2 3 4; do
	sed "s/type=\"1\"/type=\"$i\"/g" ../tests/case2.xml > tmp.xml
	sed -i "s/convergence=\"10\"/convergence=\"2\"/g" tmp.xml
	sed -i "s/land=\"1\"/land=\"2\"/g" tmp.xml
	./ballisticpgo tmp.xml out2-rk-$i-2-1
done
for i in 3 4; do
	sed "s/type=\"1\"/type=\"$i\"/g" ../tests/case2.xml > tmp.xml
	sed -i "s/convergence=\"10\"/convergence=\"2\"/g" tmp.xml
	sed -i "s/land=\"1\"/land=\"3\"/g" tmp.xml
	./ballisticpgo tmp.xml out2-rk-$i-3-1
done
cp ../tests/case2ms.xml tmp.xml
sed -i "s/convergence=\"9\"/convergence=\"2\"/g" tmp.xml
./ballisticpgo tmp.xml out2-ms-2-1-1-1
sed -i "s/type=\"1\"/type=\"2\"/g" tmp.xml
./ballisticpgo tmp.xml out2-ms-2-2-1-1
sed -i "s/land=\"1\"/land=\"2\"/g" tmp.xml
./ballisticpgo tmp.xml out2-ms-2-2-2-1
sed -i "s/type=\"2\"/type=\"1\"/g" tmp.xml
./ballisticpgo tmp.xml out2-ms-2-1-2-1
sed -i "s/type=\"1\"/type=\"2\"/g" tmp.xml
./ballisticpgo tmp.xml out2-ms-3-1-2-1
sed -i "s/type=\"1\"/type=\"2\"/g" tmp.xml
./ballisticpgo tmp.xml out2-ms-3-2-2-1
sed -i "s/land=\"2\"/land=\"1\"/g" tmp.xml
./ballisticpgo tmp.xml out2-ms-3-2-1-1
sed -i "s/type=\"2\"/type=\"1\"/g" tmp.xml
./ballisticpgo tmp.xml out2-ms-3-1-1-1
for i in 2 3; do
	sed "s/type=\"1\"/type=\"$i\"/g" ../tests/case2.xml > tmp.xml
	sed -i "s/convergence=\"10\"/convergence=\"2\"/g" tmp.xml
	sed -i "s/time-step=\"1\"/time-step=\"0\"/g" tmp.xml
	sed -i "s/kt=\"0.6\"/dt=\"0.05\"/g" tmp.xml
	sed -i "s/time-step=\"0\"/time-step=\"1\"/g" tmp.xml
	./ballisticpgo tmp.xml out2-rk-$i-1-2
done
for i in 2 3; do
	sed "s/type=\"1\"/type=\"$i\"/g" ../tests/case2.xml > tmp.xml
	sed -i "s/convergence=\"10\"/convergence=\"2\"/g" tmp.xml
	sed -i "s/time-step=\"1\"/time-step=\"0\"/g" tmp.xml
	sed -i "s/kt=\"0.6\"/dt=\"0.05\"/g" tmp.xml
	sed -i "s/land=\"1\"/land=\"2\"/g" tmp.xml
	sed -i "s/time-step=\"0\"/time-step=\"1\"/g" tmp.xml
	./ballisticpgo tmp.xml out2-rk-$i-2-2
done
for i in 3; do
	sed "s/type=\"1\"/type=\"$i\"/g" ../tests/case2.xml > tmp.xml
	sed -i "s/convergence=\"10\"/convergence=\"2\"/g" tmp.xml
	sed -i "s/time-step=\"1\"/time-step=\"0\"/g" tmp.xml
	sed -i "s/kt=\"0.6\"/dt=\"0.05\"/g" tmp.xml
	sed -i "s/land=\"1\"/land=\"3\"/g" tmp.xml
	sed -i "s/time-step=\"0\"/time-step=\"1\"/g" tmp.xml
	./ballisticpgo tmp.xml out2-rk-$i-3-2
done
#rm tmp.xml
gnuplot plot
