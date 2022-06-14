clear
count=5
for i in $(seq $count); do
    ./mp_linux.out 10 matrix/matrix10000.mtx 4
done