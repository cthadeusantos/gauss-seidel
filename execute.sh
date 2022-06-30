clear
count=5
for i in $(seq $count); do
    mpirun -np 2 ./mpi.out 10 matrix/matrix16.mtx
done