sequential: 
	clear
	gcc -std=c99  gseidel_sequential.c -o sequential.out -lm

schedule:
	clear
	gcc  gseidel_mp_schedule.c -std=c99 -lm -fopenmp -o openmp_schedule.out  
