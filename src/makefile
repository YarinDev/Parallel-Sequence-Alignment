build:
	mpicxx -c main.c
	mpicxx -o project main.c 

clear:
	rm -rf *.o project

run:
	mpiexec -n 4 ./project < input.txt
	
run2:
	mpiexec -n 4 ./project < input2.txt
