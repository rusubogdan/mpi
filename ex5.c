// Exercise 5 - http://www.hpc.cineca.it/content/exercise-5
// Distribute a global square NxN matrix (with N fixed) over P processors, so that each task has a local portion of
// it in its memory. Initialize such portion with the task's rank.
// Each task sends its first and last columns (if Fortran) or rows (if C) to its neighbours (i.e.: the 
// first column/row to the left processor, the last column/row to the right processor). 
// Note that each task has to actually allocate a larger portion of its submatrix (ghost cells), with two extra
// columns/rows at the extremities to hold the received values. 

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> 

int main(int argc, char *argv[])  
{

	if (argc < 3) return -1;

	MPI_Init(&argc, &argv);

	int rank, size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int** rows;
	int** matrix;
	int n = atoi(argv[1]);
	int tag = 200;

	printf("%d\n", n);

	// max rows to be distributed over processes 
		int rowsSize = n / size;
		if (n % size > 0)
			rowsSize++;

	// dynamic partial matrix allocation
	rows = malloc(rowsSize * sizeof(int*));
	for (i = 0; i < rowsSize; i++) {
		rows[i] = malloc(n * sizeof(int));
	}

	if (rank == 0) {
		// dynamic matrix allocation
		matrix = malloc(n * sizeof(int*));
		int i;
		for (i = 0; i < n; i++) {
			matrix[i] = malloc(n * sizeof(int));
		}	

		// first part is kept in master
		// build and distribute the partial matrices
		int rankIt = 1; // start from node 1
		for (i = rowsSize; i < n;) {
			// send starting from line i of size rowSize
			MPI_Send(&matrix[i], rowsSize * n, MPI_Int, rankIt, tag, MPI_COMM_WORLD);
		}


	} else {

		MPI_Recv(&rows, rowsSize * n, MPI_Int, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
	}

	return 0;
}