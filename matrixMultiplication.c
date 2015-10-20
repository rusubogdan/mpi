#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
 
int main(int argc, char** argv)
{
	int rank, size;
	MPI_Request *requestList,requestNull;
	MPI_Status  status;

	// Start MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
    int m, n, p;

	if (argc < 4 ) return -1;

	m = atoi(argv[1]);
	n = atoi(argv[2]);    
	p = atoi(argv[3]);

	int b[n][p]; // visible to all processes

	int slaves = size - 1;
	int maxSize;

	if( rank == 0 )
	{
		int a[m][n], c[m][p]; // visible only to the master process

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				a[i][j] = rand() % 3 + 1; 				
			}
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < p; j++) {
				b[i][j] = rand() % 3 + 1; 				
			}
		}

		// calculate fragment size of the matrix
		// x slaves = x fragments of max size N/x 
		// rows per slave = N / x = maxSize
		maxSize = m / slaves;
		if (m % slaves > 0) {
			maxSize++;
		}

		// bcast the matrix b - used by all nodes in multiplication
		MPI_Bcast(&b, n * p, MPI_INT, 0, MPI_COMM_WORLD);

		printf("Max size of each block sent to a slave %d \n", maxSize);
		
		// send size of aux matrix to all nodes
		MPI_Bcast(&maxSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

		requestList = (MPI_Request*)malloc((size-1)*sizeof(MPI_Request));

		int aux[maxSize][n];
		int d[maxSize][n];

		// extract rows from a and send them to slaves 
		for (int i = 1; i < size && rank * (maxSize-1) + i <= m; i++) {
			int z = 0;
			// create aux matrix
			for (int j = maxSize * (i - 1); j < maxSize * i && j < m; j++) {
			    for (int k = 0; k < n; k++) {
					aux[z][k] = a[j][k];
				}
				z++;
			}

			MPI_Isend(&aux, maxSize * n, MPI_INT, i, 0, MPI_COMM_WORLD, &requestNull);
		}

	    for(int pr = 1; pr < size; pr++) {
	      MPI_Irecv(&d, maxSize * n, MPI_INT, pr, 1, MPI_COMM_WORLD, &(requestList[pr - 1])); //&(requestList[pr - 1])
	    }
	    
	    for(int prW = 1; prW < size; prW++) {
	    	// the finished process's index
	       int index; 

	       MPI_Waitany(size - 1, requestList, &index, &status);

	       // mytodo
	       for (int i = 0; i < maxSize && index * (maxSize-1) + i <= m; i++) {
	       		for (int j = 0; j < n; j++) {
	       			c[(index) * maxSize + i][j] = d[i][j];
	       		}
	       }
	    }

	    printf("The result of a x b = \n");
	    for (int i = 0; i < m; i++) {
	    	for (int j = 0; j < n; j++) {
	    		printf("%d ", a[i][j]);
	    	}

	    	printf("   ");
	    	for (int j = 0; j < n; j++) {
	    		printf("%d ", b[i][j]);
	    	}

	    	printf("   ");
	    	for (int j = 0; j < n; j++) {
	    		printf("%d ", c[i][j]);
	    	}
	    	printf("\n");
	    }
	} else {
		MPI_Bcast(&b, n * p, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&maxSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

		int aux[maxSize][n];
		int d[maxSize][n];

	    MPI_Recv(&aux, maxSize * n, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

		for (int i = 0; i < maxSize && rank * (maxSize-1) + i <= m; i++) {
		    // multiply aux and b
			int z = 0, s;
			for (int j = 0; j < n ; j++) {
				s = 0;

				for (int z = 0; z < n ; z++) {
					s += aux[i][z] * b[z][j];					
 				}
 				
				d[i][j] = s;
			}
		}

		MPI_Send(&d, maxSize * n, MPI_INT, 0, 1, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	
	return 0;
}
