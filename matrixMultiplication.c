#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>

using namespace std;

int c[100][100];

// give m n p from command line
int main(int argc, char** argv) {
    int m, n, p;

    if (argc < 4 ) return -1;

    m = atoi(argv[1]);
    n = atoi(argv[2]);    
    p = atoi(argv[3]);

    int a[m][n], b[n][p]; // col.a = row.b for multiplication
	// int c[m][p]; // holding the result

    MPI_Init(&argc, &argv);

    int worldRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

	int tag = 201; // for mpi send
	int fragments = 0;

	int aux[m][n];

	// number of slaves
	int slaves = worldSize - 1;
    int maxSize;

	if (worldRank == 0) {
		// build a and b
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

		// bcast the matrix b - used by all nodes in multiplication
		MPI_Bcast(&b, n * p, MPI_INT, 0, MPI_COMM_WORLD);		

		// calculate fragment size of the matrix
		// x slaves = x fragments of max size N/x 
		// rows per slave = N / x = maxSize
		maxSize = m / slaves; // exclude the master
		if (m % slaves > 0) {
			maxSize++;
		}

		printf("Max size of each block sent to a slave %d \n", maxSize);
		
		// send size of aux matrix to all nodes
		MPI_Bcast(&maxSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

		MPI_Request request;
    	MPI_Status status;
		
		// exclude the master		
		for (int i = 1; i < worldSize && worldRank * (maxSize-1) + i <= m; i++) {
			int z = 0;
			for (int j = maxSize * (i - 1); j < maxSize * i && j < m; j++) {
			    for (int k = 0; k < n; k++) {
					aux[z][k] = a[j][k];
				}
				z++;
			}

			MPI_Isend(&aux[0][0], maxSize * n, MPI_INT, i, tag, MPI_COMM_WORLD, &request);
			MPI_Wait (&request, &status);
			//MPI_Send(&z, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
			//sleep(1);
		
		}

		sleep(5);
		printf("C:______________ \n");
		for (int i = 0; i < m; i++) {
		
			for (int j = 0; j < n; j++) {
				printf("%d ", c[i][j]);
			}
			printf("\n");
		}
		
			
	} else {
		int m_size;
		MPI_Request request;
    	MPI_Status status;

		MPI_Bcast(&b[0][0], n * p, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&maxSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
				
		MPI_Irecv(&aux[0][0], maxSize * n, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);
		//MPI_Recv(&m_size, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Wait (&request, &status);
		int y = 0;
		printf("Node %d of %d \n", worldRank, worldSize);
		for (int i = 0; i < maxSize && worldRank * (maxSize-1) + i <= m; i++) {

		    // multiply aux and b
			int z = 0; int s = 0;
			int j = 0;			
			for (j = 0; j < n ; j++) {
				s = 0;
				for (int z = 0; z < n ; z++) {
					s += aux[i][z] * b[z][j];					
 			}
				printf("c[%d][%d]=%d ",(worldRank-1)*maxSize + i,j,s);
				c[(worldRank-1)*maxSize + i][j] = s;

				
			}
			sleep(1);
			printf("\n");printf("\n");printf("\n");
			for (int i = 0; i < m; i++) {
		
				for (int j = 0; j < n; j++) {
					printf("%d ", c[i][j]);
				}
				printf("\n");
			}	
		}
	}
	
	// if (worldRank == 0) {
	// 	sleep(5);
	// 	printf("C:______________ \n");
	// 	for (int i = 0; i < m; i++) {
		
	// 		for (int j = 0; j < n; j++) {
	// 			printf("%d ", c[i][j]);
	// 		}
	// 		printf("\n");
	// 	}
	// }


	MPI_Finalize();	
    return 0;
}
