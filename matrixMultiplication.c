#include <mpi.h>
#include <cstdio>
#include <cstdlib>

using namespace std;

// give m n p from command line
int main(int argc, char** argv) {
    int m, n, p;

    if (argc < 4 ) return -1;

    m = atoi(argv[1]);
    n = atoi(argv[2]);    
    p = atoi(argv[3]);

    int a[m][n], b[n][p]; // col.a = row.b for multiplication
	int c[m][p]; // holding the result

    MPI_Init(&argc, &argv);

    int worldRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

	int tag = 201; // for mpi send
	int fragments = 0;

	printf("Node %d of %d \n", worldRank, worldSize);

	int aux[m][n];

	// number of slaves
	int slaves = worldSize - 1;
    int maxSize;

	if (worldRank == 0) {
		// build a and b
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				a[i][j] = rand() % 10 + 1; 				
			}
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < p; j++) {
				b[i][j] = rand() % 10 + 1; 				
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

		// exclude the master		
		for (int i = 1; i < worldSize && worldRank * (maxSize-1) + i <= m; i++) {
			int z = 0;
			for (int j = maxSize * (i - 1); j < maxSize * i && j < m; j++) {
			    for (int k = 0; k < n; k++) {
					aux[z][k] = a[j][k];
				}
				z++;
			}

			MPI_Send(&aux[0][0], maxSize * n, MPI_INT, i, tag, MPI_COMM_WORLD);
	
		
		}
		
			
	} else {
		MPI_Bcast(&b[0][0], n * p, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&maxSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
				
		MPI_Recv(&aux[0][0], maxSize * n, MPI_INT, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		 int y = 0;

		for (int i = 0; i < maxSize && worldRank * (maxSize-1) + i <= m; i++) {
			// multiply aux and b
			int z = 0; int s = 0;
			int j = 0;			
			for (j = 0; j < n ; j++) {
				for (int z = 0; z < n ; z++) {
			  		s += aux[i][z] * b[z][maxSize*(worldRank-1)+j];
					//printf("%d ", aux[i][j]);
				}
				printf("%d ", s);
				c[maxSize*(worldRank-1)+i][maxSize*(worldRank-1)+j] = s;
			}
printf("\n");
	

		}
		printf("\n");
		}


 	MPI_Finalize();          
	
	if (worldRank == 0) {	

	for (int i = 0; i < m; i++) {
		printf("%d -> ", i);
		for (int j = 0; j < n; j++) {
			printf("%d ", b[i][j]);
		}		
		printf("\n");
	}	
}
	
    return 0;
}
