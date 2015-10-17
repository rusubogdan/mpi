// Task 1. Sum(1/n!)
// Author: Anton Danshin
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
 
int main(int argc, char** argv) {
  double result;
  int n = 0;
  if(argc>1)
    // atoi convert string to integer
    n = atoi(argv[1]);
 
  double startwtime = 0.0;
  double endwtime;
 
  // Initialize the MPI environment - cmd args not required
  MPI_Init(&argc, &argv);
  // Find out rank, size - rank of the process - 0 is master
  int worldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  int worldSize;
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
 
  // this is the master process
  if(worldRank == 0)
    // double - time in seconds 
    startwtime = MPI_Wtime();

  // void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm 
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
 
  double a = 0.0;
  double f = 1.0;
  double s = 0.0;
 
  int i;
  // n = 0 ...
  int t = n/worldSize;

  printf("Execution of process %d of %d \n", worldRank + 1, worldSize);
  printf("n = %s ", argv[0]);

  // -np N - N means number pf processes - the number of cores to be used is system dependant
  if(n % worldSize)
	t++;

  

  for(i = (n / worldSize) * worldRank; i < worldRank * (n / worldSize) + t; i++) {
    printf("...");
    if(i + 1 > n)
       break;
    f /= (double)(i + 1);
    s += f;
  }
 
  // Receive A and F from previous node
  double a_p, f_p = 1.0;
  if(worldRank > 0) {
    MPI_Recv(&f_p, 1, MPI_DOUBLE, worldRank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    f *= f_p;
  }
  s *=f_p;
 
  // Send my A and F to the next node
  if(worldRank < worldSize-1) {
    MPI_Send(&f, 1, MPI_DOUBLE, worldRank+1, 0, MPI_COMM_WORLD);
  }
 
  MPI_Reduce(&s, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
 
  if(worldRank==0) {
    endwtime = MPI_Wtime();
    printf("Result: %.10f\n", result+1);
    printf("Time elapsed: %f seconds \n", (endwtime-startwtime)*10);
  }
 
  MPI_Finalize();
}
