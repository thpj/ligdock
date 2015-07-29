#include <stdio.h>
#include <unistd.h>
#include "mpi.h"

int main(int argc, char **argv)

{
int *buf, i, rank, nints, len;
  char hostname[256];
  char command_string[256];

   MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Barrier(MPI_COMM_WORLD);
  //set i to number of receptor conformations. Rank is used as the ligand number. For example, HIV0_452 corresponds to HIV conformation 0 with ligand 452.
  //Be cetain to set up the directory file locations appropriate to where your files are for your docking experiments and your autogrid4 executable.
  for (i = 0; i < 10; i++) {
      sprintf(command_string, "cd Z_%d; ../../autogrid4 -p yourReceptorBasename%d_%d.gpf -l yourReceptorBasename%d_%d.glg", rank, i, rank, i, rank);
      system(command_string);
      MPI_Barrier(MPI_COMM_WORLD);
      }
      MPI_Finalize();
      printf("All Done\n");
      return 0;
} 

