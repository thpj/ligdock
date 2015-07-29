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
  //set i to the number of receptor conformations. Rank is set to the ligand number. For example, hiv0_230 means HIV conformation 0, with ligand 230.
  for (i = 0; i < 10; i++) {
      sprintf(command_string, "cd Z_%d; ../../autodock4 -p yourReceptorBasename%d_%d.dpf -l yourReceptorBasename%d_%d.dlg", rank, i, rank, i, rank);
      system(command_string);
      MPI_Barrier(MPI_COMM_WORLD);
      }
      MPI_Finalize();
      printf("All Done\n");
      return 0;
} 

