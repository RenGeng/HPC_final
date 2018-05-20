#include <stdlib.h>
#include <shalw.h>
#include <parse_args.hpp>
#include <memory.h>
#include <init.h>
#include <forward.h>
#include <export.h>
#include <mpi.h>
#include <immintrin.h>

/* A EFFACER*/
#include <iostream>
#include <unistd.h>

double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy, *hFil_global;
double elaps_time;
int size_x, size_y, nb_steps;
int height_bloc,width_bloc;
int global_size_x, global_size_y;
double dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
bool file_export;
std::string export_path;
int my_rank, np;

double start,end; // Pour mesurer le temps

int main(int argc, char **argv) {
   // Initiation des processus
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&np);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if(my_rank==0) start = MPI_Wtime();
 

    parse_args(argc, argv);
    // printf("Command line options parsed my_rank=%d\n",my_rank);

    alloc();
    // printf("Memory allocated my_rank=%d\n",my_rank);

    

    gauss_init();
    // printf("State initialised my_rank=%d\n",my_rank);

  forward();
  // printf("State computed my_rank = %d\n",my_rank);

  dealloc();
  // printf("Memory freed\n");

  if (my_rank==0)
  {
    end = MPI_Wtime();
    printf("Temps d'exécution:%lf\n",end-start);
  }
  MPI_Finalize();
  return EXIT_SUCCESS;
}
