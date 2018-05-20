#include <stdlib.h>
#include <shalw.h>


#include <stdio.h>

void alloc(void) {

   height_bloc = global_size_x/sqrt(np);
   width_bloc = global_size_y/sqrt(np);
   size_x=(my_rank<sqrt(np) || my_rank>= np - sqrt(np)) ? (height_bloc+1) : (height_bloc+2);
   size_y=(my_rank%(int)sqrt(np)==0 || my_rank%(int)sqrt(np)==sqrt(np)-1) ? (width_bloc+1) : (width_bloc+2);

   // printf("Process %d et size_x=%d size_y=%d\n",my_rank,size_x,size_y);

   	if(my_rank==0) hFil_global = (double *) calloc(2*global_size_x*global_size_y, sizeof(double));

	  hFil = (double *) calloc(2*height_bloc*width_bloc, sizeof(double));
	  uFil = (double *) calloc(2*size_x*size_y, sizeof(double));
	  vFil = (double *) calloc(2*size_x*size_y, sizeof(double));
	  hPhy = (double *) calloc(2*size_x*size_y, sizeof(double));
	  uPhy = (double *) calloc(2*size_x*size_y, sizeof(double));
	  vPhy = (double *) calloc(2*size_x*size_y, sizeof(double));

}

void dealloc(void) {
	if(my_rank==0) free(hFil_global);
  free(hFil);
  free(uFil);
  free(vFil);
  free(hPhy);
  free(uPhy);
  free(vPhy);
}
