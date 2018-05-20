#include <math.h>
#include <shalw.h>
#include <stdlib.h>
#include <stdio.h>

void gauss_init(void) {
  double gmx, gmy, gsx, gsy;

  gmx = global_size_x * dx / 2 ;
  gmy = global_size_y * dy / 2 ;
  gsx = 25000 ;
  gsy = 25000 ;

  // printf("process %d (my_rank<sqrt(np))*height_bloc + (my_rank>=sqrt(np))*(height_bloc+1) = %d\n",my_rank,(my_rank<sqrt(np))*height_bloc + (my_rank>=sqrt(np))*(height_bloc+1));


  for (int i = 0; i <height_bloc;  i++) {
    for (int j = 0; j < width_bloc; j++) {
      HFIL(0, i, j) = height *
	(exp(- pow(((i+((int)my_rank/(int)sqrt(np)*height_bloc)) * dx - gmx) / gsx, 2) / 2.)) *
	(exp(- pow(((j+(my_rank%(int)sqrt(np)*width_bloc)) * dy - gmy) / gsy, 2) / 2.)) ;

    }
  }
}
