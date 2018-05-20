#include <stdio.h>
#include <math.h>
#include <shalw.h>
#include <export.h>
#include <mpi.h>
#include <stdlib.h>


#define TAG_HPHY_LIGNE 1
#define TAG_UPHY_LIGNE 2
#define TAG_VPHY_LIGNE 3
#define TAG_HPHY_COLONNE 4
#define TAG_UPHY_COLONNE 5
#define TAG_VPHY_COLONNE 6


double hFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //HPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return HPHY(t, i, j);
  return HPHY(t - 1, i, j) +
    alpha * (HFIL(t - 1, i-1*(my_rank>=sqrt(np)), j-1*(my_rank%(int)sqrt(np)!=0)) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
}

double uFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //UPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return UPHY(t, i, j);
  return UPHY(t - 1, i, j) +
    alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
}

double vFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //VPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return VPHY(t, i, j);
  return VPHY(t - 1, i, j) +
    alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
}

double hPhy_forward(int t, int i, int j) {
  double c, d;
  
  c = 0.;
  if (i > 0)
    c = UPHY(t - 1, i - 1, j);

  d = 0.;
  if (j < size_y - 1)
    d = VPHY(t - 1, i, j + 1);

  return HFIL(t - 1, i-1*(my_rank>=sqrt(np)), j-1*(my_rank%(int)sqrt(np)!=0)) -
    dt * hmoy * ((UPHY(t - 1, i, j) - c) / dx +
		 (d - VPHY(t - 1, i, j)) / dy);
}

double uPhy_forward(int t, int i, int j) {
  double b, e, f, g;
  
  if (i == size_x - 1)
    return 0.;

  b = 0.;
  if (i < size_x - 1)
    b = HPHY(t - 1, i + 1, j);

  e = 0.;
  if (j < size_y - 1)
    e = VPHY(t - 1, i, j + 1);

  f = 0.;
  if (i < size_x - 1)
    f = VPHY(t - 1, i + 1, j);

  g = 0.;
  if (i < size_x - 1 && j < size_y - 1)
    g = VPHY(t - 1, i + 1, j + 1);

  return UFIL(t - 1, i, j) +
    dt * ((-grav / dx) * (b - HPHY(t - 1, i, j)) +
	  (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) -
	  (dissip * UFIL(t - 1, i, j)));
}

double vPhy_forward(int t, int i, int j) {
  double c, d, e, f;

  if (j == 0)
    return 0.;

  c = 0.;
  if (j > 0)
    c = HPHY(t - 1, i, j - 1);

  d = 0.;
  if (i > 0 && j > 0)
    d = UPHY(t - 1, i -1, j -1);

  e = 0.;
  if (i > 0)
    e = UPHY(t - 1, i - 1, j);

  f = 0.;
  if (j > 0)
    f = UPHY(t - 1, i, j - 1);

  return VFIL(t - 1, i, j) +
    dt * ((-grav / dy) * (HPHY(t - 1, i, j) - c) -
	  (pcor / 4.) * (d + e + f + UPHY(t - 1, i, j)) -
	  (dissip * VFIL(t - 1, i, j)));
}

void forward(void) {
  FILE *file = NULL;
  double svdt = 0.;
  int t = 0;
  MPI_Status status;
  MPI_Request send,HPHY_ligne,VPHY_ligne,VPHY_colonne;

  /* Création de data type pour envoyer les colonnes*/
  MPI_Datatype colonne;
  MPI_Type_vector(height_bloc,1,size_y,MPI_DOUBLE,&colonne);
  MPI_Type_commit(&colonne);

  /* Data type pour le gather */
  MPI_Datatype bloc_gather,bloc_test;
  if(my_rank==0)
    {
      MPI_Type_vector(height_bloc,width_bloc,global_size_y,MPI_DOUBLE,&bloc_test);
      MPI_Type_create_resized(bloc_test, 0, sizeof(double), &bloc_gather);
      MPI_Type_commit(&bloc_gather);
    }

  int disp[np],count[np];
  if(my_rank==0)
    {
      disp[0] = 0;
      count[0] = 1;

      for(int i = 1;i<np;i++)
	{
	  count[i] = 1;
	  if(i%(int)sqrt(np)==0) disp[i] = i/(int)sqrt(np) * global_size_y * height_bloc;
	  else disp[i] = disp[i-1] + width_bloc;
	}

    }

  MPI_Gatherv(&HFIL(t,0,0),height_bloc*width_bloc,MPI_DOUBLE,&HFIL_global(t,0,0),count,disp,bloc_gather,0,MPI_COMM_WORLD);

  if (file_export && my_rank == 0) 
    {
      file = create_file();
      export_step(file, t);
    }
   
  for (t = 1; t < nb_steps; t++) 
    {  
      
      if (t == 1) 
  	{
  	  svdt = dt;
  	  dt = 0;
  	}
      if (t == 2)
  	{
  	  dt = svdt / 2.;
  	}
    
      for (int i = (my_rank<sqrt(np)) ? 0:1; i <(my_rank<sqrt(np))*height_bloc + (my_rank>=sqrt(np))*(height_bloc+1);  i++)
  	{
	  for (int j = (my_rank%(int)sqrt(np)==0) ? 0:1; j < (my_rank%(int)sqrt(np)==0)*width_bloc + (my_rank%(int)sqrt(np)!=0)*(width_bloc+1); j++)
  	    { 
          
	      /* Receive bloquant car dans les calculs on a besoin de i-1 de UPHY et j-1 de UPHY et HPHY */
	      if(i==0+(my_rank>=sqrt(np))*1 && j==0+(my_rank%(int)sqrt(np)!=0)*1 && my_rank>=sqrt(np))
		{
		  MPI_Recv(&UPHY(t,0,(my_rank%(int)sqrt(np)==0) ? 0:1),width_bloc,MPI_DOUBLE,my_rank-(int)sqrt(np),TAG_UPHY_LIGNE,MPI_COMM_WORLD,&status);
		}

	      /* Receive bloquant car dans les calculs on a besoin de j-1 de UPHY et HPHY */
	      if(i==0+(my_rank>=sqrt(np))*1 && j==0+(my_rank%(int)sqrt(np)!=0)*1 && my_rank%(int)sqrt(np)!=0)
		{
		  MPI_Recv(&HPHY(t,(my_rank<sqrt(np) ? 0:1),0),1,colonne,my_rank-1,TAG_HPHY_COLONNE,MPI_COMM_WORLD,&status);
		  MPI_Recv(&UPHY(t,(my_rank<sqrt(np) ? 0:1),0),1,colonne,my_rank-1,TAG_UPHY_COLONNE,MPI_COMM_WORLD,&status);
		}

	      /* On bloque quand on arrive à la dernière ligne car HPHY et VPHY ont besoin de i+1 */
	      if(i==((my_rank<sqrt(np))*height_bloc + (my_rank>=sqrt(np))*(height_bloc+1))-1 && j==0+(my_rank%(int)sqrt(np)!=0)*1 && t>1 && my_rank<np-sqrt(np))
		{
		  MPI_Wait(&HPHY_ligne,NULL);
		  MPI_Wait(&VPHY_ligne,NULL);
		}

	      /* On bloque quand on arrive à la dernière ligne car VPHY a besoin de j+1 */
	      if(j==((my_rank%(int)sqrt(np)==0)*width_bloc + (my_rank%(int)sqrt(np)!=0)*(width_bloc+1))-1 && i==0+(my_rank>=sqrt(np))*1 && t>1 && my_rank%(int)sqrt(np)!=sqrt(np)-1)
		{
		  MPI_Wait(&VPHY_colonne,NULL);
		}

  	      HPHY(t, i, j) = hPhy_forward(t, i, j);
  	      UPHY(t, i, j) = uPhy_forward(t, i, j);
  	      VPHY(t, i, j) = vPhy_forward(t, i, j);
  	      // HFIL(t, i, j) = hFil_forward(t, i, j);
  	      HFIL(t, i-1*(my_rank>=sqrt(np)), j-1*(my_rank%(int)sqrt(np)!=0))= hFil_forward(t, i, j);

  	      UFIL(t, i, j) = uFil_forward(t, i, j);
  	      VFIL(t, i, j) = vFil_forward(t, i, j);
  	    }
  	}


      /* Envoie des lignes */
      if(my_rank>=sqrt(np))
  	{
	  MPI_Isend(&HPHY(t,1,(my_rank%(int)sqrt(np)==0) ? 0:1),width_bloc,MPI_DOUBLE,my_rank-(int)sqrt(np),TAG_HPHY_LIGNE,MPI_COMM_WORLD,&send);
	  MPI_Isend(&VPHY(t,1,(my_rank%(int)sqrt(np)==0) ? 0:1),width_bloc,MPI_DOUBLE,my_rank-(int)sqrt(np),TAG_VPHY_LIGNE,MPI_COMM_WORLD,&send);
  		
	  // printf("Process %d a envoyé à %d HPHY:%lf\n",my_rank,my_rank-(int)sqrt(np),HPHY(t,1,(my_rank%(int)sqrt(np)==0) ? 0:1));
  	}

      if(my_rank<np-sqrt(np))
  	{
	  MPI_Isend(&UPHY(t,size_x-2,(my_rank%(int)sqrt(np)==0) ? 0:1),width_bloc,MPI_DOUBLE,my_rank+(int)sqrt(np),TAG_UPHY_LIGNE,MPI_COMM_WORLD,&send);
	  MPI_Irecv(&HPHY(t,size_x-1,(my_rank%(int)sqrt(np)==0) ? 0:1),width_bloc,MPI_DOUBLE,my_rank+(int)sqrt(np),TAG_HPHY_LIGNE,MPI_COMM_WORLD,&HPHY_ligne);
	  MPI_Irecv(&VPHY(t,size_x-1,(my_rank%(int)sqrt(np)==0) ? 0:1),width_bloc,MPI_DOUBLE,my_rank+(int)sqrt(np),TAG_VPHY_LIGNE,MPI_COMM_WORLD,&VPHY_ligne);

	  // printf("Process %d a reçu HPHY:%lf\n",my_rank,HPHY(t,size_x-1,(my_rank%(int)sqrt(np)==0) ? 0:1));
  	}


      /* Envoie des colonnes */
      if(my_rank%(int)sqrt(np)!=sqrt(np)-1)
  	{
	  MPI_Isend(&HPHY(t,(my_rank<sqrt(np) ? 0:1),size_y-2),1,colonne,my_rank+1,TAG_HPHY_COLONNE,MPI_COMM_WORLD,&send);
	  MPI_Isend(&UPHY(t,(my_rank<sqrt(np) ? 0:1),size_y-2),1,colonne,my_rank+1,TAG_UPHY_COLONNE,MPI_COMM_WORLD,&send);
	  MPI_Irecv(&VPHY(t,(my_rank<sqrt(np) ? 0:1),size_y-1),1,colonne,my_rank+1,TAG_VPHY_COLONNE,MPI_COMM_WORLD,&VPHY_colonne);

	  // printf("Process %d a envoyé à process %d HPHY:%lf\n",my_rank,my_rank+1,HPHY(t,30,size_y-2));
  	}

      if(my_rank%(int)sqrt(np)!=0)
  	{
	  MPI_Isend(&VPHY(t,(my_rank<sqrt(np) ? 0:1),1),1,colonne,my_rank-1,TAG_VPHY_COLONNE,MPI_COMM_WORLD,&send);

	  // printf("Process %d a reçu HPHY:%lf \n",my_rank,HPHY(t,30,0));
  	}

      MPI_Gatherv(&HFIL(t,0,0),height_bloc*width_bloc,MPI_DOUBLE,&HFIL_global(t,0,0),count,disp,bloc_gather,0,MPI_COMM_WORLD);
      if (file_export && my_rank==0) 
	{
	  export_step(file, t);
	}
        
      if (t == 2) 
	{
	  dt = svdt;
	}
    }

  if (file_export&&my_rank==0) 
    {
      finalize_export(file);
    }
  // printf("Fin forward process %d\n",my_rank);
}
