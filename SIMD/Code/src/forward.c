#include <stdio.h>
#include <math.h>
#include <shalw.h>
#include <export.h>
#include <mpi.h>
#include <stdlib.h>
#include <omp.h>
#include <immintrin.h> 


#define TAG_HPHY_LIGNE 1
#define TAG_UPHY_LIGNE 2
#define TAG_VPHY_LIGNE 3
#define TAG_HPHY_COLONNE 4
#define TAG_UPHY_COLONNE 5
#define TAG_VPHY_COLONNE 6


void hFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //HPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    {
      __m256d hphy;
      hphy = _mm256_loadu_pd(&HPHY(t, i, j*4));
      _mm256_storeu_pd(&HFIL(t,i,j*4),hphy);
      // return HPHY(t, i, j);
    }
  else
    {
      __m256d hphy_t1,hphy_t2,hfil_t,resultat;
      hphy_t1 = _mm256_loadu_pd(&HPHY(t - 1, i, j*4));
      hphy_t2 = _mm256_loadu_pd(&HPHY(t, i, j*4));
      hfil_t = _mm256_loadu_pd(&HFIL(t - 1, i-1*(my_rank>=sqrt(np)), (j-1*(my_rank%(int)sqrt(np)!=0)) * 4));

      resultat = _mm256_mul_pd(_mm256_set1_pd(-2.0),hphy_t1);
      resultat = _mm256_add_pd(hfil_t,resultat);
      resultat = _mm256_add_pd(resultat,hphy_t2);
      resultat = _mm256_mul_pd(_mm256_set1_pd(alpha),resultat);
      resultat = _mm256_add_pd(resultat,hphy_t1);

      _mm256_storeu_pd(&HFIL(t, i-1*(my_rank>=sqrt(np)), (j-1*(my_rank%(int)sqrt(np)!=0)) * 4),resultat);

      // return HPHY(t - 1, i, j) +
      // alpha * (HFIL(t - 1, i-1*(my_rank>=sqrt(np)), j-1*(my_rank%(int)sqrt(np)!=0)) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
    }
  
}

void uFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //UPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    {
      __m256d uphy;
      uphy = _mm256_loadu_pd(&UPHY(t, i, j*4));
      _mm256_storeu_pd(&UFIL(t, i, j*4),uphy);
      // return UPHY(t, i, j);
    }
  else
    {
      __m256d uphy_t1, uphy_t2, ufil_t,resultat;
      uphy_t1 = _mm256_loadu_pd(&UPHY(t - 1, i, j*4));
      uphy_t2 = _mm256_loadu_pd(&UPHY(t, i, j*4));
      ufil_t = _mm256_loadu_pd(&UFIL(t - 1, i, j*4));

      resultat = _mm256_mul_pd(_mm256_set1_pd(-2.0),uphy_t1);
      resultat = _mm256_add_pd(resultat,ufil_t);
      resultat = _mm256_add_pd(resultat,uphy_t2);
      resultat = _mm256_mul_pd(_mm256_set1_pd(alpha),resultat);
      resultat = _mm256_add_pd(resultat,uphy_t1);

      _mm256_storeu_pd(&UFIL(t, i, j*4),resultat);
    }

  // return UPHY(t - 1, i, j) +
  //   alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
}

void vFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //VPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    {
      __m256d vphy;
      vphy = _mm256_loadu_pd(&VPHY(t, i, j*4));
      _mm256_storeu_pd(&VFIL(t,i,j*4),vphy);
      // return VPHY(t, i, j);
    }
  else
    {
      __m256d vphy_t1, vphy_t2, vfil_t, resultat;
      vphy_t1 = _mm256_loadu_pd(&VPHY(t - 1, i, j*4));
      vphy_t2 = _mm256_loadu_pd(&VPHY(t, i, j*4));
      vfil_t = _mm256_loadu_pd(&VFIL(t - 1, i, j*4));

      resultat = _mm256_mul_pd(_mm256_set1_pd(-2.0),vphy_t1);
      resultat = _mm256_add_pd(resultat,vfil_t);
      resultat = _mm256_add_pd(resultat,vphy_t2);
      resultat = _mm256_mul_pd(_mm256_set1_pd(alpha),resultat);
      resultat = _mm256_add_pd(resultat,vphy_t1);

      _mm256_storeu_pd(&VFIL(t,i,j*4),resultat);
    }


  // return VPHY(t - 1, i, j) +
  //   alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
}

void hPhy_forward(int t, int i, int j) {
  // double c, d;

  __m256d c;  
  if (i > 0)
    { 
      c = _mm256_loadu_pd(&UPHY(t - 1, i - 1, j));
      // c = UPHY(t - 1, i - 1, j);
    }
  else
    {
      c = _mm256_set1_pd(0.0);
      // c = 0.;
    }

  __m256d d;
  if (j < size_y - 3)
    {
      d = _mm256_loadu_pd(&VPHY(t - 1, i, (j + 1) * 4));
      // d = VPHY(t - 1, i, j + 1);
    }
  else
    {
      d = _mm256_set1_pd(0.0);
      // d = 0.;
    }

  __m256d hfil_t, uphy_t, vphy_t,res1,res2;
  hfil_t = _mm256_loadu_pd(&HFIL(t - 1, i-1*(my_rank>=sqrt(np)), (j-1*(my_rank%(int)sqrt(np)!=0)) * 4));
  uphy_t = _mm256_loadu_pd(&UPHY(t - 1, i, j*4));
  vphy_t = _mm256_loadu_pd(&VPHY(t - 1, i, j*4));

  res1 = _mm256_sub_pd(uphy_t,c);
  res1 = _mm256_mul_pd(res1,_mm256_set1_pd(1.0/dx));

  res2 = _mm256_sub_pd(d,vphy_t);
  res2 = _mm256_mul_pd(res2,_mm256_set1_pd(1.0/dy));
  res2 = _mm256_add_pd(res1,res2);
  res2 = _mm256_mul_pd(res2,_mm256_set1_pd(-dt*hmoy));
  res2 = _mm256_add_pd(res2,hfil_t);

  _mm256_storeu_pd(&HPHY(t,i,j),res2);

  // return HFIL(t - 1, i-1*(my_rank>=sqrt(np)), j-1*(my_rank%(int)sqrt(np)!=0)) - dt * hmoy * ((UPHY(t - 1, i, j) - c) / dx + (d - VPHY(t - 1, i, j)) / dy);
}

void uPhy_forward(int t, int i, int j) {
  // double b, e, f, g;
  
  if (i == size_x - 1)
    {
      __m256d uphy;
      uphy = _mm256_set1_pd(0.0);
      _mm256_storeu_pd(&UPHY(t,i,j*4),uphy);
      // return 0.;
    }
  else
    {

      __m256d b;
      if (i < size_x - 1)
	{
      
	  b = _mm256_loadu_pd(&HPHY(t - 1, i + 1, j*4));
	  // b = HPHY(t - 1, i + 1, j);
	}
      else
	{
	  b = _mm256_set1_pd(0.0);
	  // b = 0.;
	}

      __m256d e;
      if (j < size_y - 3)
	{
	  e = _mm256_loadu_pd(&VPHY(t - 1, i, (j + 1) * 4));
	  // e = VPHY(t - 1, i, j + 1);
	}
      else
	{
	  e = _mm256_set1_pd(0.0);
	}

      __m256d f;
      if (i < size_x - 1)
	{
	  f = _mm256_loadu_pd(&VPHY(t - 1, i + 1, j*4));
	  // f = VPHY(t - 1, i + 1, j);
	}
      else
	{
	  f = _mm256_set1_pd(0.0);
	  // f = 0.;
	}

      __m256d g;
      if (i < size_x - 1 && j < size_y - 3)
	{
	  g = _mm256_loadu_pd(&VPHY(t - 1, i + 1, (j + 1) * 4));
	  // g = VPHY(t - 1, i + 1, j + 1);
	}
      else
	{
	  g = _mm256_set1_pd(0.0);
	  // g = 0.;
	}

      __m256d ufil_t,hphy_t,vphy_t,res1,res2,res3;

      ufil_t = _mm256_loadu_pd(&UFIL(t - 1, i, j*4));
      hphy_t = _mm256_loadu_pd(&HPHY(t - 1, i, j*4));
      vphy_t = _mm256_loadu_pd(&VPHY(t - 1, i, j*4));

      res1 = _mm256_sub_pd(b,hphy_t);
      res1 = _mm256_mul_pd(_mm256_set1_pd(-grav / dx),res1);

      res2 = _mm256_add_pd(vphy_t,e);
      res2 = _mm256_add_pd(res2,f);
      res2 = _mm256_add_pd(res2,g);
      res2 = _mm256_mul_pd(_mm256_set1_pd(pcor/4.0),res2);
      res2 = _mm256_add_pd(res1,res2);

      res3 = _mm256_mul_pd(_mm256_set1_pd(dissip),ufil_t);
      res3 = _mm256_sub_pd(res2,res3);
      res3 = _mm256_mul_pd(_mm256_set1_pd(dt),res3);
      res3 = _mm256_add_pd(res3,ufil_t);

      _mm256_storeu_pd(&UPHY(t,i,j),res3);
      // return UFIL(t - 1, i, j) + dt * ((-grav / dx) * (b - HPHY(t - 1, i, j)) + (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) - (dissip * UFIL(t - 1, i, j)));
    }
}

void vPhy_forward(int t, int i, int j) {
  // double c, d, e, f;

  if (j == 0)
    {
      VPHY(t,i,j) = 0;
      // return 0.;
    }
  
  else
    {
      __m256d c;
      if (j > 0 && j < size_y - 3)
	{
	  c = _mm256_loadu_pd(&HPHY(t - 1, i, (j - 1) * 4));
	  // c = HPHY(t - 1, i, j - 1);
	}
      else
	{
	  c = _mm256_set1_pd(0.0);
	  // c = 0.;
	}

      __m256d d;
      if (i > 0 && j > 0 && j < size_y - 3)
	{
	  d = _mm256_loadu_pd(&UPHY(t - 1, i -1, (j -1) * 4));
	  // d = UPHY(t - 1, i -1, j -1);
	}
      else
	{
	  d = _mm256_set1_pd(0.0);
	  // d = 0.;
	}

      __m256d e;
      if (i > 0)
	{
	  d = _mm256_loadu_pd(&UPHY(t - 1, i - 1, j * 4));
	  // e = UPHY(t - 1, i - 1, j);
	}
      else
	{
	  e = _mm256_set1_pd(0.0);
	  // e = 0.;
	}

      __m256d f;
      if (j > 0 && j < size_y - 3)
	{
	  f = _mm256_loadu_pd(&UPHY(t - 1, i, (j - 1) * 4));
	  // f = UPHY(t - 1, i, j - 1);
	}
      else
	{
	  f = _mm256_set1_pd(0.0);
	  // f = 0.;
	}

      __m256d vfil_t,hphy_t,uphy_t,res1,res2,res3;

      vfil_t = _mm256_loadu_pd(&VFIL(t - 1, i, j*4));
      hphy_t = _mm256_loadu_pd(&HPHY(t - 1, i, j*4));
      uphy_t = _mm256_loadu_pd(&UPHY(t - 1, i, j*4));

      res1 = _mm256_sub_pd(hphy_t,c);
      res1 = _mm256_mul_pd(_mm256_set1_pd(-grav / dy),res1);

      res2 = _mm256_add_pd(uphy_t,e);
      res2 = _mm256_add_pd(res2,f);
      res2 = _mm256_add_pd(res2,d);
      res2 = _mm256_mul_pd(_mm256_set1_pd(pcor/4.0),res2);
      res2 = _mm256_sub_pd(res1,res2);

      res3 = _mm256_mul_pd(_mm256_set1_pd(dissip),vfil_t);
      res3 = _mm256_sub_pd(res2,res3);
      res3 = _mm256_mul_pd(_mm256_set1_pd(dt),res3);
      res3 = _mm256_add_pd(res3,vfil_t);

      _mm256_storeu_pd(&UPHY(t,i,j),res3);

      // return VFIL(t - 1, i, j) + dt * ((-grav / dy) * (HPHY(t - 1, i, j) - c) - (pcor / 4.) * (d + e + f + UPHY(t - 1, i, j)) - (dissip * VFIL(t - 1, i, j)));
    }
  
}

void forward(void) {
  FILE *file = NULL;
  double svdt = 0.;
  int t = 0;
  MPI_Status status;

  MPI_Datatype colonne,bloc_envoie;
  MPI_Type_vector(height_bloc,1,size_y,MPI_DOUBLE,&colonne);
  MPI_Type_commit(&colonne);


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

	  // disp[i] = i*128;
	  if(i%(int)sqrt(np)==0) disp[i] = i/(int)sqrt(np) * global_size_y * height_bloc;
	  else disp[i] = disp[i-1] + width_bloc;

	}

    }

  	
  if(my_rank%(int)sqrt(np)==0 || my_rank%(int)sqrt(np)==sqrt(np)-1)
    {
      MPI_Type_vector(height_bloc,width_bloc,1,MPI_DOUBLE,&bloc_test);
      MPI_Type_create_resized(bloc_test, 0,  sizeof(double), &bloc_envoie);
      MPI_Type_commit(&bloc_envoie);
    }
  else
    {
      MPI_Type_vector(height_bloc,width_bloc,2,MPI_DOUBLE,&bloc_test);
      MPI_Type_create_resized(bloc_test, 0, sizeof(double), &bloc_envoie);
      MPI_Type_commit(&bloc_envoie);
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
    
#pragma omp parallel for
      for (int i = (my_rank<sqrt(np)) ? 0:1; i <((my_rank<sqrt(np))*height_bloc + (my_rank>=sqrt(np))*(height_bloc+1))/4;  i++)    
  	{
	  // #pragma omp parallel for
	  for (int j = (my_rank%(int)sqrt(np)==0) ? 0:1; j < ((my_rank%(int)sqrt(np)==0)*width_bloc + (my_rank%(int)sqrt(np)!=0)*(width_bloc+1))/4; j++)
  	    {   
	      {
		hPhy_forward(t, i, j);
		uPhy_forward(t, i, j);
		vPhy_forward(t, i, j);
		// HFIL(t, i, j) = hFil_forward(t, i, j);
		hFil_forward(t, i, j);
		uFil_forward(t, i, j);
		vFil_forward(t, i, j);
	      }

  	    }
  	}

      if(my_rank>=sqrt(np))
  	{
  		
	  MPI_Send(&HPHY(t,1,(my_rank%(int)sqrt(np)==0) ? 0:1),width_bloc,MPI_DOUBLE,my_rank-(int)sqrt(np),TAG_HPHY_LIGNE,MPI_COMM_WORLD);
	  MPI_Send(&VPHY(t,1,(my_rank%(int)sqrt(np)==0) ? 0:1),width_bloc,MPI_DOUBLE,my_rank-(int)sqrt(np),TAG_VPHY_LIGNE,MPI_COMM_WORLD);
	  MPI_Recv(&UPHY(t,0,(my_rank%(int)sqrt(np)==0) ? 0:1),width_bloc,MPI_DOUBLE,my_rank-(int)sqrt(np),TAG_UPHY_LIGNE,MPI_COMM_WORLD,&status);

	  // printf("Process %d a envoyé à %d HPHY:%lf\n",my_rank,my_rank-(int)sqrt(np),HPHY(t,1,(my_rank%(int)sqrt(np)==0) ? 0:1));
  	}

      if(my_rank<np-sqrt(np))
  	{
  		
	  MPI_Send(&UPHY(t,size_x-2,(my_rank%(int)sqrt(np)==0) ? 0:1),width_bloc,MPI_DOUBLE,my_rank+(int)sqrt(np),TAG_UPHY_LIGNE,MPI_COMM_WORLD);
	  MPI_Recv(&HPHY(t,size_x-1,(my_rank%(int)sqrt(np)==0) ? 0:1),width_bloc,MPI_DOUBLE,my_rank+(int)sqrt(np),TAG_HPHY_LIGNE,MPI_COMM_WORLD,&status);
	  MPI_Recv(&VPHY(t,size_x-1,(my_rank%(int)sqrt(np)==0) ? 0:1),width_bloc,MPI_DOUBLE,my_rank+(int)sqrt(np),TAG_VPHY_LIGNE,MPI_COMM_WORLD,&status);

	  // printf("Process %d a reçu HPHY:%lf\n",my_rank,HPHY(t,size_x-1,(my_rank%(int)sqrt(np)==0) ? 0:1));
  	}

      if(my_rank%(int)sqrt(np)!=sqrt(np)-1)
  	{
	  MPI_Send(&HPHY(t,(my_rank<sqrt(np) ? 0:1),size_y-2),1,colonne,my_rank+1,TAG_HPHY_COLONNE,MPI_COMM_WORLD);
	  MPI_Send(&UPHY(t,(my_rank<sqrt(np) ? 0:1),size_y-2),1,colonne,my_rank+1,TAG_UPHY_COLONNE,MPI_COMM_WORLD);
	  MPI_Recv(&VPHY(t,(my_rank<sqrt(np) ? 0:1),size_y-1),1,colonne,my_rank+1,TAG_VPHY_COLONNE,MPI_COMM_WORLD,&status);

	  // printf("Process %d a envoyé à process %d HPHY:%lf\n",my_rank,my_rank+1,HPHY(t,30,size_y-2));
  	}

      if(my_rank%(int)sqrt(np)!=0)
  	{
	  MPI_Send(&VPHY(t,(my_rank<sqrt(np) ? 0:1),1),1,colonne,my_rank-1,TAG_VPHY_COLONNE,MPI_COMM_WORLD);
	  MPI_Recv(&HPHY(t,(my_rank<sqrt(np) ? 0:1),0),1,colonne,my_rank-1,TAG_HPHY_COLONNE,MPI_COMM_WORLD,&status);
	  MPI_Recv(&UPHY(t,(my_rank<sqrt(np) ? 0:1),0),1,colonne,my_rank-1,TAG_UPHY_COLONNE,MPI_COMM_WORLD,&status);

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
