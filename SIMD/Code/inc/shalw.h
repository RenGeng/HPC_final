#include <string>
#include <math.h>
extern double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy, *hFil_global;
extern int size_x, size_y, nb_steps; // variable locale
extern int global_size_x, global_size_y; // variable globale
extern int height_bloc,width_bloc;
extern double dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
extern bool file_export;
extern std::string export_path;
extern int my_rank, np;



#define HFIL(t, i, j) hFil[ (j) +			\
			    (i) * width_bloc +		\
			    ((t)%2) * height_bloc * width_bloc ]

#define HFIL_global(t, i, j) hFil_global[ (j) +			\
			    (i) * global_size_y +		\
			    ((t)%2) * global_size_x * global_size_y ]

#define UFIL(t, i, j) uFil[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define VFIL(t, i, j) vFil[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define HPHY(t, i, j) hPhy[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define UPHY(t, i, j) uPhy[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define VPHY(t, i, j) vPhy[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
