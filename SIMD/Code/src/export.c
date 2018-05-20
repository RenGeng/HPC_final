#include <stdio.h>
#include <shalw.h>
#include <immintrin.h> 
FILE *create_file(void) {
  FILE *f;
  char fname[256];

  sprintf(fname, "%s/shalw_%dx%d_T%d.sav", export_path.c_str(), global_size_x, global_size_y, nb_steps);

  f = fopen(fname, "w+b");

  return f;
}

void export_step(FILE *f, int t) {
  fwrite((void *)&HFIL_global(t, 0, 0), sizeof(double), global_size_x * global_size_y, f);
}

void finalize_export(FILE *f) {
  fclose(f);
}
