all:
	make -C Découpage\ par\ bande/Non_bloquant/Code/
	make -C Découpage\ par\ bande/Bloquant/Code
	make -C Découpage\ par\ bloc/Non\ bloquant/Code/
	make -C Découpage\ par\ bloc/Bloquant/Code
	make -C OpenMPI/Code
	make -C SIMD/Code

clean:
	make clean -C Découpage\ par\ bande/Non_bloquant/Code/
	make clean -C Découpage\ par\ bande/Bloquant/Code
	make clean -C Découpage\ par\ bloc/Non\ bloquant/Code/
	make clean -C Découpage\ par\ bloc/Bloquant/Code
	make clean -C OpenMPI/Code
	make clean -C SIMD/Code