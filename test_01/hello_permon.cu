/* include petsc */
#include "petsc.h"
#include "mpi.h"
#include "fllopqp.h"

int main( int argc, char *argv[] )
{
	/* initialize Petsc */
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

	TRY( PetscPrintf(PETSC_COMM_WORLD,"Hello from PERMON!\n") );
	
	/* give info about MPI */
	int size, rank; /* size and rank of communicator */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	TRY( PetscPrintf(MPI_COMM_WORLD,"- number of processors: %d\n",size) );
	TRY( PetscSynchronizedPrintf(MPI_COMM_WORLD," - hello from processor: %d\n",rank) );
	TRY( PetscSynchronizedFlush(MPI_COMM_WORLD,NULL) );
	TRY( PetscPrintf(PETSC_COMM_WORLD,"-------------------------------\n") );

	/* finalize Petsc */
	PetscFinalize();

	return 0;
}


