/* include petsc */
#include "petsc.h"
#include "mpi.h"
#include "fllopqp.h"

/* default values */
#define PERMON_EXAMPLES_DEFAULT_N 10
#define PERMON_EXAMPLES_DEFAULT_M 2

void assemble_problem(int n, int m, Mat K, Vec f, Vec lb){
  /* error handler */
  PetscErrorCode ierr;  

  PetscInt i,j,II,JJ,Istart,Iend; /* iterators */
  PetscReal value; 
  
  /* assemble the system matrix */
  ierr = MatCreate(PETSC_COMM_WORLD,&K);CHKERRV(ierr);
  ierr = MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,m*n,m*n);CHKERRV(ierr);
  ierr = MatSetFromOptions(K);CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(K,5,NULL,5,NULL);CHKERRV(ierr);
  ierr = MatSeqAIJSetPreallocation(K,5,NULL);CHKERRV(ierr);
  ierr = MatGetOwnershipRange(K,&Istart,&Iend);CHKERRV(ierr);
  for (II = Istart; II < Iend; II++) { 
    value = -1.0; i = II/n; j = II - i*n;  
    if (i>0)   {JJ = II - n; ierr = MatSetValues(K,1,&II,1,&JJ,&value,INSERT_VALUES);CHKERRV(ierr);}
    if (i<m-1) {JJ = II + n; ierr = MatSetValues(K,1,&II,1,&JJ,&value,INSERT_VALUES);CHKERRV(ierr);}
    if (j>0)   {JJ = II - 1; ierr = MatSetValues(K,1,&II,1,&JJ,&value,INSERT_VALUES);CHKERRV(ierr);}
    if (j<n-1) {JJ = II + 1; ierr = MatSetValues(K,1,&II,1,&JJ,&value,INSERT_VALUES);CHKERRV(ierr);}
    value = 4.0; ierr = MatSetValues(K,1,&II,1,&II,&value,INSERT_VALUES);CHKERRV(ierr);
  }
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
  ierr = MatSetOption(K,MAT_SYMMETRIC,PETSC_TRUE);CHKERRV(ierr);

}

int main( int argc, char *argv[] )
{
	/* initialize Petsc */
	FllopInitialize(&argc,&argv,PETSC_NULL);

	/* error handler */
	PetscErrorCode ierr;  
  
	/* say hello */
	ierr = PetscPrintf(PETSC_COMM_WORLD,"This is STRING example\n");CHKERRQ(ierr);

	/* prepare variables */
	PetscInt n; /* number of discretisation nodes (size of the problem) */
	PetscInt m;
	Mat K; /* stiffness matrix */
	Vec f; /* vector of external forces */
	Vec lb; /* lower bound (rigid obstacle - non-penetration condition) */ 
	Vec u; /* solution vector */

	/* load console parameters */
	PetscBool set; /* the variable was provided in console parameters or not */
	ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,&set);CHKERRQ(ierr);
	if(!set) n = PERMON_EXAMPLES_DEFAULT_N;
	
	ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,&set);CHKERRQ(ierr);
	if(!set) m = PERMON_EXAMPLES_DEFAULT_M;


	/* print parameters loaded from console */
	ierr = PetscPrintf(PETSC_COMM_WORLD," n             = %d\n", n);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," m             = %d\n", m);CHKERRQ(ierr);
	
	
	/* finalize Petsc */
	FllopFinalize();

	return 0;
}


