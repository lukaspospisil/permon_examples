static char help[] = "";

#include <permonqp.h> /* manipulation with quadratic programming problems (QP) */
#include <permonqps.h> /* manipulation with solvers (QPS) */
#include "mpi.h" /* get rank of the procesor */

#include <iostream> /* std::cout */
#include <string> /* manipulation with strings */
#include <fstream> /* manipulation with files */
#include <sstream> /* manipulation with string streams */

/* default values */
#define PERMON_EXAMPLE_POLYTOPE_DEFAULT_M 3

void assemble_problem(int m, Mat *A, Vec *b, Vec *lb, Mat *B, Vec *c){
	PetscErrorCode ierr; /* error handler */

	PetscScalar *C_data;
	C_data = (PetscScalar *)malloc(4*m*sizeof(PetscScalar));

	PetscScalar t;
	for(int i=0;i<m;i++){
		t = i*2*M_PI/(double)m;
		
		/* P points */
		C_data[2*i+0] = -2.0 + cos(t);
		C_data[2*i+1] = 0.0 + sin(t);

		/* -Q points */
		C_data[2*(m+i)+0] = - (2.0 + cos(M_PI - t));
		C_data[2*(m+i)+1] = - (0.0 + sin(M_PI - t));
	}

	Mat C;
	
	ierr = MatCreateDense(PETSC_COMM_WORLD,2,2*m,2,2*m,C_data, &C);CHKERRV(ierr);
	ierr = PetscObjectSetName((PetscObject)C,"points");CHKERRV(ierr);
	ierr = MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
	ierr = MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
	
	/* A = 2*C'*C */
	ierr = MatTransposeMatMult(C, C, MAT_INITIAL_MATRIX , PETSC_DEFAULT, A);CHKERRV(ierr);
	ierr = MatTransposeMatMult(C, C, MAT_INITIAL_MATRIX , PETSC_DEFAULT, A);CHKERRV(ierr);
	ierr = MatScale(*A,2.0);CHKERRV(ierr);
	ierr = MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
	ierr = MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
	ierr = PetscObjectSetName((PetscObject)*A,"Hessian matrix");CHKERRV(ierr);
	
//	ierr = MatView(*A, PETSC_VIEWER_STDOUT_WORLD);CHKERRV(ierr);

	/* b = 0 */
	ierr = VecCreate(PETSC_COMM_WORLD, b);CHKERRV(ierr);
	ierr = VecSetSizes(*b,2*m,2*m);CHKERRV(ierr);
	ierr = VecSetFromOptions(*b);CHKERRV(ierr);
	ierr = VecSet(*b,0.0);CHKERRV(ierr);
	ierr = VecAssemblyBegin(*b);CHKERRV(ierr);
	ierr = VecAssemblyEnd(*b);CHKERRV(ierr);
	ierr = PetscObjectSetName((PetscObject)*b,"linear term");CHKERRV(ierr);

//	ierr = VecView(*b, PETSC_VIEWER_STDOUT_WORLD);CHKERRV(ierr);

	/* lb = 0 */
	ierr = VecDuplicate(*b,lb);CHKERRV(ierr);
	ierr = PetscObjectSetName((PetscObject)*lb,"lower bound");CHKERRV(ierr);
	
//	ierr = VecView(*lb, PETSC_VIEWER_STDOUT_WORLD);CHKERRV(ierr);

	/* B = [1 1 1, 0 0 0; 0 0 0 1 1 1]; */
	ierr = MatCreateDense(PETSC_COMM_WORLD,2,2*m,2,2*m,NULL, B);CHKERRV(ierr);
	ierr = PetscObjectSetName((PetscObject)*B,"matrix of EQ");CHKERRV(ierr);

	PetscScalar *B_data;
	ierr = MatDenseGetArray(*B,&B_data);CHKERRV(ierr);
	for(int i=0;i<m;i++){
		B_data[2*i] = 1.0;
		B_data[2*i+1] = 0.0;

		B_data[2*(m+i)] = 0.0;
		B_data[2*(m+i)+1] = 1.0;
	}
	ierr = MatDenseRestoreArray(*B,&B_data);CHKERRV(ierr);
	ierr = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
	ierr = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
	
//	ierr = MatView(*B, PETSC_VIEWER_STDOUT_WORLD);CHKERRV(ierr);

	/* c = 1 */
	ierr = VecCreate(PETSC_COMM_WORLD, c);CHKERRV(ierr);
	ierr = VecSetSizes(*c,2,2);CHKERRV(ierr);
	ierr = VecSetFromOptions(*c);CHKERRV(ierr);
	ierr = VecSet(*c,1.0);CHKERRV(ierr);
	ierr = VecAssemblyBegin(*c);CHKERRV(ierr);
	ierr = VecAssemblyEnd(*c);CHKERRV(ierr);
	ierr = PetscObjectSetName((PetscObject)*c,"rhs of EQ");CHKERRV(ierr);

//	ierr = VecView(*c, PETSC_VIEWER_STDOUT_WORLD);CHKERRV(ierr);

	free(C_data);
}

int main( int argc,char **args )
{
	/* initialize Petsc */
	PermonInitialize(&argc,&args,(char *)0,help);

	/* error handler */
	PetscErrorCode ierr;
  
	/* say hello */
	ierr = PetscPrintf(PETSC_COMM_WORLD,"This is POLYTOPE example\n");CHKERRQ(ierr);

	/* prepare variables */
	PetscInt m; /* number of discretisation nodes (size of the problem) */
	Vec y; /* solution vector */
	Mat A; /* Hessian matrix */
	Vec b; /* linear term */
	Vec lb; /* lower bound (rigid obstacle - non-penetration condition) */ 
	Mat B; /* matrix of equality constraints */
	Vec c; /* rhs of equality constraints */

	/* load console parameters */
	PetscBool set; /* the variable was provided in console parameters or not */
	ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,&set);CHKERRQ(ierr);
	if(!set) m = PERMON_EXAMPLE_POLYTOPE_DEFAULT_M;
	
	/* print parameters loaded from console */
	ierr = PetscPrintf(PETSC_COMM_WORLD," m                    = %d\n", m);CHKERRQ(ierr);

	/* assemble problem */
	assemble_problem(m, &A, &b, &lb, &B, &c);

	/* test */
//	ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//	ierr = VecView(b, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//	ierr = VecView(lb, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//	ierr = MatView(B, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//	ierr = VecView(c, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	/* prepare QP problem */
	QP qp;
	ierr = QPCreate(PETSC_COMM_WORLD, &qp);CHKERRQ(ierr);
	ierr = QPSetOperator(qp, A);CHKERRQ(ierr); /* set Hessian matrix */
	ierr = QPSetRhs(qp, b);CHKERRQ(ierr); /* set linear term vector */
	ierr = QPAddEq(qp,B, c);CHKERRQ(ierr); /* add equality constraints */
	ierr = QPSetBox(qp, NULL,lb, NULL);CHKERRQ(ierr); /* add lowerbound from contact */
	
	/* print some infos about QP */
//	ierr = QPView(qp, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	/* prepare QP solver */
	QPS qps;
	ierr = QPSCreate(PETSC_COMM_WORLD, &qps);CHKERRQ(ierr);
	ierr = QPSSetQP(qps, qp);CHKERRQ(ierr); /* Insert the QP problem into the solver. */
//	ierr = QPSSetTolerances(qps, setting.rtol, setting.atol, setting.dtol, setting.maxit);CHKERRQ(ierr); /* Set QPS options from settings */
//	ierr = QPSMonitorSet(qps,QPSMonitorDefault,NULL,0);CHKERRQ(ierr); /* Set the QPS monitor */
	ierr = QPTFromOptions(qp);CHKERRQ(ierr); /* Perform QP transforms */
	ierr = QPSSetFromOptions(qps);CHKERRQ(ierr); /* Set QPS options from the options database (overriding the defaults). */
	ierr = QPSSetUp(qps); CHKERRQ(ierr); /* Set up QP and QPS. */
	
	/* print some infos about QPS */
//	ierr = QPSView(qps, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	/* solve QP using QPS */
	ierr = QPSSolve(qps);CHKERRQ(ierr);

	/* get solution */
	ierr = QPGetSolutionVector(qp, &y);CHKERRQ(ierr);
	
//	ierr = VecView(y, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	/* compute the distance (function value) */
	Vec Ay;
	ierr = VecDuplicate(y,&Ay);CHKERRQ(ierr);
	ierr = MatMult(A,y,Ay);
	PetscScalar my_distance;
	ierr = VecDot(y,Ay,&my_distance);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Final distance       = %f (should be %f)\n",0.25*my_distance, 2.0);CHKERRQ(ierr);
	ierr = VecDestroy(&Ay);CHKERRQ(ierr);

	/* eigen value test */
	PetscScalar lambda;
	ierr = MatGetMaxEigenvalue(A, NULL, &lambda, 1e-10, 1e5);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Max eigenvalue test  = %f (should be %f)\n",lambda/(2.0*(double)m), 9.0);CHKERRQ(ierr);

	/* destroy the mess */
	ierr = QPDestroy(&qp);CHKERRQ(ierr);
	ierr = QPSDestroy(&qps);CHKERRQ(ierr);
	
	ierr = VecDestroy(&y);CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = VecDestroy(&b);CHKERRQ(ierr);
	ierr = VecDestroy(&lb);CHKERRQ(ierr);
	ierr = MatDestroy(&B);CHKERRQ(ierr);
	ierr = VecDestroy(&c);CHKERRQ(ierr);


	/* finalize Petsc */
	PermonFinalize();

	return 0;
}


