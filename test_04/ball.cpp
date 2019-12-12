static char help[] = "";

#include <permonqp.h> /* manipulation with quadratic programming problems (QP) */
#include <permonqps.h> /* manipulation with solvers (QPS) */
#include "mpi.h" /* get rank of the procesor */

#include <iostream> /* std::cout */
#include <string> /* manipulation with strings */
#include <fstream> /* manipulation with files */
#include <sstream> /* manipulation with string streams */

/* default values */
#define PERMON_EXAMPLE_BALL_DEFAULT_D 2
#define PERMON_EXAMPLE_BALL_DEFAULT_M 100

void assemble_problem(int d, int m, Mat *C, Mat *A, Vec *b, Vec *lb, Mat *B, Vec *c){
	PetscErrorCode ierr; /* error handler */

	/* b = dot */
	ierr = VecCreate(PETSC_COMM_WORLD, b);CHKERRV(ierr);
	ierr = VecSetSizes(*b,m,m);CHKERRV(ierr);
	ierr = VecSetFromOptions(*b);CHKERRV(ierr);

	/* C */
	ierr = MatCreateDense(PETSC_COMM_WORLD,d,m,d,m,NULL, C);CHKERRV(ierr);
	
	PetscScalar *b_array;
	ierr = VecGetArray(*b,&b_array);CHKERRV(ierr);

	PetscScalar *C_data;
	ierr = MatDenseGetArray(*C, &C_data);CHKERRV(ierr);

	PetscScalar num[d];
	PetscScalar dist;
	PetscScalar mydot;
	for(int j=0;j<m;j++){
		do {
			dist = 0;
			for(int i=0;i<d;i++){
				num[i] = 0.0 + (double)rand()/RAND_MAX*2.5;
				dist += (num[i] - 1.0)*(num[i] - 1.0);
			}
			dist = sqrt(dist);
		} while(dist > 1);

		mydot = 0.0;
		for(int i=0;i<d;i++){
			C_data[d*j+i] = num[i];
			mydot += num[i]*num[i];
		}
		b_array[j] = mydot;
	}

	ierr = VecRestoreArray(*b,&b_array);CHKERRV(ierr);
	ierr = VecAssemblyBegin(*b);CHKERRV(ierr);
	ierr = VecAssemblyEnd(*b);CHKERRV(ierr);
	ierr = PetscObjectSetName((PetscObject)*b,"linear term");CHKERRV(ierr);

	ierr = MatDenseRestoreArray(*C, &C_data);CHKERRV(ierr);
	ierr = PetscObjectSetName((PetscObject)*C,"points");CHKERRV(ierr);
	ierr = MatAssemblyBegin(*C,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
	ierr = MatAssemblyEnd(*C,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
	
	/* A = 2*C'*C */
	ierr = MatTransposeMatMult(*C, *C, MAT_INITIAL_MATRIX , PETSC_DEFAULT, A);CHKERRV(ierr);
	ierr = MatScale(*A,2.0);CHKERRV(ierr);
	ierr = MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
	ierr = MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
	ierr = PetscObjectSetName((PetscObject)*A,"Hessian matrix");CHKERRV(ierr);

	/* lb = 0 */
	ierr = VecDuplicate(*b,lb);CHKERRV(ierr);
	ierr = PetscObjectSetName((PetscObject)*lb,"lower bound");CHKERRV(ierr);
	
//	ierr = VecView(*lb, PETSC_VIEWER_STDOUT_WORLD);CHKERRV(ierr);

	/* B = [1 1 1, 0 0 0; 0 0 0 1 1 1]; */
	ierr = MatCreateDense(PETSC_COMM_WORLD,1,m,1,m,NULL, B);CHKERRV(ierr);
	ierr = PetscObjectSetName((PetscObject)*B,"matrix of EQ");CHKERRV(ierr);

	PetscScalar *B_data;
	ierr = MatDenseGetArray(*B,&B_data);CHKERRV(ierr);
	for(int i=0;i<m;i++){
		B_data[i] = 1.0;
	}
	ierr = MatDenseRestoreArray(*B,&B_data);CHKERRV(ierr);
	ierr = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
	ierr = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
	
//	ierr = MatView(*B, PETSC_VIEWER_STDOUT_WORLD);CHKERRV(ierr);

	/* c = 1 */
	ierr = VecCreate(PETSC_COMM_WORLD, c);CHKERRV(ierr);
	ierr = VecSetSizes(*c,1,1);CHKERRV(ierr);
	ierr = VecSetFromOptions(*c);CHKERRV(ierr);
	ierr = VecSet(*c,1.0);CHKERRV(ierr);
	ierr = VecAssemblyBegin(*c);CHKERRV(ierr);
	ierr = VecAssemblyEnd(*c);CHKERRV(ierr);
	ierr = PetscObjectSetName((PetscObject)*c,"rhs of EQ");CHKERRV(ierr);

//	ierr = VecView(*c, PETSC_VIEWER_STDOUT_WORLD);CHKERRV(ierr);
}

int main( int argc,char **args )
{
	/* initialize Petsc */
	PermonInitialize(&argc,&args,(char *)0,help);

	/* prepare ranodm generator */
//	srand ( time ( NULL));
	srand ( 0 );

	/* error handler */
	PetscErrorCode ierr;
  
	/* say hello */
	ierr = PetscPrintf(PETSC_COMM_WORLD,"This is SMALLEST ENCLOSING BALL example\n");CHKERRQ(ierr);

	/* prepare variables */
	PetscInt d; /* dimension of points */
	PetscInt m; /* number of points */
	Vec y; /* solution vector */
	Mat C; /* matrix of points */
	Mat A; /* Hessian matrix */
	Vec b; /* linear term */
	Vec lb; /* lower bound (rigid obstacle - non-penetration condition) */ 
	Mat B; /* matrix of equality constraints */
	Vec c; /* rhs of equality constraints */


	/* load console parameters */
	PetscBool set_d; /* the variable was provided in console parameters or not */
	ierr = PetscOptionsGetInt(NULL,NULL,"-d",&d,&set_d);CHKERRQ(ierr);
	if(!set_d) d = PERMON_EXAMPLE_BALL_DEFAULT_D;

	PetscBool set_m; /* the variable was provided in console parameters or not */
	ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,&set_m);CHKERRQ(ierr);
	if(!set_m) m = PERMON_EXAMPLE_BALL_DEFAULT_M;
	
	/* print parameters loaded from console */
	ierr = PetscPrintf(PETSC_COMM_WORLD," d                    = %d\n", d);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," m                    = %d\n", m);CHKERRQ(ierr);

	/* assemble problem */
	assemble_problem(d, m, &C, &A, &b, &lb, &B, &c);

	/* test */
//	ierr = MatView(C, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
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

	/* compute the center */
	Vec p;
	ierr = VecCreate(PETSC_COMM_WORLD, &p);CHKERRQ(ierr);
	ierr = VecSetSizes(p,d,d);CHKERRQ(ierr);
	ierr = VecSetFromOptions(p);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(p);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(p);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)p,"center");CHKERRQ(ierr);
	
	ierr = MatMult(C,y,p);

	ierr = PetscPrintf(PETSC_COMM_WORLD," center               = [");CHKERRQ(ierr);
	PetscScalar *p_array;
	ierr = VecGetArray(p,&p_array);CHKERRQ(ierr);
	for(int i=0;i<d;i++){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"%f", p_array[i]);CHKERRQ(ierr);
		if(i < d-1){
			ierr = PetscPrintf(PETSC_COMM_WORLD,", ");CHKERRQ(ierr);
		}
	}	
	ierr = VecRestoreArray(p,&p_array);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"] (should be 1)\n");CHKERRQ(ierr);
	
	PetscScalar yAy, by;
	ierr = VecDot(p,p,&yAy);CHKERRQ(ierr);
	ierr = VecDot(b,y,&by);CHKERRQ(ierr);
	
	PetscScalar r = yAy - by;
	if(r < 0){
		r = sqrt(-r);
	} else {
		r = sqrt(r);
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD," radius               = %f (should be %f)\n",r, 1.0);CHKERRQ(ierr);


	/* destroy the mess */
	ierr = QPDestroy(&qp);CHKERRQ(ierr);
	ierr = QPSDestroy(&qps);CHKERRQ(ierr);

	ierr = VecDestroy(&p);CHKERRQ(ierr);
	ierr = MatDestroy(&C);CHKERRQ(ierr);
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


