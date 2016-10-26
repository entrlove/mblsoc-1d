/*
 C++ header
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>
/*
 SLEPc header
*/
#include <slepceps.h>
/*
 Project header
*/
#include "dec.h"
int main(int argc,char **argv)
{
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                          Locals for Hamiltonian
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
     Dimension 
    */
    size_t site_num             = 8;
    size_t reduce_num           = 4;
    size_t remain_num           = 4;
    size_t particle_num         = 4;
    /* 
     Kinetic and soc
    */
    double t_1                  = 1.0;
    double t_2                  = 0.5;
    double lambda_z_1           = 1.0;
    double lambda_z_2           = 0.5;
    double lambda_y_1           = 1.0;
    double lambda_y_2           = 0.5;
    size_t hopping_distance_1   = 1;
    size_t hopping_distance_2   = 2;
    bool   boundary             = 1;
    /*
     Interaction and magnetic field
    */
    double hz                   = 0.0;
    double U                    = 5.0;
    double W                    = 40.0;
    double hb                   = 0.1;
    double ub                   = 0.1;
    int    seed                 = 2;
    /* 
     Disorder realization
    */
    std::vector<double> dis(site_num);
    /*
     Space and subspace
    */
    std::vector<dim_soc_subspace>     dim;
    std::vector<basis_soc_with_index> space;
    std::vector<size_t>               position_reduce(reduce_num);
    std::vector<size_t>               position_remain(remain_num);
    /*
     Set the reduced and remained position
    */
    for (size_t i = 0; i < reduce_num; i++) position_reduce[i] = i;
    for (size_t i = 0; i < remain_num; i++) position_remain[i] = i + reduce_num;
    /*
     Generate the space for Hamiltonian
    */
    Space_soc_index(site_num, particle_num, position_reduce, position_remain, dim, space);
    
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                         Locals for SLEPc and PETSc
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    Mat            A;                    /* problem matrix */
    EPS            eps_lm,eps_sm,eps_tm; /* eigenproblem solver context */
    PetscScalar    kr,ki;
    PetscReal      tol;
    Vec            xr,xi;
    PetscInt       n=space.size(),Istart,Iend,nev=50,its,maxit,nconv;
    PetscInt       d_nz=particle_num*(4+4+1),o_nz=particle_num*(4+4+1);
    PetscErrorCode ierr;
    ST             st;
    KSP            ksp;
    PC             pc;
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                               Locals for MPI
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    int            comm_sz, my_rank;
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                               Locals for spectrum
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    double         eval_max,eval_min,eval_target,pos=0.5,r_avg;
    std::vector<double> eval;
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                         Locals for entanglement entropy
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    Vec            evec_wr;
    VecScatter     scatter;
    IS             from_w,to_w,from_n,to_n;
    int            *idx_from,*idx_to;
    PetscScalar    *pevec_wr;
    double         ee,ee_avg;
    std::ofstream  output;
    std::vector<std::complex<double>> evec;
    
    /*
     Initialize the SLEPc context
    */
    SlepcInitialize(&argc,&argv,NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    /*
     Print information from the first process
    */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSLEPc initialized with %d processes\n", comm_sz);CHKERRQ(ierr);
    
    for (W = 0; W < 40.01; W += 2)
    {
        for (pos = 0.1; pos < 0.99; pos += 0.1)
        {
            for (seed = 1; seed < 201; seed++)
            {
                /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                 Compute the operator matrix that defines the eigensystem, Ax=kx
                 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
                ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
                ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
                ierr = MatSetType(A,MATMPIAIJ);CHKERRQ(ierr);
                ierr = MatMPIAIJSetPreallocation(A,d_nz,NULL,o_nz,NULL);CHKERRQ(ierr);

                ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
               
                /*
                 Set the matrix elements
                */
                Hamiltonian_hd(space,              /* Space for Hamiltonian */
                               position_reduce,    /* Reduced position */
                               position_remain,    /* Remained position */
                               A,                  /* Hamiltonian */
                               Istart,             /* Starting row of the Hamiltonian for the process */
                               Iend,               /* Ending (+1) row of the Hamiltonian for the process */
                               site_num,           /* Site number */
                               particle_num,       /* Particle number */
                               hopping_distance_1, /* Distance for nearest hopping */
                               t_1,                /* Kinetic strength for nearest hopping */
                               lambda_y_1,         /* Dresselhaus soc strength for nearest hopping */
                               boundary            /* Boundary condition */
                               );
                Hamiltonian_hd(space,              /* Space for Hamiltonian */
                               position_reduce,    /* Reduced position */
                               position_remain,    /* Remained position */
                               A,                  /* Hamiltonian */
                               Istart,             /* Starting row of the Hamiltonian for the process */
                               Iend,               /* Ending (+1) row of the Hamiltonian for the process */
                               site_num,           /* Site number */
                               particle_num,       /* Particle number */
                               hopping_distance_2, /* Distance for next-nearest hopping */
                               t_2,                /* Kinetic strength for next-nearest hopping */
                               lambda_y_2,         /* Dresselhaus soc strength for next-nearest hopping */
                               boundary            /* Boundary condition */
                               );
                Hamiltonian_sr(space,              /* Space for Hamiltonian */
                               position_reduce,    /* Reduced position */
                               position_remain,    /* Remained position */
                               A,                  /* Hamiltonian */
                               Istart,             /* Starting row of the Hamiltonian for the process */
                               Iend,               /* Ending (+1) row of the Hamiltonian for the process */
                               site_num,           /* Site number */
                               particle_num,       /* Particle number */
                               hopping_distance_1, /* Distance for nearest hopping */
                               lambda_z_1,         /* Rashba soc strength for nearest hopping */
                               boundary            /* Boundary condition */
                               );
                Hamiltonian_sr(space,              /* Space for Hamiltonian */
                               position_reduce,    /* Reduced position */
                               position_remain,    /* Remained position */
                               A,                  /* Hamiltonian */
                               Istart,             /* Starting row of the Hamiltonian for the process */
                               Iend,               /* Ending (+1) row of the Hamiltonian for the process */
                               site_num,           /* Site number */
                               particle_num,       /* Particle number */
                               hopping_distance_2, /* Distance for next-nearest hopping */
                               lambda_z_2,         /* Rashba soc strength for next-nearest hopping */
                               boundary            /* Boundary condition */
                               );
                srand(seed);
                for (size_t i = 0; i < site_num; i++) dis[i] = W*((double)rand()/RAND_MAX - 0.5);
                Hamiltonian_dg(space,              /* Space for Hamiltonian */
                               A,                  /* Hamiltonian */
                               Istart,             /* Starting row of the Hamiltonian for the process */
                               Iend,               /* Ending (+1) row of the Hamiltonian for the process */
                               site_num,           /* Site number */
                               dis,                /* Disorder realization */
                               hz,                 /* Out-of-plane Zeeman field */
                               U,                  /* Interaction strength */
                               hb,                 /* Magnetic field on left boundary */
                               ub                  /* Magnetic field on right boundary */
                               );
                /*
                 Assemble the matrix
                */
                ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
                ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

                ierr = MatCreateVecs(A,NULL,&xr);CHKERRQ(ierr);
                ierr = MatCreateVecs(A,NULL,&xi);CHKERRQ(ierr);

                
                /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                            Create the eigensolver and set various options
                 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
                /*
                 Create eigensolver context
                */
                ierr = EPSCreate(PETSC_COMM_WORLD,&eps_lm);CHKERRQ(ierr);
                ierr = EPSCreate(PETSC_COMM_WORLD,&eps_sm);CHKERRQ(ierr);
                ierr = EPSCreate(PETSC_COMM_WORLD,&eps_tm);CHKERRQ(ierr);
                /*
                 Set the common properties of operators.
                */
                ierr = EPSSetOperators(eps_lm,A,NULL);CHKERRQ(ierr);
                ierr = EPSSetOperators(eps_sm,A,NULL);CHKERRQ(ierr);
                ierr = EPSSetOperators(eps_tm,A,NULL);CHKERRQ(ierr);
                ierr = EPSSetProblemType(eps_lm,EPS_HEP);CHKERRQ(ierr);
                ierr = EPSSetProblemType(eps_sm,EPS_HEP);CHKERRQ(ierr);
                ierr = EPSSetProblemType(eps_tm,EPS_HEP);CHKERRQ(ierr);
                /*
                 Set the eigensolver context for compute the largest eigenvalue
                */
                ierr = EPSSetWhichEigenpairs(eps_lm,EPS_LARGEST_REAL);CHKERRQ(ierr);
                /*
                 Set the eigensolver context for compute the smallest eigenvalue
                 */
                ierr = EPSSetWhichEigenpairs(eps_sm,EPS_SMALLEST_REAL);CHKERRQ(ierr);
                /*
                 Set the eigensolver context for compute the target eigenvalue
                 */
                ierr = EPSSetDimensions(eps_tm,nev,PETSC_DEFAULT,PETSC_DEFAULT);
                ierr = EPSSetWhichEigenpairs(eps_tm,EPS_TARGET_MAGNITUDE);CHKERRQ(ierr);
                ierr = EPSGetST(eps_tm,&st);CHKERRQ(ierr);
                ierr = STGetKSP(st,&ksp);CHKERRQ(ierr);
                ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
                ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
                ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
                ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);CHKERRQ(ierr);
                ierr = STSetType(st,STSINVERT);CHKERRQ(ierr);
               
                /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                  Solve the eigensystem
                 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
                /*
                 Compute the largest and smallest eigenvalue
                */
                ierr = EPSSolve(eps_lm);CHKERRQ(ierr);
                ierr = EPSSolve(eps_sm);CHKERRQ(ierr);
                
                ierr = EPSGetConverged(eps_lm,&nconv);CHKERRQ(ierr);
                if (nconv>0) {
                    ierr = EPSGetEigenpair(eps_lm,0,&kr,&ki,NULL,NULL);CHKERRQ(ierr);
                    eval_max = PetscRealPart(kr);
                } else {
                    eval_max = std::numeric_limits<double>::quiet_NaN();
                }

                ierr = EPSGetConverged(eps_sm,&nconv);CHKERRQ(ierr);
                if (nconv>0) {
                    ierr = EPSGetEigenpair(eps_sm,0,&kr,&ki,NULL,NULL);CHKERRQ(ierr);
                    eval_min = PetscRealPart(kr);
                } else {
                    eval_min = std::numeric_limits<double>::quiet_NaN();
                }
                /*
                 Compute the target eigenvalues
                */
                nconv = 0;
                eval_target = pos*(eval_min-eval_max) + eval_max;
                /*
                 If target is not NaN, compute the targeted eigenvalues
                */
                if (std::isnan(eval_target) == 0)
                {
                    ierr = EPSSetTarget(eps_tm,eval_target);CHKERRQ(ierr);
                    ierr = EPSSolve(eps_tm);CHKERRQ(ierr);
                    ierr = EPSGetConverged(eps_tm,&nconv);CHKERRQ(ierr);
                }
                
                /*
                 Optional: Get some information from the solver
                 */
                ierr = EPSGetIterationNumber(eps_tm,&its);CHKERRQ(ierr);
                ierr = EPSGetDimensions(eps_tm,&nev,NULL,NULL);CHKERRQ(ierr);
                ierr = EPSGetTolerances(eps_tm,&tol,&maxit);CHKERRQ(ierr);
                
                MPI_Barrier(MPI_COMM_WORLD);
                /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                     Compute the average ratio
                 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
                eval.resize(nconv);
                for (int i=0;i<nconv;i++){
                    ierr = EPSGetEigenpair(eps_tm,i,&kr,&ki,NULL,NULL);CHKERRQ(ierr);
                    eval[i] = PetscRealPart(kr);
                }
                ravg(eval, r_avg);
                
                MPI_Barrier(MPI_COMM_WORLD);
                /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                    Compute the entangle entropy
                 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
                idx_to   = (int *)malloc(n*sizeof(int));
                idx_from = (int *)malloc(n*sizeof(int));
                /*
                 Initialize the indices for the source and destination
                 */
                for (int i = 0; i < n; i++) { idx_to[i] = idx_from[i] = i; }
                /*
                 Apply the storage for eivenvector
                 */
                VecCreateSeq(PETSC_COMM_SELF,n,&evec_wr);
                /*
                 Create the scatter for null index and whole index
                 */
                ISCreateGeneral(PETSC_COMM_SELF,0,idx_from,PETSC_COPY_VALUES,&from_n);
                ISCreateGeneral(PETSC_COMM_SELF,0,idx_to,PETSC_COPY_VALUES,&to_n);
                ISCreateGeneral(PETSC_COMM_SELF,n,idx_from,PETSC_COPY_VALUES,&from_w);
                ISCreateGeneral(PETSC_COMM_SELF,n,idx_to,PETSC_COPY_VALUES,&to_w);
                /*
                 Gather eigenvectors
                */
                evec.clear();
                for (int i=0;i<nconv;i++){
                    ierr = EPSGetEigenpair(eps_tm,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
                    /* 
                     Create appropriate scatter for i-th eigenvector on each process
                    */
                    if (i%comm_sz == my_rank) {
                        VecScatterCreate(xr,from_w,evec_wr,to_w,&scatter);
                    }
                    else {
                        VecScatterCreate(xr,from_n,evec_wr,to_n,&scatter);
                    }
                    /* 
                     Do gather
                    */
                    VecScatterBegin(scatter,xr,evec_wr,INSERT_VALUES,SCATTER_FORWARD);
                    VecScatterEnd(scatter,xr,evec_wr,INSERT_VALUES,SCATTER_FORWARD);
                    VecScatterDestroy(&scatter);
                    /*
                     Push the eigenvector into container
                     */
                    if (i%comm_sz == my_rank) {
                        VecGetArray(evec_wr, &pevec_wr);
                        evec.insert(evec.end(), pevec_wr, pevec_wr+n);
                        VecRestoreArray(evec_wr, &pevec_wr);
                     }
                }
                /*
                 Initialize the average value of entangle entropy
                 */
                ee_avg=0;
                /*
                 Sum entangle entropies on first process
                */
                if (my_rank!=0){
                    for (int i=0;i<evec.size()/n;i++) {
                        Entanglement_entropy(dim,&evec[i*n],ee);
                        MPI_Send(&ee,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
                    }
                } else {
                    for (int i=0;i<evec.size()/n;i++) {
                        Entanglement_entropy(dim,&evec[i*n],ee);
                        ee_avg+=ee;
                    }
                    for (int i=0;i<nconv-evec.size()/n;i++) {
                        MPI_Recv(&ee,1,MPI_DOUBLE,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                        ee_avg+=ee;
                    }
                    /*
                     If nconv equals to zero, the expression generates NaN
                    */
                    ee_avg /= nconv;
                    ee_avg /= site_num;
                }
                
                /*
                 Output results
                 */
                if (my_rank == 0)
                {
                    std::cout << seed << std::endl << pos << std::endl << W << std::endl;
                    std::cout << eval_min << std::endl << eval_max << std::endl;
                    std::cout << its << std::endl;
                    output.open("res.txt",std::ios::app);
                    output << std::setprecision(8);
                    output << std::setw(16) << std::fixed
                           << std::setw(16) << std::fixed << eval_min
                           << std::setw(16) << std::fixed << eval_max
                           << std::setw(16) << std::fixed << eval_target
                           << std::setw(16) << std::fixed << r_avg
                           << std::setw(16) << std::fixed << ee_avg
                           << std::setw(16) << std::fixed << its
                           << std::setw(16) << std::fixed << maxit
                           << std::setw(16) << std::fixed << tol
                           << std::endl;
                    output.close();
                }

                /*
                 Free work space
                */
                free(idx_to);
                free(idx_from);
                ierr = VecDestroy(&evec_wr);CHKERRQ(ierr);
                ierr = ISDestroy(&from_n);  CHKERRQ(ierr);
                ierr = ISDestroy(&from_w);  CHKERRQ(ierr);
                ierr = ISDestroy(&to_n);    CHKERRQ(ierr);
                ierr = ISDestroy(&to_w);    CHKERRQ(ierr);
                
                ierr = EPSDestroy(&eps_lm); CHKERRQ(ierr);
                ierr = EPSDestroy(&eps_sm); CHKERRQ(ierr);
                ierr = EPSDestroy(&eps_tm); CHKERRQ(ierr);
                ierr = MatDestroy(&A);      CHKERRQ(ierr);
                ierr = VecDestroy(&xr);     CHKERRQ(ierr);
                ierr = VecDestroy(&xi);     CHKERRQ(ierr);
            }
        }
    }
    ierr = SlepcFinalize();
    return ierr;
}
