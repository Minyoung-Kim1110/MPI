/*
Minyoung Kim (Apr. 19, 2022)
2022 Spring PSC class 
Subject: MPI
MPI softward: msmpi // msmpi use gcc to compile mpi code 
compile:mpicc -g ${file} -o ${fileDirname}/${fileBasenameNoExtension}.exe
*/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

double func_integrand(double x){
    return 1/cosh(x);
}

int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);

    int rank, num_processor;

    MPI_Comm_size(MPI_COMM_WORLD, &num_processor);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int tag;
    MPI_Status status; // Structure to store information 
    double a=0.0, b=1.0; //integration from 0 to 1  
    double tolerance = 1e-4;
    double y_now = 0.0; 
    double y_prev ; 
    double y_rank_now =0.0 ;
    double y_rank_prev ; 
    int converged = 0; // not converged = 0, if converged 1 

    // number of bin 
    int nbin = 1; 
    //split the bin for each processor 
    double bin_min = a + (b-a) * (double) rank / (double) num_processor;
    double bin_max = bin_min +(b-a)/(double) num_processor;
    double delta = fabs(bin_max - bin_min)/(double) nbin; // fabs = float number absolute value 
    while (converged ==0){
        y_prev = y_now; 
        y_rank_prev = y_rank_now;  
        y_now = 0; 
        y_rank_now = 0 ; 
        for(int i =0; i<nbin; i++){
            double start = bin_min + delta * (double) i; 
            double end = start + delta ; 
            double mean = (func_integrand(start) + func_integrand(end))/2.0;
            y_rank_now += delta * mean ; 
        }
        // operation for each y_rank_now and store the result in y_now 
        // type of y_rank_now is MPI_DOUBLE, save y_now value in rank=0
        // In here, operation is MPI_SUM = summation 
        MPI_Reduce(&y_rank_now, &y_now, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank ==0){
            fprintf(stdout, "Number of bins in total = %d,  y = %.12f\n", nbin*num_processor, y_now);
            double err = fabs(y_now - y_prev) / (fabs(y_now + y_prev)/2.0);
            if (err < tolerance){
                converged = 1 ; 
            }
         }
        //Broadcast converged value from rank 0 to other processors in MPI_COMM_WORLD
        MPI_Bcast(&converged, 1, MPI_INT, 0, MPI_COMM_WORLD);
        nbin *=2;
        delta /=2;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0 ) {
        fprintf(stdout, "Numerical integrated value = %.8f\n", y_now);
        fprintf(stdout, "It should be %.8f\n", 2.0 * (atan(exp(1))-atan(exp(0))));
    }
    //Finalize mpi 
    MPI_Finalize();
    return 0;
}