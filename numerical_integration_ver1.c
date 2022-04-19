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
        
        if (rank == 0){
            y_now += y_rank_now ; 
            for(int rank_num = 1; rank_num < num_processor; rank_num ++){
                //tag should be non-negative integer 
               tag = rank_num; 
               double sum_rank ;
               // for each processor except rank=0, 
               // get the y_rank_now value which is integrated value for each portion (from bin_min to bin_max)
               MPI_Recv(&sum_rank, 1, MPI_DOUBLE, rank_num, tag, MPI_COMM_WORLD, &status);
               y_now += sum_rank ; 
           }
           fprintf(stdout, "Number of bins in total = %d,  y = %.8f\n", nbin*num_processor, y_now);

         }// rank != 0  -> they should send y_rank_now values to rank 0 
         else{
             tag = rank; 
             MPI_Send(&y_rank_now, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
         }
        //Update converged variable in root(rank = 0)
        if (rank ==0){
            double err = fabs(y_now - y_prev) / (fabs(y_now + y_prev)/2.0);
            if (err < tolerance){
                converged = 1 ; 
            }
            for(int rank_num = 1; rank_num < num_processor; rank_num ++){
                tag = num_processor + rank_num; 
                MPI_Send(&converged, 1, MPI_INT, rank_num, tag, MPI_COMM_WORLD);
            }
         }
        //Send message to each processor that if it is converged or not
        else{
            tag = num_processor + rank; 
            MPI_Recv(&converged, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        }
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