/*
Minyoung Kim (Apr. 27, 2022)
2022 Spring PSC class 
Subject: MPI
MPI softward: msmpi // msmpi use gcc to compile mpi code 
compile:mpicc -g ${file} -o ${fileDirname}/${fileBasenameNoExtension}.exe
*/

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{   
    // start mpi 
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
 
    // Time to wait before processing, in seconds
    const double waiting_time = 1 / (double) size;
    int message = 0;
    double start;
    double end;
 
    if(rank == 0)
    {        
        //Set start 
        start = MPI_Wtime();
 
        //Simulate latency 
        while(MPI_Wtime() - start < waiting_time)
        {
            
        }
        //time = start + waiting_time 

        //Send message to process 1 
        MPI_Send(&message, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);        
    }
    else
    {
        //Recieve message from previous processor 
        MPI_Recv(&message, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
 
        //Start checking time 
        start = MPI_Wtime();
 
        //Simulate the latency
        while(MPI_Wtime() - start < waiting_time)
        {
            // We keep looping until <waiting_time> seconds have elapsed
        }
 
        //Except final processor 
        if(rank != size - 1)
        {
            //Send message to next processor 
            MPI_Send(&message, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
        }
    }
 
    // Wait for the very last one to 
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
 
    printf("[MPI process %2d] time elapsed during the job: %.2fs.\n", rank, end - start);
 
    MPI_Finalize();
 
}