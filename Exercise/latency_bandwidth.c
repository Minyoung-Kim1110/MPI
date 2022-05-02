/*
Minyoung Kim (Apr. 27, 2022)
2022 Spring PSC class 
Subject: MPI
MPI software: mpich  
compile:mpicc -g ${file} -o ${fileDirname}/${fileBasenameNoExtension}.exe
*/

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{   
    MPI_Init(&argc, &argv);
    int rank, size; 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int tag = 0 ;
    MPI_Status status;
    int count, array[100000];
    int n = 1000; 
    double start, end ; 
    double time[n];

    if (rank == 0){
        // If I am the first process, send message with increasing size to next processor 
        for (count=1; count< n; count ++){
            start = MPI_Wtime(); // in second
            MPI_Send(&array[0], count, MPI_INT, 1, count, MPI_COMM_WORLD);
            //Receive the message from next processor 
            MPI_Recv(&array[0], count, MPI_INT, 1, count, MPI_COMM_WORLD, &status);
            end = MPI_Wtime();
            time[count-1] = end - start; //elapsed time for ping-pong
            printf("Message size = [%3d], elapsed time = [%0.5f] ms\n", count, time[count-1]*1000.0);
        }
    }else if(rank == 1){
        //If I am the second processor, receive message from the first processor and send back. 
        //With increasing size of message 
        for(count=1; count < n; count ++){
            MPI_Recv(&array[0], count, MPI_INT, 0, count, MPI_COMM_WORLD, &status);
            MPI_Send(&array[0], count, MPI_INT, 0, count, MPI_COMM_WORLD);
        } 
    }else{}

    MPI_Finalize();
}
