/*
Minyoung Kim (May. 26, 2022)
2022 Spring PSC class 
Subject: MPI
MPI softward: mpich
compile:mpicc -g ${file} -o ${fileDirname}/${fileBasenameNoExtension}.exe -lm
*/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int main(int argc, char*argv[]){
    MPI_Init(&argc, &argv);
    int rank, size; 
    // 현재 프로세스 ID
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // 전체 프로세스 개수 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    //If I am root, become a sender
    if (rank ==0 ){
        int message; 
        MPI_Request request;

        //Prepare send request 
        MPI_Send_init(&message, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &request);
        for(int i=0; i<3; i++){
            message = i*100; 
            //Actual send starts 
            MPI_Start(&request);
            //Wait for the send to complete 
            MPI_Wait(&request, MPI_STATUS_IGNORE);
            printf("[Processor %d]: Send message %d\n", rank, message);
        }
    }else if(rank == 1){
        //If I am rank 2, become a receiver 
        int received;
        for(int i = 0; i < 3; i++){
            MPI_Recv(&received, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("[Processor %d]: Got message %d\n", rank, received);
        }
    }
    MPI_Finalize();
}
