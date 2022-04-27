/*
Minyoung Kim (Apr. 27, 2022)
2022 Spring PSC class 
Subject: MPI
MPI softward: msmpi // msmpi use gcc to compile mpi code 
compile:mpicc -g ${file} -o ${fileDirname}/${fileBasenameNoExtension}.exe
*/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    int rank, size;
    // 정보 저장하는 곳 
    MPI_Status status;
    // 현재 프로세스 ID
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // 전체 프로세스 개수 
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Request request; 
    if (rank == 0 ){
        int message = 10 ; 
        double start = MPI_Wtime();
        double waiting_time = 0.5; //waiting time in seconds
        while (MPI_Wtime()- start < waiting_time){
            //If I am first processor, wait for a while.. 
        }
        // If I am the first processor, send message to second processor
        MPI_Isend(&message, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &request);

    }
        
    else if (rank ==1){
        int message = 0 ; 
        //If I am second processor, I will receive message. 
        // Before receive, message = 0 
        printf("Before receive message = %d\n", message);
        // Receive message with non- blocking communication 
        MPI_Irecv(&message,1,MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
        // message may not received but next codes are executed! Hence, message = 0 
        printf("Receive code executed but non blocking message = %d\n", message);
        // Wait until message is received
        MPI_Wait(&request, &status);
        // message has been received -> message = 10 
        printf("After blocking message = %d\n", message);

    }

    MPI_Finalize();
}