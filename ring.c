/*
Minyoung Kim (Apr. 27, 2022)
2022 Spring PSC class 
Subject: MPI
MPI software: mpich  
compile:mpicc -g ${file} -o ${fileDirname}/${fileBasenameNoExtension}.exe
*/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    int rank, size;
    int baton;
    // 정보 저장하는 곳 
    MPI_Status status;
    // 현재 프로세스 ID
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // 전체 프로세스 개수 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    

    if (rank==0){
        baton = 1; 
        // baton 1개, MPI_INT 타입을 1번 프로세스에 tag = 999 로 전달 
        MPI_Send(&baton, 1, MPI_INT, 1, 999, MPI_COMM_WORLD);
        //baton 1개, MPI_INT 타입을 size-1 번 프로세스에서 받으려고 대기 
        MPI_Recv(&baton, 1, MPI_INT, size - 1, 999, MPI_COMM_WORLD, &status);
        printf("Baton %d: Process %d --> Process 0\n", baton, size-1);
    }else{
        // baton 받고 다음 프로세스에 던지기 
        MPI_Recv(&baton, 1, MPI_INT, rank-1, 999, MPI_COMM_WORLD, &status);
        printf("Baton %d: Process %d --> Process %d\n", baton, rank-1, rank);
        MPI_Send(&baton, 1, MPI_INT, (rank+1)%size, 999, MPI_COMM_WORLD);
    }

    MPI_Finalize();
}