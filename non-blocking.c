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
    
    

    MPI_Finalize();
}