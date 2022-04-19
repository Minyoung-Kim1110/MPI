/*
Minyoung Kim (Apr. 19, 2022)
2022 Spring PSC class 
Subject: MPI
MPI softward: msmpi // msmpi use gcc to compile mpi code 
compile:
*/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>



int main(int argc, char *argv[]){

    int processor_size, rank; 
    // MPI 초기화 
    MPI_Init(&argc, &argv);

    // MPI 에서 사용하는 프로세스와 각 프로세스의 랭크를 정수형 인자의 포인터를 넘겨줘서 거기에 저장함. 
    //MPI_COMM_WORLD : 커뮤니케이터, 모든 프로세스의 집합 
    MPI_Comm_size(MPI_COMM_WORLD, &processor_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank ==0){
        fprintf(stdout, "Processor num %d\n", processor_size);

    }

    MPI_Barrier(MPI_COMM_WORLD);

    fprintf(stdout, "[PROCESSOR %d]: Hello World! \n", rank);

    MPI_Finalize();
 

}