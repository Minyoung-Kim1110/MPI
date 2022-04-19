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

//개인 함수 만들기 
void my_function(int *in, int *out, int *len, MPI_Datatype *dt){
    *out = 0 ; 
    for (int i=0; i<*len; i++){
        *out+= in[i]+1;
    }
    printf("%d\n", *out);
}

int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);

    // 만든 함수를 사용자 정의 operator 에 저장 
    MPI_Op myop;
    MPI_Op_create((MPI_User_function *) my_function, 1, &myop);

    int rank, num_processor;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processor);

    int out[num_processor];
    MPI_Reduce(&rank, &out, 1,  MPI_INT, myop, 0, MPI_COMM_WORLD);
    if (rank == 0){
        for (int i=0; i<num_processor; i++){
            fprintf(stdout, "output value %d\n", out[i]);

        }

    }
    


    MPI_Finalize();
}
        