/*
Minyoung Kim (Apr. 19, 2022)
2022 Spring PSC class 
Subject: MPI
MPI softward: mpich
compile:mpicc -g ${file} -o ${fileDirname}/${fileBasenameNoExtension}.exe
*/
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

//개인 함수 만들기 
void my_function(int *in, int *out, int *len, MPI_Datatype *dt){
    for (int i=0; i<*len; i++){
        out[i]+= in[i];
    }
}

int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);

    // 만든 함수를 사용자 정의 operator 에 저장 
    MPI_Op myop;
    MPI_Op_create((MPI_User_function *) my_function, 1, &myop);

    int rank, num_processor;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processor);
    int input[num_processor];

    if(rank == 0){
        // At root, input={100,100,...}
        for(int i =0; i< num_processor; i++){
            input[i]=100;
        }
    }
    else{
        for(int i =0; i< num_processor; i++){
            input[i]=0;
        }input[rank]+=1;
    }
    int out[100];
    // output array is initialized with input array in root  
    MPI_Reduce(&input, &out, num_processor,  MPI_INT, myop, 0, MPI_COMM_WORLD);
    if (rank == 0){
        for (int i=0; i<num_processor; i++){
            fprintf(stdout, "output value %d\n", out[i]);

        }

    }
    MPI_Finalize();
}
        