/*
Minyoung Kim (Apr. 28, 2022)
2022 Spring PSC class 
Subject: MPI
MPI software: mpich  
compile:mpicc -g ${file} -o ${fileDirname}/${fileBasenameNoExtension}.exe
*/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>


int main(int argc, char* argv[])
{
    //Initialize mpi 
    MPI_Init(&argc, &argv);
 
    // Get my rank and number of processors 
    int  rank, size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    if(size!=3){
        printf("You have to run this code with 3 processors\n");
        MPI_Abort(comm, EXIT_FAILURE);
    }
    if(rank == 0){
        
        int value = 100; 

        // Get 1, 1, 2 values from 3 processors
        int count_array[3]={1,1,2};
        // Receive data at [1], [3], [5]
        int displacements[3] = {1,3,5};
        //Determine buffer size which should be bigger than 5+2
        int n=displacements[2]+count_array[2]; 
        int *buffer =(int*)malloc(n * sizeof(int));
        printf("[Processor %d] my value = %d\n", rank, value);
        // gather value with receive buffer, receive count array, displacements 
        MPI_Gatherv(&value, 1, MPI_INT, buffer, count_array, displacements, MPI_INT, 0, comm);
        printf("[Processor %d] Value gathered\n", rank);
        for(int i = 0; i<n; i++){
            printf("[Processor %d] value %d = %d\n", rank, i, buffer[i]);
        }
        free(buffer);

    }else if(rank == 1){
        // Define my value
        int value = 101;
        printf("[Processor %d] my value = %d\n", rank, value);
        // Gather code to send value with NULL arguments for receive buffer, receive count array, displacements 
        MPI_Gatherv(&value, 1, MPI_INT, NULL, NULL, NULL, MPI_INT, 0, comm);
    }else if(rank == 2){
        //Define my value. If I am second processor, prepare int array with size 2 
        int values[2] = {102, 103};
        printf("[Processor %d] my value = %d & %d\n", rank, values[0], values[1]);
        MPI_Gatherv(&values[0], 2, MPI_INT, NULL, NULL, NULL, MPI_INT, 0 ,comm);

    }
 
    MPI_Finalize();
}