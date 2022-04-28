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

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    int myid, numprocs ; 
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Status status;
    int numcol, numrow, row;
    double ans, x[numcol]; 
    double *buffer; 
    

    if(myid==0){
        int sender;
        double A[numrow][numcol],  y[numrow];
        

        //Initialize matrix 
        numrow = numprocs-1; 
        numcol = 10 ; 
        MPI_Bcast(&numrow, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&numcol, 1, MPI_INT, 0, MPI_COMM_WORLD);
        for (int i=0; i< numrow; i++){
            y[i]=0;
            for ( int j =0 ; j< numcol; j++){
                A[i][j] = 1; 
                x[j]=1;
            }
        }
        // //print matrix 
        // for(int i=0; i< numrow; i++){
        //     printf("A[%d]= ",i);
        //     for(int j=0; j<numcol; j++){
        //         printf("%f ", A[i][j]);
        //     }
        //     printf("\n");

        // }
        // for(int i=0; i<numcol; i++){
        //     printf("X[%d]= %f\n", i,  x[i]);
        // }

        int numsent = 0; 
        // 모든 프로세스에 x 전달 
        MPI_Bcast(&x[0], numcol, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // A 행렬의 i 번째 행을 i+1 번째 process 에 전달
        // tag = i 
        for (int i=0; i< numprocs-1 && i<numrow; i++){
            MPI_Send(&A[i][0], numcol, MPI_DOUBLE, i+1, i, MPI_COMM_WORLD);
            numsent++;
        }
    
        //계산 결과 수합 
        for(int i=0; i<numrow; i++){
            MPI_Recv(&ans, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            sender = status.MPI_SOURCE; 
            row = status.MPI_TAG; 
            y[row]=ans; 
            //numproc>numrow 일 경우 계산 끝난 process 에 다음 행을 전달해주기! 
            if(numsent<numrow){
                MPI_Send(&A[numsent][0], numcol, MPI_DOUBLE, sender, numsent, MPI_COMM_WORLD);
                numsent ++;
            }
            // MPI_BOTTOM : indicate the bottom of the address space 
            else MPI_Send(MPI_BOTTOM, 0, MPI_DOUBLE, sender, numrow, MPI_COMM_WORLD);
        }

        //계산 결과 프린트 
        for(int i=0; i<numrow; i++){
            printf("result = %f\n", y[i]);
        }
    }
    else{ 
        
        MPI_Bcast(&numrow, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&numcol, 1, MPI_INT, 0, MPI_COMM_WORLD);
        double x[numcol];
        // x 수령 
        MPI_Bcast(&x[0], numcol, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        buffer = (double*)malloc(sizeof(double)*numcol);
        int po=1; // calculation is end or not 
        while(myid <= numrow && po){
            // buffer 에 A[i] 저장 
            MPI_Recv(buffer, numcol, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            row = status.MPI_TAG; 
            if (row < numrow){
                ans = 0.0; 
                // ans = A[i] dot x
                for(int j=0; j<numcol; j++){
                    ans += buffer [j] * x[j];
                }
                // 계산 결과 하나를 0번 프로세스에 보내기, tag = row 
                MPI_Send(&ans, 1, MPI_DOUBLE, 0, row, MPI_COMM_WORLD);
            }
            // 다 계산했으면 break!
            else po=0;
        }
        free(buffer);
    }
    
    MPI_Finalize();

}
