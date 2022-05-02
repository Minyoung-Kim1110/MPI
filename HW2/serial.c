/*
Minyoung Kim (May. 1, 2022)
Subject: MPI
MPI software: mpich  
compile:gcc -g ${file} -o ${fileDirname}/${fileBasenameNoExtension}.exe
execution: ./${fileBasenameNoExtension}.exe
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

//Give the value of A_{ij} for fixed n where i,j = 1, ..., n 
double element(int i, int j, int n){
    //Initalize s with 0 where s = A_{ij}
    double s = 0.0;
    //Add with loop of k 
    for (int k=1; k<n+1; k++){
        //l=0
        s+=1/(sqrt(pow((double)(i-k), 2) + (double) pow((double) j, 2)));
        // |l|= 1, ... , 5n
        for(int l=1; l<5*n+1; l++){
            //There's no condition for j == -l since j is positive. Hence, always add
            s+=1/(sqrt(pow((double)(i-k), 2) + (double) pow((double) (j+l), 2)));
            if( (i-k)==0 && (j-l)==0){ //If denominator is 0 
                //skip execution 
                continue;
            }
            else{
                //add for l 
                s+=1/(sqrt(pow((double)(i-k), 2) + (double) pow((double) (j-l), 2)));
            }
            
        }
    }
    return s; 
}

//Generate n by n matrix A using pointer to pointer
double** generate_matrix(int n){
    //generate pointer to pointer to allocate 2D matrix
    double** A= (double**)malloc(n * sizeof(double*));
    for(int i=0; i<n; i++){
        //allocate for each row
        A[i] = (double*)malloc(n * sizeof(double));
    }
    //for each element, generate element by A_element function 
    for(int i=1; i<n+1; ++i){
        for(int j=1; j<n+1; ++j){
            printf("(%d, %d) =%f\n", i, j, A_element(i, j, n));
            //matrix index should start from 0 hence cast i,j to i-1, j-1
            A[i-1][j-1] = element(i, j, n);
        }
    }
    return A; 
}

//Transpose matrix 
double** transpose(double** mat, int n){
    //generate pointer to pointer to allocate 2D matrix
    double** tr= (double**)malloc(n * sizeof(double*));
    for(int i=0; i<n; i++){
        tr[i] = (double*)malloc(n * sizeof(double));
    }
    //Transpose by flipping index 
    for(int i=0; i<n; ++i){
        for(int j=0; j<n; ++j){
            tr[j][i]=mat[i][j];
        }
    }
    return tr;
}

//Print matrix for debugging 
void print_matrix(double** A, int n){
    for(int i=1; i<n+1; ++i){
        for(int j=1; j<n+1; ++j){
            printf("%.3f ", A[i-1][j-1]);
        }//print element with same row in the same line
        printf("\n");
    }
}
int main(){
    int n=2;
    double** mat = generate_matrix(n);
    print_matrix(mat, n);
    mat = transpose(mat, n);
    printf("matrix transposed!\n");
    print_matrix(mat, n);
}