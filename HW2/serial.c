/*
Minyoung Kim (May. 1, 2022)
compile:gcc -g ${file} -o ${fileDirname}/${fileBasenameNoExtension}.exe
execution: ./${fileBasenameNoExtension}.exe

This file contains generating matrix with serial execution 
Get n by terminal save matrix in Matrix_serial_dim={n}.txt file 

For n=2, 5, these results are used to test parallel version 
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

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

//Print matrix in file  
void print_matrix(FILE* fp ,double** A, int n){
    for(int i=1; i<n+1; ++i){
        for(int j=1; j<n+1; ++j){
            fprintf(fp, "%.3f ", A[i-1][j-1]);
        }//print element with same row in the same line
        fprintf(fp, "\n");
    }
}
int main(){
    
    //Get matrix dimension
    int n;
    printf("Enter matrix dimension: ");
    scanf("%d", &n);

    //Generate matrix
    double** mat = generate_matrix(n);

    //Save in file 
    char filename[100] ;
    sprintf(filename, "Matrix_serial_dim=%d.txt", n);
    FILE* fp = fopen(filename, "w");
    print_matrix(fp, mat, n);
    
    //Transpose matrix
    mat = transpose(mat, n);

    //write indicator line 
    fprintf(fp, "\nTransposed\n\n");

    //Save in file 
    print_matrix(fp, mat, n);
    fclose(fp);
}
// int main(){
//     //Get matrix dimension
//     int dims[8] = {2,4,8,16,32,64,128,256};
//     int n;
//     clock_t start, end;
//     for(int i=0; i<8; i++){
//         n=dims[i];
//         start = clock(); 
//         //Generate matrix
//         double** mat = generate_matrix(n);

//         //Save in file 
//         char filename[100] ;
//         sprintf(filename, "Matrix_serial_dim=%d.txt", n);
//         FILE* fp = fopen(filename, "w");
//         print_matrix(fp, mat, n);

//         //Transpose matrix
//         mat = transpose(mat, n);

//         //write indicator line 
//         fprintf(fp, "\nTransposed\n\n");

//         //Save in file 
//         print_matrix(fp, mat, n);
//         fclose(fp);
//         end = clock();
//         double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
//         printf("DIM =%d, time spent %f\n", n, time_spent);
//     }
// }
