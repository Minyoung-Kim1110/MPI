/*
Minyoung Kim (May. 10, 2022)
Subject: MPI
MPI software: mpich  
compile:mpicc -g ${file} -o ${fileDirname}/${fileBasenameNoExtension}.exe
execution: mpiexec -n {number of processors} ./${fileBasenameNoExtension}.exe

In this code, I want to generate A_ij = sum_{k,l} {(i-k)^2 + (j-l)^2}^{-1/2}
where i, j, k = 1, 2, ... n and l = -5n, -5n+1, ..., 5n
and transpose it using MPI coolective communication(explicitly MPI_Alltoall)

If number of processors is bigger than matrix dimension n, use only n processors 
to generate and transpose matrix. 
Else if the number of processors d is smaller than matrix dimension n, use all processors by 
making d by d block matrix. 
Each processor contains a row of d blocks and each block has d/n elements. (if d%n is not 0, than d/n +1)

Use MPI_Alltoall to transpose block matrix and in each processor, transpose elements in block. 

I used scanf function to get matrix dimension n. 
Also, after generating or transposing matrix, collect all data in root to save in a single file. 

*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

//Function for generating matrix 
double element(int i, int j, int n);

//Function for transposing small matrix 
MPI_Comm generate_new_communicator(int n);
void print_vector_in_file(FILE* fp, double* vec, int n);
void print_with_multiprocessor(FILE* fp, double* row, int rank, int n, MPI_Comm new_comm);
void make_and_transpose_small(FILE* fp, int rank, int n);

//Function for transposing big matrix 
void print_row_in_file(FILE* fp, int rank, double* row, int num_processors, int quotient, int n);
void print_with_multiprocessor_block(FILE* fp, double* row, int rank, int num_processors, int quotient, int n);
double* generate_row_with_block(double* row, int rank, int num_processors, int quotient, int n);
double* transpose_block(double* transposed_row, int block_index, int quotient);
void make_and_transpose_big(FILE* fp, int n, int num_processors, int rank);

//Function for test 
int test_element();
int compare_two_txt(char* filename1, char* filename2);
int test_print_with_multiprocessor(int rank, int num_processors);
int test_make_and_transpose_small(int rank, int num_processors);
int test_transpose_block();
int test_make_and_transpose_big(int rank, int num_processors);


//Give the value of A_{ij} for fixed n where i, j = 1, ..., n 
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
//Generate new communicator to transpose small matrix
MPI_Comm generate_new_communicator(int n){
    //Get the group of the default communicator(MPI_COMM_WORLD)
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    
    //Generate new group 
    MPI_Group new_group;
    // ranks: rank of processors to be included in new communicator
    int *ranks=(int*) malloc(n * sizeof(int));
    for(int i=0; i<n; i++){
        ranks[i]=i;
    }
    //Generate new group that contains n processors
    MPI_Group_incl(world_group, n, ranks, &new_group);

    //Generate new communicator from new group 
    MPI_Comm new_communicator; 
    MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_communicator);
    return new_communicator;
}
//print function that prints a row in a line
void print_vector_in_file(FILE* fp, double* vec, int n){
    for(int i=0; i<n; i++){
        char buf[sizeof(double)];
        fprintf(fp, "%.3f ", vec[i]);
    }fprintf(fp, "\n");
}
//print each row in order
void print_with_multiprocessor(FILE* fp, double* row, int rank, int n, MPI_Comm new_comm){
    MPI_Status status;
    // Collect data in root to save in txt file 
    if (rank==0){
        print_vector_in_file(fp, row, n);
        //allocate memory to save data temporary(just for print)
        double* row_tmp= (double*)malloc(n*sizeof(double));
        //recieve data in order 
        for(int i=1; i<n; i++){
            MPI_Recv(row_tmp, n, MPI_DOUBLE, i, i, new_comm, &status);
            //print! i th row is printed at i th turn
            print_vector_in_file(fp, row_tmp, n);     
        } 
        free(row_tmp);
    }else if(rank<n){
        //Send data to root
        MPI_Send(row, n, MPI_DOUBLE, 0, rank, new_comm);        
    }
}
//If matrix dimension n is less than processor number, use this function to make and transpose matrix 
void make_and_transpose_small(FILE* fp, int rank, int n){  
    //Generate communicator with n processors
    MPI_Comm new_comm = generate_new_communicator(n);
    //Generate a single row for each processor 
    if(rank<n){
        //Generate pointer to allocate a row and transposed data
        double* row = (double*) malloc(n*sizeof(double));
        double* transposed_row = (double*) malloc(n*sizeof(double));
        // For each rank(0, ... , n-1) row_index = (1, 2, ..., n)
        int row_index = rank +1;  
        //From 1 to n generate each element
        for(int j=1; j<n+1; ++j){
            row[j-1] = element(row_index, j, n);
        }

        //Save data in txt file 
        print_with_multiprocessor(fp, row, rank, n, new_comm);
        if(rank ==0 ){
            //write indicator line 
            fprintf(fp, "\nTransposed\n\n");
        }
        //Transpose Matrix  
        MPI_Alltoall(row, 1, MPI_DOUBLE, transposed_row, 1, MPI_DOUBLE, new_comm);
        //Save data in txt file 
        print_with_multiprocessor(fp, transposed_row, rank, n, new_comm);
        
        //free data
        free(row);
        free(transposed_row);
    }
}

//Print function for block matrix 
void print_row_in_file(FILE* fp, int rank, double* row, int num_processors, int quotient, int n){
    for(int i=0; i<quotient; i++){
        for (int block_index = 0; block_index<num_processors; block_index++){
            for(int j=0; j<quotient; j++){
                //If the element is matrix element not padded then save in file
                if(rank*quotient+i<n && block_index*quotient+j<n )
                    fprintf(fp, "%.3f ", row[block_index*quotient*quotient+i*quotient+j]);
            }
        }
        fprintf(fp, "\n");
    }
}
//revised version of print_with_multiprocessor with print_row_in_file function 
void print_with_multiprocessor_block(FILE* fp, double* row, int rank, int num_processors, int quotient, int n){
    
    MPI_Status status;
    int length = num_processors*quotient*quotient;
    // Collect data in root to save in txt file 
    if (rank==0){
        print_row_in_file(fp, rank, row, num_processors, quotient, n);
        //allocate memory to save data temporary(just for print)
        double* row_tmp= (double*)malloc(length*sizeof(double));
        //recieve data in order 
        for(int i=1; i<num_processors; i++){
            MPI_Recv(row_tmp, length, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
            //print! i th row is printed at i th turn
            print_row_in_file(fp, i,  row_tmp, num_processors, quotient, n);
        
        } 
        free(row_tmp);
    }else if(rank<n){
        //Send data to root
        MPI_Send(row, length, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);        
    }

}
//Generate elements obeying the rule and save data as follow
double* generate_row_with_block(double* row, int rank, int num_processors, int quotient, int n){
    //each row has part of matrix(A). row=A[rank*quotient:(rank+1)*quotient, 0:quotient*num_processor]
    //Save data as follow
    //row = [block_1_[0,0], ... ,block_1_[0, quotient-1], block_1_[1,0], ... , block_1_[quotient-1,quotient-1], 
    //       block_2_[0,0], ... ,block_2_[0, quotient-1], block_2_[1,0], ... , block_2_[quotient-1,quotient-1], 
    //       , ... ,
    //      block_(num_processors)[quotient-1, quotient-1]]
    int block_index, row_index=rank*quotient, column_index;
    // Generate elements until n(matrix dimension)
    
    //Generate elements following the rule
    for (block_index = 0; block_index<num_processors; block_index++){
        column_index = block_index*quotient;
        for(int i=0; i<quotient; i++){
            for(int j=0; j<quotient; j++){
                if(column_index+j+1<=n && row_index+i+1<=n){
                    row[block_index*quotient*quotient+i*quotient+j] = element(row_index+i+1, column_index+j+1, n);
                }
                //for index > n, pad 0.0
                else{
                    row[block_index*quotient*quotient+i*quotient+j] =0.0;
                }
            }
        } 
    }
    return row;
}

double* transpose_block(double* transposed_row, int block_index, int quotient){
    int start_pt = block_index * quotient*quotient;
    //generate tmp for temporary save transposed data 
    double* tmp = (double*)malloc(quotient*quotient*sizeof(double));
    for(int i=0; i<quotient; i++){
        for(int j=0; j<quotient; j++){
            tmp[j*quotient+i]=transposed_row[start_pt+i*quotient+j];
        }
    }
    //write in transposed row 
    for(int i=0; i<quotient; i++){
        for(int j=0; j<quotient; j++){
            transposed_row[start_pt+i*quotient+j]=tmp[i*quotient+j];
        }
    }
    return transposed_row;
}

//If n > number of processors, distribute task by making block
void make_and_transpose_big(FILE* fp, int n, int num_processors, int rank){
    //Calculate block size
    int quotient =  n/num_processors+ ((n % num_processors) != 0);
    //I want to make num_processors by num_processors matrix with blocks. 
    //Each block is quotient by quotient matrix 

    //Generate row of blocks for each processors. 
    double* row = (double*) malloc(num_processors*quotient*quotient*sizeof(double));
    
    //Generate elements 
    row = generate_row_with_block(row, rank, num_processors, quotient, n);

    //Save data in txt file 
    print_with_multiprocessor_block(fp, row, rank, num_processors, quotient, n);
    //allocate space to save transposed data
    double* transposed_row = (double*) malloc(num_processors*quotient*quotient*sizeof(double));
    //Transpose block matrix  
    MPI_Alltoall(row, quotient*quotient, MPI_DOUBLE, transposed_row, quotient*quotient, MPI_DOUBLE, MPI_COMM_WORLD);
    free(row);
    if(rank ==0 ){
        //write indicator line 
        fprintf(fp, "Transposed\n\n");
    }
    //Transpose each block 
    for (int block_index = 0; block_index<num_processors; block_index++){
        transposed_row=transpose_block(transposed_row, block_index, quotient);
    }
    //Save data in txt file 
    print_with_multiprocessor_block(fp, transposed_row, rank, num_processors, quotient, n);
    free(transposed_row);
} 

//This function tests that generated element is correct for i, j, n = 1,1,1 case
int test_element(){
    //we can calculate the element when n =1 
    double answer =1.0/6.0+ 1.0/5.0+1.0/4.0+1.0/3.0+1.0/2.0+1.0+1.0+1.0/2.0+1.0/3.0+1.0/4.0;
    double err = abs(element(1,1,1)-answer);
    if (err> 0.0001){
        return 0; 
    }else{
        return 1; 
    }
}
//This function tests that two txt file is same or not
int compare_two_txt(char* filename1, char* filename2){
    FILE *fp1, *fp2;
    //Open two files 
    fp1 = fopen(filename1, "r");
    fp2 = fopen(filename2, "r");
    //Test they are same or not 
    while(1){
        int c1 = fgetc(fp1);
        int c2 = fgetc(fp2);
        if(c1 == EOF || c2 == EOF)
            break;
        //There is different element
        if(c1 != c2){
            return 0;
        }
    }
    //They are same!
    fclose(fp1);
    fclose(fp2);
    return 1; 
}
//This tests print_with_multiprocessor function that if it prints in order or not
int test_print_with_multiprocessor(int rank, int num_processors){   
    if (rank == 0){
        //Generate answer file
        FILE* fp = fopen("test_answer.txt", "w");
        for (int i=0; i<num_processors; i++){
            for(int j=0; j<num_processors; j++){
                fprintf(fp, "%.3f ", (double) i);
            }fprintf(fp, "\n");
        }
        fclose(fp);   
    }

    //Generate result file 
    FILE* fp_ret = fopen("test_result.txt", "w");
    double* row = (double*)malloc(num_processors*sizeof(double));
    for(int i=0; i<num_processors; i++){
        row[i]=(double) rank;
    }

    print_with_multiprocessor(fp_ret, row, rank, num_processors,  MPI_COMM_WORLD);
    fclose(fp_ret);

    int pass; 
    if (rank == 0 ){
        //Test in root; 
        char* filename1 = "test_answer.txt";
        char* filename2 = "test_result.txt";
        pass =  compare_two_txt(filename1, filename2);
    }
    //Broadcast the test result
    MPI_Bcast(&pass, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return pass; 
}
//This function test make_and_transpose_small function with 3 processors and dimension =2 
//Compare result with serial results 
int test_make_and_transpose_small(int rank, int num_processors){
    if(num_processors!=3){
        printf("test_make_and_transpose_small function requires 3 processors\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    char* filename1 = "Matrix_serial_dim=2.txt";
    char* filename2 = "Matrix_parallel_dim=2.txt";
    FILE* fp  = fopen(filename2, "w");
    make_and_transpose_small(fp, rank, 2);
    fclose(fp);

    int pass; 
    if (rank == 0 ){
        //Test in root; 
        pass =  compare_two_txt(filename1, filename2);
    }
    //Broadcast the test result
    MPI_Bcast(&pass, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return pass; 
}
//This function tests transpose operation of matrix with simple test matrix 
int test_transpose_block(){
    int quotient = 2 ;
    int block_index = 0 ; 
    //Transpose with simple matrix
    /* | 1, 2 |           | 1, 3 |
       | 3, 4 | ->        | 2, 4 | */
    double* row =(double*)malloc(4*sizeof(double));
    row[0]=1.0; row[1]=2.0; row[2]=3.0; row[3]=4.0;

    double* answer =(double*)malloc(4*sizeof(double));
    answer[0]=1.0; answer[1]=3.0; answer[2]=2.0; answer[3]=4.0; 
    row = transpose_block(row, block_index, quotient);

    int i=0;
    while(1){
        //If there exist different element, return 0 to indicate false(failed)
        if(row[i]!=answer[i]){
            return 0;
        }
        i++;
        if(i>=4){
            break;
        }
    }
    //All elements are same return 1 to indicate true(pass)
    return 1 ; 

}
//This function test make_and_transpose_small function with 3 processors and dimension =5 
//Compare result with serial results 
int test_make_and_transpose_big(int rank, int num_processors){
    if(num_processors!=3){
        printf("test_make_and_transpose_small function requires 3 processors\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    int n= 5; //dimension = 5 
    char* filename1 = "Matrix_serial_dim=5.txt";
    char* filename2 = "Matrix_parallel_dim=5.txt";
    FILE* fp  = fopen(filename2, "w");
    make_and_transpose_big(fp, n, num_processors, rank);
    fclose(fp);

    int pass; 
    if (rank == 0 ){
        //Test in root; 
        pass =  compare_two_txt(filename1, filename2);
    }
    //Broadcast the test result
    MPI_Bcast(&pass, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return pass; 
}






int main(int argc, char* argv[]){
    //Initialize MPI
    MPI_Init(&argc, &argv);

    //Get rank(id) of current processor and number of total processors.
    int rank, num_processors;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processors);

    //Get matrix dimension n 
    int n;
    if(rank == 0 ){
        printf("Enter matrix dimension: ");
        scanf("%d", &n);
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    char filename[100];   
    sprintf(filename, "Matrix_parallel_dim=%d.txt", n);
 
    FILE* fp = fopen(filename, "w");
    
    if(n<=num_processors){
        make_and_transpose_small(fp, rank, n);
        
    }else{
        make_and_transpose_big(fp, n, num_processors, rank);
    }
    fclose(fp);
    MPI_Finalize();
}
// int main(int argc, char* argv[]){
//     //Initialize MPI
//     MPI_Init(&argc, &argv);

//     //Get rank(id) of current processor and number of total processors.
//     int rank, num_processors;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &num_processors);

//     int dims[8] = {2,4,8,16,32,64,128,256};
//     int n;
//     double start, end ; 
//     for(int i=0; i<8; i++){
//         n=dims[i];
//         start = MPI_Wtime();
//         MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

//         char filename[100];   
//         sprintf(filename, "Matrix_parallel_dim=%d.txt", n);
    
//         FILE* fp = fopen(filename, "w");
        
//         if(n<=num_processors){
//             make_and_transpose_small(fp, rank, n);
            
//         }else{
//             make_and_transpose_big(fp, n, num_processors, rank);
//         }
//         fclose(fp);
//         end = MPI_Wtime(); 
//         if (rank == 0 )
//             printf("DIM =%d, time spent %f\n", n, end - start);
//     }
//     MPI_Finalize();
// }
// int main(int argc, char* argv[]){
//     //Initialize MPI
//     MPI_Init(&argc, &argv);

//     //Get rank(id) of current processor and number of total processors.
//     int rank, num_processors;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &num_processors);

//     if(rank == 0 ){
//         //Test for generating function
//         int pass = test_element(); 
//         if(pass == 0 ){
//             printf("!!!element function failed the test!!!\n");
//         }else{
//             printf("element function passed the test\n");
//         }
//         pass = test_transpose_block();
//         if(pass == 0 ){
//             printf("!!!transpose_block function failed the test!!!\n");
//         }else{
//             printf("transpose_block function passed the test\n");
//         }
//     }

//     int pass = test_print_with_multiprocessor(rank, num_processors);    
//     if (rank == 0 ){
//         if(pass == 0 ){
//             printf("!!!print_with_multiprocessor function failed the test!!!\n");
//         }else{
//             printf("print_with_multiprocessor function passed the test\n");
//         }
//     }
//     pass = test_make_and_transpose_small(rank, num_processors);    
//     if (rank == 0 ){
//         if(pass == 0 ){
//             printf("!!!make_and_transpose_small function failed the test!!!\n");
//         }else{
//             printf("make_and_transpose_small function passed the test\n");
//         }
//     }
//     pass = test_make_and_transpose_big(rank, num_processors);    
//     if (rank == 0 ){
//         if(pass == 0 ){
//             printf("!!!make_and_transpose_big function failed the test!!!\n");
//         }else{
//             printf("make_and_transpose_big function passed the test\n");
//         }
//     }
//     MPI_Finalize();
// }