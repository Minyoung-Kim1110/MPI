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



int main(int argc, char *argv[]){