

// this header file will include all the functions and classes 
// that I will need to implement my codes. 


#include<iostream> 
#include<fstream>
#include <stdio.h>

#include<String> 
#include<math.h> 
#include<cstdlib>
#include<ctime> 
#include<assert.h>
#include<vector> 

#ifndef MISCELL
#define MISCELL

using namespace std;  


//void sayHello() ;  //declartion 
int partition(int *a,int *index, int left, int right) ;  // declartion 
void quicksort(int *a, int *index, int left, int right) ;  // declartion 
void printMatrix(int** anMatrix, int ro, int co, string MatrixName) ; 
void printArray(int anArray[], int left, int right, string ArrayName) ; 
int facotrial(int n);
bool reverseWords (char str[]); 
//void charindexarray( char str [], vector<int>  D[] );

#include "miscell.cpp"
#endif 
