#ifndef ALGO_H
#define ALGO_H
#include <iostream>

template<class T>
struct coord{
       T x0;
       T y0;
       coord(T x, T y);
       
};

// Print matrix
void printMatrix(int ** A, int N);

// Sum matrix 
void sum(int **A, int **B, coord<int> R_A, coord<int> R_B, int **C, coord<int> R_C, coord<int> sz);

// Sub matrix 
void sub(int **A, int **B, coord<int> R_A, coord<int> R_B, int **C, coord<int> R_C, coord<int> sz);

// O(N^3) Algorithms 
void matrixMultiplication_Square(int ** A, int ** B, int N, int **C);

// O(N^{log2(7)}) Algorithms - Strassen's Algorithms
void matrixMultiplication_Strassen(int ** A, int ** B, int N, coord<int> stA, coord<int> stB, int ** C, coord<int> stC, int ***S, int ***P);

// User Interface
void matrixMultiplication_SA(int ** A, int ** B, int N, int **C);

#endif