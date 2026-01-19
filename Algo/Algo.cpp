#include "Algo.h"
using namespace std;

template<class T>
coord<T>::coord(T x, T y):x0(x), y0(y) {
}


void sum(int **A, int **B, coord<int> R_A, coord<int> R_B, int **C, coord<int> R_C, coord<int> sz) {
    for (int i = R_A.x0; i < R_A.x0 + sz.x0; i++) {
        for (int j = R_A.y0; j< R_A.y0 + sz.y0; j++) {
            C[R_C.x0+i-R_A.x0][R_C.y0+j-R_A.y0] = A[i][j] + B[R_B.x0+i-R_A.x0][R_B.y0+j-R_A.y0];
        }
    }
}

void sub(int **A, int **B, coord<int> R_A, coord<int> R_B, int **C, coord<int> R_C, coord<int> sz) {
    for (int i = R_A.x0; i < R_A.x0 + sz.x0; i++) {
        for (int j = R_A.y0; j< R_A.y0 + sz.y0; j++) {
            C[R_C.x0+i-R_A.x0][R_C.y0+j-R_A.y0] = A[i][j] - B[R_B.x0+i-R_A.x0][R_B.y0+j-R_A.y0];
        }
    }
}

void matrixMultiplication_Square(int ** A, int ** B, int N, int **C) {
     
     for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i][j] = 0;
            for (int k = 0; k < N; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
     }
     
}

void printMatrix(int ** A, int N) {
    for (int i = 0; i < N; i++) {
        cout<<((i==0)?"":"\n");
        cout<<A[i][0];
        for (int j = 1; j < N; j++) {
            cout<<" "<<A[i][j];
        }
     }
}

void matrixMultiplication_Strassen(int ** A, int ** B, int N, coord<int> stA, coord<int> stB, int ** C, coord<int> stC, int ***S, int ***P) {
    
    if (N == 1) {
        
        C[stC.x0][stC.y0] = A[stA.x0][stA.y0] * B[stB.x0][stB.y0];
        
        return;
    }
    coord<int> size(N/2, N/2);
    coord<int> start_00(0, 0);
    coord<int> start_11(stA.x0, stA.y0);
    coord<int> start_12(stA.x0, stA.y0 + N/2);
    coord<int> start_21(stA.x0 + N/2, stA.y0);
    coord<int> start_22(stA.x0 + N/2, stA.y0 + N/2);

    coord<int> start_11B(stB.x0, stB.y0);
    coord<int> start_12B(stB.x0, stB.y0 + N/2);
    coord<int> start_21B(stB.x0 + N/2, stB.y0);
    coord<int> start_22B(stB.x0 + N/2, stB.y0 + N/2);

    coord<int> start_11C(stC.x0, stC.y0);
    coord<int> start_12C(stC.x0, stC.y0 + N/2);
    coord<int> start_21C(stC.x0 + N/2, stC.y0);
    coord<int> start_22C(stC.x0 + N/2, stC.y0 + N/2);
    int *** S_new = new int**[11];
    for (int i = 1; i < 11; i++)
    {
        S_new[i] = new int*[N/2];
        for (int j = 0; j < N/2; j++)
        {
            S_new[i][j] = new int[N/2];
        }
        
    }
    
    int *** P_new = new int**[8];
    for (int i = 1; i < 8; i++)
    {
        P_new[i] = new int*[N/2];
        for (int j = 0; j < N/2; j++)
        {
            P_new[i][j] = new int[N/2];
        }
        
    }
    
    // S1 B_12 - B_22
    sub(B, B, start_12B, start_22B, S[1], start_00, size);
    // S2 A_11 + A_12
    sum(A, A, start_11, start_12, S[2], start_00, size);
    // S3 A_21 + A_22
    sum(A, A, start_21, start_22, S[3], start_00, size);
    // S4 B_21 - B_11
    sub(B, B, start_21B, start_11B, S[4], start_00, size);
    // S5 A_11 + A_22
    sum(A, A, start_11, start_22, S[5], start_00, size);
    // S6 B_11 + B_22
    sum(B, B, start_11B, start_22B, S[6], start_00, size);
    // S7 A_12 - A_22
    sub(A, A, start_12, start_22, S[7], start_00, size);
    // S8 B_21 + B_22
    sum(B, B, start_21B, start_22B, S[8], start_00, size);
    // S9 A_11 - A_21
    sub(A, A, start_11, start_21, S[9], start_00, size);
    // S10 B_11 + B_12
    sum(B, B, start_11B, start_12B, S[10], start_00, size);
    
    // P1 A_11 x S_1
    matrixMultiplication_Strassen(A, S[1], N/2, start_11, start_00, P[1], start_00,S_new, P_new);
    
    // P2 S_2 x B_22
    matrixMultiplication_Strassen(S[2], B, N/2, start_00, start_22B, P[2], start_00,S_new, P_new);

    // P3 S_3 x B_11
    matrixMultiplication_Strassen(S[3], B, N/2, start_00, start_11B, P[3], start_00,S_new, P_new);
    
    // P4 A_22 x S_4
    matrixMultiplication_Strassen(A, S[4], N/2, start_22, start_00, P[4], start_00, S_new, P_new);
    
    // P5 S_5 x S_6
    matrixMultiplication_Strassen(S[5], S[6], N/2, start_00, start_00, P[5], start_00, S_new, P_new);
    
    // P6 S_7 x S_8
    matrixMultiplication_Strassen(S[7], S[8], N/2, start_00, start_00, P[6], start_00, S_new, P_new);
    
    // P7 S_9 x S_10
    matrixMultiplication_Strassen(S[9], S[10], N/2, start_00, start_00, P[7], start_00, S_new, P_new);
    
    // C_11 P_5 + P_4 - P_2 + P_6
    sum(P[5], P[4], start_00, start_00, C, start_11C, size);
    sub(C, P[2], start_11C, start_00, C, start_11C, size);
    sum(C, P[6], start_11C, start_00, C, start_11C, size);
    // C_12 P_1 + P_2
    sum(P[1], P[2], start_00, start_00, C, start_12C, size);
    // C_21 P_3 + P_4
    sum(P[3], P[4], start_00, start_00, C, start_21C, size);
    // C_22 P_5 + P_1 - P_3 - P_7
    sum(P[5], P[1], start_00, start_00, C, start_22C, size);
    sub(C, P[3], start_22C, start_00, C, start_22C, size);
    sub(C, P[7], start_22C, start_00, C, start_22C, size);
}

void matrixMultiplication_SA(int ** A, int ** B, int N, int **C) {
    coord<int> st(0, 0);
    int *** S_new = new int**[11];
    for (int i = 0; i < 11; i++)
    {
        S_new[i] = new int*[N];
        for (int j = 0; j < N; j++)
        {
            S_new[i][j] = new int[N];
        }
        
    }
    
    int *** P_new = new int**[8];
    for (int i = 0; i < 8; i++)
    {
        P_new[i] = new int*[N];
        for (int j = 0; j < N; j++)
        {
            P_new[i][j] = new int[N];
        }
        
    }
    matrixMultiplication_Strassen(A, B, N, st, st, C, st, S_new, P_new);
    
}