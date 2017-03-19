#include <iostream>
#include <cassert>
#include <cmath>

#include "matrix.hpp"


int main()
{
    matrix A, B, C, X, Y;


    // TEST: MTIMES
    // ----------------------------------------
    A = {{1.0, 2,  0},
         {4, 3.0, -1.0}};

    B = {{5, 1},
         {2.0, 3},
         {3, 4}};

    X = {{9, 7.0},
         {23, 9}};

    C = A*B;

    assert(C == X);

    std::cout << C << "\n" << std::endl;
    // ----------------------------------------




    // TEST: INV (A square matrix N-by-N)
    // ----------------------------------------
    A = {{4, 7},
         {2, 6}};

    Y = matrix::inv(A);

    X = {{0.6, -0.7},
         {-0.2, 0.4}};

    assert(Y == X);



    std::cout << Y << "\n" << std::endl;
    // ----------------------------------------




    // TEST: MLDIVIDE (A square matrix N-by-N)
    // ----------------------------------------
    /*
        EXEMPLE 1:
        Solve the following linear system of equations (Ax = B):
            x + y = 2
            x - y = 0
    */

    A = {{1,  1},
         {1, -1}};

    B = {{2},
         {0}};

    Y = matrix::mldivide(A, B);

    X = {{1},
         {1}};

    assert(Y == X);
    assert(A*X == B);


    /*
        EXEMPLE 2:
    */
    A = {{ 6.80, -6.05, -0.45,  8.32, -9.67},
         {-2.11, -3.30,  2.58,  2.71, -5.14},
         { 5.66,  5.36, -2.70,  4.35, -7.26},
         { 5.97, -4.44,  0.27, -7.17,  6.08},
         { 8.23,  1.08,  9.04,  2.14, -6.87}};

    B = {{ 4.02, -1.56,  9.81},
         { 6.19,  4.00, -4.09},
         {-8.22, -8.67, -4.57},
         {-7.57,  1.75, -8.61},
         {-3.03,  2.86,  8.99}};

    X = matrix::mldivide(A, B);

    assert(A*X == B);

    std::cout << X << "\n" << std::endl;
    // ----------------------------------------


    return 0;
}


/*
    #include <chrono>


    std::chrono::time_point<std::chrono::system_clock> start, end;

    start = std::chrono::system_clock::now();
    Z = X*Y;
    end = std::chrono::system_clock::now();

    int64_t elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << elapsed_time << std::endl;
*/

