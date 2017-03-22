#include <iostream>
#include <cassert>
#include <cmath>

#include "matrix.hpp"


int main()
{
    matrix A, B, C, X;


    // TEST: PLUS, MINUS & TIMES
    // ----------------------------------------
    A = {{1.0, 2.0, 3.0, 4.0},
         {5.0, 6.0, 7.0, 8.0}};

    X = 20.0 + A;

    B = {{21.0, 22.0, 23.0, 24.0},
         {25.0, 26.0, 27.0, 28.0}};

    std::cout << X << "\n" << std::endl;

    assert(X == B);


    X = 20.0 - B;

    B = {{-1.0, -2.0, -3.0, -4.0},
         {-5.0, -6.0, -7.0, -8.0}};

    std::cout << X << "\n" << std::endl;

    assert(X == B);


    X =  2.5 * A;

    B = {{ 2.5,  5.0,  7.5, 10.0},
         {12.5, 15.0, 17.5, 20.0}};

    std::cout << X << "\n" << std::endl;

    assert(X == B);
    // ----------------------------------------




    // TEST: TRANSPOSE
    // ----------------------------------------
    A = {{1.0, 2.0, 3.0, 4.0},
         {5.0, 6.0, 7.0, 8.0}};

    X = matrix::transpose(A);

    B = {{1, 5},
         {2, 6},
         {3, 7},
         {4, 8}};

    std::cout << X << "\n" << std::endl;

    assert(X == B);
    assert(A == matrix::transpose(X));
    // ----------------------------------------




    // TEST: MTIMES
    // ----------------------------------------
    A = {{1.0,   2,    0},
         {  4, 3.0, -1.0}};

    B = {{  5, 1},
         {2.0, 3},
         {  3, 4}};

    X = A*B;

    C = {{ 9, 7.0},
         {23,   9}};

    std::cout << X << "\n" << std::endl;

    assert(X == C);
    // ----------------------------------------




    // TEST: INV (Square matrix)
    // ----------------------------------------
    /*
        EXEMPLE 1:
    */
    A = {{4, 7},
         {2, 6}};

    X = matrix::inv(A);

    B = {{ 0.6, -0.7},
         {-0.2,  0.4}};

    std::cout << X << "\n" << std::endl;

    assert(X == B);
    assert(A*X == X*A);
    assert(A == matrix::inv(X));


    /*
        EXEMPLE 2:
    */
    A = {{2.0, 3.0},
         {4.0, 5.0}};

    X = matrix::inv(A);

    B = {{-2.5,  1.5},
         { 2.0, -1.0}};

    std::cout << X << "\n" << std::endl;

    assert(X == B);
    assert(A*X == X*A);
    assert(A == matrix::inv(X));
    // ----------------------------------------




    // TEST: MLDIVIDE (Square matrix)
    // ----------------------------------------
    /*
        EXEMPLE 1: (m = n)
        Solve the following linear system of equations (Ax = B):
            x + y = 2
            x - y = 0
    */

    A = {{1,  1},
         {1, -1}};

    B = {{2},
         {0}};

    X = matrix::mldivide(A, B);

    C = {{1},
         {1}};

    assert(X == C);
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

    std::cout << X << "\n" << std::endl;

    assert(A*X == B);
    // ----------------------------------------




    // TEST: MLDIVIDE (Least-Squares Solution of Overdetermined System)
    // ----------------------------------------
    /*
        EXEMPLE 1: (m > n)
    */
    A = {{1, 1, 1},
         {2, 3, 4},
         {3, 5, 2},
         {4, 2, 5},
         {5, 4, 3}};

    B = {{-10, -3},
         { 12, 14},
         { 14, 12},
         { 16, 16},
         { 18, 16}};

    X = matrix::mldivide(A, B);

    C = {{2, 1},
         {1, 1},
         {1, 2}};

    std::cout << X << "\n" << std::endl;

    assert(X == C);
    //assert(A*X == B);


    /*
        EXEMPLE 2: (m > n)
    */
    A = {{0.0,  2.0},
         {2.0, -1.0},
         {2.0, -1.0},
         {0.0,  1.5},
         {2.0, -1.0},
         {2.0, -1.0}};

    B = {{1.0,  4.0, 1.0},
         {1.0,  1.0, 2.0},
         {1.0, -1.0, 1.0},
         {1.0,  3.0, 2.0},
         {1.0,  1.0, 1.0},
         {1.0, -1.0, 1.0}};

    X = matrix::mldivide(A, B);

    C = {{0.78, 1.0, 1.025},
         {0.56, 2.0, 0.800}};

    std::cout << X << "\n" << std::endl;

    assert(X == C);
    //assert(A*X == B);
    // ----------------------------------------




    // TEST: MLDIVIDE (Least-Squares Solution of Underdetermined System)
    // ----------------------------------------
    /*
        EXEMPLE 1: (m < n)
    */
    A = {{0.0,  2.0,  2.0, 0.0,  2.0,  2.0},
         {2.0, -1.0, -1.0, 1.5, -1.0, -1.0}};
    
    B = {{1.0},
         {1.0}};

    X = matrix::mldivide(A, B);

    C = {{0.480},
         {0.125},
         {0.125},
         {0.360},
         {0.125},
         {0.125}};

    std::cout << X << "\n" << std::endl;

    assert(X == C);
    assert(A*X == B);


    /*
        EXEMPLE 2: (m < n)
    */
    A = {{1, 1, 1},
         {1, 1, 2}};

    B = {{1},
         {3}};

    X = matrix::mldivide(A, B);

    C = {{-0.5},
         {-0.5},
         { 2.0}};

    std::cout << X << "\n" << std::endl;

    assert(X == C);
    assert(A*X == B);


    /*
        EXEMPLE 3: (m < n)
    */

    A = {{0.0,  2.0},
         {2.0, -1.0},
         {2.0, -1.0},
         {0.0,  1.5},
         {2.0, -1.0},
         {2.0, -1.0}};

    B = {{1.0,  4.0,  1.0},
         {1.0,  1.0,  2.0},
         {1.0, -1.0,  1.0},
         {1.0,  3.0,  2.0},
         {1.0,  1.0,  1.0},
         {1.0, -1.0,  1.0}};

    X = matrix::mldivide(A, B);

    C = {{0.7800, 1.0000, 1.0250},
         {0.5600, 2.0000, 0.8000}};

    std::cout << X << "\n" << std::endl;

    assert(X == C);
    // ----------------------------------------


    // TEST: PINV
    // ----------------------------------------
    /*
        EXEMPLE 1: (m > n)
    */
    A = {{64,  2,  3, 61, 60,  6},
         { 9, 55, 54, 12, 13, 51},
         {17, 47, 46, 20, 21, 43},
         {40, 26, 27, 37, 36, 30},
         {32, 34, 35, 29, 28, 38},
         {41, 23, 22, 44, 45, 19},
         {49, 15, 14, 52, 53, 11},
         { 8, 58, 59,  5,  4, 62}};

    X = matrix::pinv(A);

    std::cout << X << "\n" << std::endl;

    assert(A == matrix::pinv(X));
    assert(A*X*A == A);
    assert(X*A*X == X);


    /*
        EXEMPLE 2: (m < n)
    */
    A = {{64,  2,  3, 61, 60,  6},
         { 9, 55, 54, 12, 13, 51}};

    X = matrix::pinv(A);

    std::cout << X << "\n" << std::endl;

    assert(A == matrix::pinv(X));
    assert(A*X*A == A);
    assert(X*A*X == X);


    /*
        EXEMPLE 3: (m = n)
    */
    A = {{ 12, -27,  36,  44,  59},
         { 26,  41, -54,  24,  23},
         { 33,  70,  59,  15, -68},
         { 43,  16,  29, -52, -61},
         {-43,  20,  71,  88,  11}};

    X = matrix::pinv(A);

    std::cout << X << "\n" << std::endl;

    assert(A == matrix::pinv(X));
    assert(A*X*A == A);
    assert(X*A*X == X);
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

