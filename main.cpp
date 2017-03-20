#include <iostream>
#include <cassert>
#include <cmath>

#include "matrix.hpp"


int main()
{
    matrix A, B, C, X;


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




    // TEST: INV (A square matrix N-by-N)
    // ----------------------------------------
    A = {{4, 7},
         {2, 6}};

    X = matrix::inv(A);

    B = {{ 0.6, -0.7},
         {-0.2,  0.4}};

    std::cout << X << "\n" << std::endl;

    assert(X == B);
    assert(A*X == X*A);
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
        EXEMPLE 1:
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
        EXEMPLE 2:
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
        EXEMPLE 1:
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
        EXEMPLE 2:
    */
    A = {{0.5,  0.5, 0.5,  0.5},
         {0.5, -1.5, 0.5, -1.5},
         {1.0,  1.0, 0.0,  1.0}};
    
    B = {{1.0,  1.0, 1.0, 0.0},
         {1.0, -1.0, 2.5, 1.0},
         {1.0,  1.0, 3.0, 0.0}};

    X = matrix::mldivide(A, B);

    C = {{1.0, 0.0,   3.5,   0.5},
         {0.0, 0.5, -0.25, -0.25},
         {1.0, 1.0,  -1.0,   0.0},
         {0.0, 0.5, -0.25, -0.25}};

std::cout << A << "\n" << std::endl; // SUPPRIMER !!!!
std::cout << B << "\n" << std::endl; // SUPPRIMER !!!!
    std::cout << X << "\n" << std::endl;
std::cout << A*X << "\n" << std::endl; // SUPPRIMER !!!!

//    assert(X == C);
//    assert(A*X == B);
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

