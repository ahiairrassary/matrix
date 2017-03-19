#include <iostream>
#include <cassert>
#include <cmath>

#include "matrix.hpp"


int main()
{
    matrix A, B, C, D, X, Y, Z;


    A = {{1.0, 2,  0},
         {4, 3.0, -1.0}};

    B = {{5, 1},
         {2.0, 3},
         {3, 4}};

    C = A*B;

    D = {{9, 7.0},
         {23, 9}};

    assert(C == D);

    std::cout << C << "\n" << std::endl;


    A = {{4, 7},
         {2, 6}};

    C = matrix::inv(A);

    D = {{0.6, -0.7},
         {-0.2, 0.4}};

    assert(C == D);

    std::cout << C << "\n" << std::endl;

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

