#include <iostream>

#include "matrix.hpp"


int main()
{
    matrix A, B, C;


    A = {{1.0, 2,  0},
         {4, 3, -1}};

    B = {{5, 1},
         {2, 3},
         {3, 4}};

    std::cout << A << "\n" << std::endl;
    std::cout << B << "\n" << std::endl;


    C = A*B;

    std::cout << C << "\n" << std::endl;


    matrix X, Y, Z;

    X = {{1, 2, 3},
         {3, 2, 1}};

    Y = {{4, 1, 1},
         {4, 1, 1}};

    std::cout << X << "\n" << std::endl;
    std::cout << Y << "\n" << std::endl;


    Z = X+Y;

    std::cout << Z << "\n" << std::endl;

    Z = Z-Y-X;

    std::cout << Z << "\n" << std::endl;


    A = {{4, 7},
         {2, 6}};

    B = matrix::inv(A);

    std::cout << B << "\n" << std::endl;


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

/*
    #include <iomanip>


    std::cout << std::fixed << std::setprecision(4);
*/

