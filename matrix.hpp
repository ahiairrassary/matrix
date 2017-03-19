#ifndef MATRIX_HPP
#define MATRIX_HPP


#include <initializer_list>
#include <iostream>
#include <vector>


class matrix
{

/*
    MATRIX

        1  2  3  4  5
      --             --
    1 | *  *  *  *  * |
    2 | *  *  *  *  * |
    3 | *  *  *  *  * |
      --             --

    Size: 3x5 --> 3 rows and 5 columns

    /!\ Elements are stored with column-major ordering (ie. column by column) !
*/

public:
    matrix();
    matrix(size_t n);
    matrix(size_t n, size_t m);
    matrix(const matrix & X);
    matrix(const std::initializer_list<std::initializer_list<double>> & list);

    matrix & operator=(const matrix & X);

    double & at(size_t row_idx, size_t column_idx);
    const double & at(size_t row_idx, size_t column_idx) const;

    double & operator()(size_t row_idx, size_t column_idx);
    const double & operator()(size_t row_idx, size_t column_idx) const;

    size_t columns() const;
    size_t rows() const;

    void print(std::ostream & os) const;

    static matrix plus(const matrix & A, const matrix & B);
    static matrix minus(const matrix & A, const matrix & B);
    static matrix mldivide(const matrix & A, const matrix & B);
    static matrix mtimes(const matrix & A, const matrix & B);
    static matrix times(const matrix & A, const matrix & B);

    static matrix inv(const matrix & X);

    static bool equal(const matrix & A, const matrix & B, const double epsilon);


private:
    std::vector<double> m_data;

    size_t m_columns;
    size_t m_rows;


};


matrix operator+(const matrix & A, const matrix & B);

matrix operator-(const matrix & A, const matrix & B);

matrix operator*(const matrix & A, const matrix & B);

bool operator==(const matrix & A, const matrix & B);

bool operator!=(const matrix & A, const matrix & B);

std::ostream & operator<<(std::ostream & os, const matrix & X);


#endif // MATRIX_HPP

