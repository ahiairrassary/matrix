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
    4 | *  *  *  *  * |
      --             --

    Size: 4x5 --> 4 rows and 5 columns

    /!\ Elements are stored with column-major ordering (ie. column by column) !
*/

public:
    matrix();
    matrix(size_t number);
    matrix(size_t r_number, size_t c_number);
    matrix(const matrix & X);
    matrix(const std::initializer_list<std::initializer_list<double>> & list);

    matrix & operator=(const matrix & X);

    double & at(size_t r_idx, size_t c_idx);
    const double & at(size_t r_idx, size_t c_idx) const;

    double & operator()(size_t r_idx, size_t c_idx);
    const double & operator()(size_t r_idx, size_t c_idx) const;

    size_t column() const;
    size_t row() const;

    void print(std::ostream & os) const;

    static matrix plus(const matrix & A, const matrix & B);
    static matrix minus(const matrix & A, const matrix & B);
    static matrix mtimes(const matrix & A, const matrix & B);
    static matrix times(const matrix & A, const matrix & B);

    static matrix inv(const matrix & X);


private:
    std::vector<double> m_data;

    size_t m_column;
    size_t m_row;


};


matrix operator+(const matrix & A, const matrix & B);

matrix operator-(const matrix & A, const matrix & B);

matrix operator*(const matrix & A, const matrix & B);

std::ostream & operator<<(std::ostream & os, const matrix & X);


#endif // MATRIX_HPP

