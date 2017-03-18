#include "matrix.hpp"

#include <cmath>

//#include <cblas.h> /* LINUX */
#include <Accelerate/Accelerate.h> /* MACOS */


matrix::matrix()
{
    m_column = 1;
    m_row = 1;

    m_data.resize(1, 0.0);
}


matrix::matrix(size_t number)
{
    if(number == 0)
    {
       throw std::runtime_error("matrix: constructor: size error");
    }

    m_column = number;
    m_row = number;

    m_data.resize(column()*row(), 0.0);
}


matrix::matrix(size_t r_number, size_t c_number)
{
    if(r_number == 0 || c_number == 0)
    {
       throw std::runtime_error("matrix: constructor: size error");
    }

    m_column = c_number;
    m_row = r_number;

    m_data.resize(column()*row(), 0.0);
}


matrix::matrix(const matrix & X)
{
    m_column = X.m_column;
    m_row = X.m_row;

    m_data = X.m_data;
}


matrix::matrix(const std::initializer_list<std::initializer_list<double>> & list)
{
    /* TODO : Tester que toutes les listes font la bonne taille ! */

    m_column = (*list.begin()).size();
    m_row = list.size();

    /* Column-major order */
    m_data.reserve(column()*row());

    for(size_t c(0); c < column(); ++c)
    {
        for(const std::initializer_list<double> & l : list)
        {
            m_data.push_back(*(l.begin() + c));
        }
    }

    /* Row-major order */
    /*
    m_data.reserve(column()*row());

    for(const std::initializer_list<double> & l : list)
    {
        m_data.insert(m_data.end(), l.begin(), l.end());
    }
    */
}


matrix & matrix::operator=(const matrix & X)
{
    m_column = X.m_column;
    m_row = X.m_row;

    m_data = X.m_data;

    return *this;
}


double & matrix::at(size_t r_idx, size_t c_idx)
{
    if(r_idx >= row() || c_idx >= column())
    {
        throw std::runtime_error("matrix: at: out of bounds");
    }

    return m_data.at(r_idx + c_idx * row());
}


const double & matrix::at(size_t r_idx, size_t c_idx) const
{
    if(r_idx >= row() || c_idx >= column())
    {
        throw std::runtime_error("matrix: at: out of bounds");
    }

    return m_data.at(r_idx + c_idx * row());
}


double & matrix::operator()(size_t r_idx, size_t c_idx)
{
    return at(r_idx, c_idx);
}


const double & matrix::operator()(size_t r_idx, size_t c_idx) const
{
    return at(r_idx, c_idx);
}


size_t matrix::column() const
{
    return m_column;
}


size_t matrix::row() const
{
    return m_row;
}


void matrix::print(std::ostream & os) const
{
    /* TODO : Ameliorer l'affichage des matrices ! */

    for(size_t r(0); r < row(); ++r)
    {
        for(size_t c(0); c < column(); ++c)
        {
            os << at(r, c) << ' ';
        }
        
        os << '\n';
    }

    os << "(" << row() << "x" << column() << ")";
}


/*
    MATRIX::PLUS
    C = A + B adds arrays A and B and returns the result in C.
*/
matrix matrix::plus(const matrix & A, const matrix & B)
{
    if(A.row() != B.row() || A.column() != B.column())
    {
        throw std::runtime_error("matrix: plus: size error");
    }

    matrix C(A.row(), A.column());

    for(size_t i(0); i < A.m_data.size(); ++i)
    {
        C.m_data[i] = A.m_data[i] + B.m_data[i];
    }

    return C;
}


/*
    MATRIX::MINUS
    C = A - B subtracts array B from array A and returns the result in C.
*/
matrix matrix::minus(const matrix & A, const matrix & B)
{
    if(A.row() != B.row() || A.column() != B.column())
    {
        throw std::runtime_error("matrix: minus: size error");
    }

    matrix C(A.row(), A.column());

    for(size_t i(0); i < A.m_data.size(); ++i)
    {
        C.m_data[i] = A.m_data[i] - B.m_data[i];
    }

    return C;
}


/*
    MATRIX::MTIMES
    C = A*B is the matrix product of A and B.
*/
matrix matrix::mtimes(const matrix & A, const matrix & B)
{
    if(A.column() != B.row())
    {
        throw std::runtime_error("matrix: mtimes: size error");
    }

    matrix C(A.row(), B.column());

    /* DGEMM  performs one of the matrix-matrix operations
           C := alpha*op( A )*op( B ) + beta*C
    */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                (int)A.row(), (int)B.column(), (int)A.column(),
                1.0, /* alpha */
                &A.m_data[0], (int)A.row(),
                &B.m_data[0], (int)B.row(),
                0.0, /* beta */
                &C.m_data[0], (int)C.row());

    return C;
}


/*
    MATRIX::TIMES
    C = times(A,B) multiplies arrays A and B element by element and returns the result in C.
*/
matrix matrix::times(const matrix & A, const matrix & B)
{
    if(A.row() != B.row() || A.column() != B.column())
    {
        throw std::runtime_error("matrix: times: size error");
    }

    matrix C(A.row(), A.column());

    for(size_t i(0); i < A.m_data.size(); ++i)
    {
        C.m_data[i] = A.m_data[i] * B.m_data[i];
    }

    return C;
}


/*
    MATRIX::INV
    Y = inv(X) computes the inverse of square matrix X.
*/
matrix matrix::inv(const matrix & X)
{
    /* TODO : Verifier les infos retournées par la varibales info ! */

    if(X.row() != X.column())
    {
        throw std::runtime_error("matrix: inv: size error");
    }

    matrix Y = X;

    std::vector<int> ipiv(std::min(Y.row(), Y.column()));
    int m = (int)Y.row();
    int n = (int)Y.column();
    int info = 0;

    /* DGETRF computes an LU factorization of a general M-by-N matrix A
       using partial pivoting with row interchanges.
    */
    dgetrf_(&m, &n,
            &Y.m_data[0], &m,
            &ipiv[0],
            &info);

    if(info != 0)
    {
        throw std::runtime_error("matrix: inv: dgetrf error");
    }


    int lwork = m*m;
    std::vector<double> work((size_t)lwork);

    /* DGETRI computes the inverse of a matrix using the LU factorization
       computed by DGETRF.
    */
    dgetri_(&m,
            &Y.m_data[0], &m,
            &ipiv[0],
            &work[0], &lwork,
            &info);

    if(info != 0)
    {
        throw std::runtime_error("matrix: inv: dgetri error");
    }

    return Y;
}




matrix operator+(const matrix & A, const matrix & B)
{
    return matrix::plus(A, B);
}


matrix operator-(const matrix & A, const matrix & B)
{
    return matrix::minus(A, B);
}


matrix operator*(const matrix & A, const matrix & B)
{
    return matrix::mtimes(A, B);
}


std::ostream & operator<<(std::ostream & os, const matrix & X)
{
    X.print(os);

    return os;
}

