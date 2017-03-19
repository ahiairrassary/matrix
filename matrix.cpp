#include "matrix.hpp"

#include <cmath>
#include <iomanip>

//#include <cblas.h> /* LINUX */
#include <Accelerate/Accelerate.h> /* MACOS */


matrix::matrix()
{
    m_columns = 1;
    m_rows = 1;

    m_data.resize(1, 0.0);
}


matrix::matrix(size_t n)
{
    if(n == 0)
    {
       throw std::runtime_error("matrix: constructor: dimension error");
    }

    m_columns = n;
    m_rows = n;

    m_data.resize(columns()*rows(), 0.0);
}


matrix::matrix(size_t n, size_t m)
{
    if(n == 0 || m == 0)
    {
       throw std::runtime_error("matrix: constructor: dimension error");
    }

    m_columns = m;
    m_rows = n;

    m_data.resize(columns()*rows(), 0.0);
}


matrix::matrix(const matrix & X)
{
    m_columns = X.m_columns;
    m_rows = X.m_rows;

    m_data = X.m_data;
}


matrix::matrix(const std::initializer_list<std::initializer_list<double>> & list)
{
    /* TODO : Tester que toutes les listes font la bonne taille ! */

    m_columns = (*list.begin()).size();
    m_rows = list.size();

    m_data.reserve(columns()*rows());

    /* Column-major order */
    for(size_t i(0); i < columns(); ++i)
    {
        for(const std::initializer_list<double> & l : list)
        {
            m_data.push_back(*(l.begin() + i));
        }
    }

    /* Row-major order */
    /*
    for(const std::initializer_list<double> & l : list)
    {
        m_data.insert(m_data.end(), l.begin(), l.end());
    }
    */
}


matrix & matrix::operator=(const matrix & X)
{
    m_columns = X.m_columns;
    m_rows = X.m_rows;

    m_data = X.m_data;

    return *this;
}


double & matrix::at(size_t row_idx, size_t column_idx)
{
    if(row_idx >= rows() || column_idx >= columns())
    {
        throw std::runtime_error("matrix: at: out of bounds");
    }

    return m_data.at(row_idx + column_idx*rows());
}


const double & matrix::at(size_t row_idx, size_t column_idx) const
{
    if(row_idx >= rows() || column_idx >= columns())
    {
        throw std::runtime_error("matrix: at: out of bounds");
    }

    return m_data.at(row_idx + column_idx*rows());
}


double & matrix::operator()(size_t row_idx, size_t column_idx)
{
    return at(row_idx, column_idx);
}


const double & matrix::operator()(size_t row_idx, size_t column_idx) const
{
    return at(row_idx, column_idx);
}


size_t matrix::columns() const
{
    return m_columns;
}


size_t matrix::rows() const
{
    return m_rows;
}


void matrix::print(std::ostream & os) const
{
    os << std::fixed << std::setprecision(3) << std::right << std::showpos;

    for(size_t r(0); r < rows(); ++r)
    {
        for(size_t c(0); c < columns(); ++c)
        {
            os << std::setw(9) << at(r, c) << ' ';
        }
        
        os << '\n';
    }

    os << "(" << rows() << "x" << columns() << ")";
}


/*
    MATRIX::PLUS
    C = A + B adds arrays A and B and returns the result in C.
*/
matrix matrix::plus(const matrix & A, const matrix & B)
{
    if(A.rows() != B.rows() || A.columns() != B.columns())
    {
        throw std::runtime_error("matrix: plus: size error");
    }

    matrix C(A.rows(), A.columns());

    for(size_t i(0); i < C.m_data.size(); ++i)
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
    if(A.rows() != B.rows() || A.columns() != B.columns())
    {
        throw std::runtime_error("matrix: minus: size error");
    }

    matrix C(A.rows(), A.columns());

    for(size_t i(0); i < C.m_data.size(); ++i)
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
    if(A.columns() != B.rows())
    {
        throw std::runtime_error("matrix: mtimes: size error");
    }

    matrix C(A.rows(), B.columns());

    /* DGEMM  performs one of the matrix-matrix operations
           C := alpha*op( A )*op( B ) + beta*C
    */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                (int)A.rows(), (int)B.columns(), (int)A.columns(),
                1.0, /* alpha */
                &A.m_data[0], (int)A.rows(),
                &B.m_data[0], (int)B.rows(),
                0.0, /* beta */
                &C.m_data[0], (int)C.rows());

    return C;
}


/*
    MATRIX::TIMES
    C = times(A,B) multiplies arrays A and B element by element and returns the result in C.
*/
matrix matrix::times(const matrix & A, const matrix & B)
{
    if(A.rows() != B.rows() || A.columns() != B.columns())
    {
        throw std::runtime_error("matrix: times: size error");
    }

    matrix C(A.rows(), A.columns());

    for(size_t i(0); i < C.m_data.size(); ++i)
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
    /* TODO : Verifier les infos retournées par la variables info ! */

    if(X.rows() != X.columns())
    {
        throw std::runtime_error("matrix: inv: size error");
    }

    matrix Y = X;

    std::vector<int> ipiv(std::min(Y.rows(), Y.columns()));
    int m = (int)Y.rows();
    int n = (int)Y.columns();
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


bool matrix::equal(const matrix & A, const matrix & B, const double epsilon)
{
    bool result = (A.rows() == B.rows() && A.columns() == B.columns());


    for(size_t i(0); i < A.m_data.size() && result; ++i)
    {
        if(!(std::abs(A.m_data[i] - B.m_data[i]) < epsilon))
        {
            result = false;
        }
    }

    return result;
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


bool operator==(const matrix & A, const matrix & B)
{
    return matrix::equal(A, B, 1.0e-12);
}


bool operator!=(const matrix & A, const matrix & B)
{
    return !matrix::equal(A, B, 1.0e-12);
}


std::ostream & operator<<(std::ostream & os, const matrix & X)
{
    X.print(os);

    return os;
}

