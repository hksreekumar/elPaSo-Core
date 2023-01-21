/* Copyright (c) 2023. Authors listed in AUTHORS.md

 * This file is part of elPaSo-Core.

 * elPaSo-Core is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your option)
 * any later version.

 * elPaSo-Core is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.

 * You should have received a copy of the GNU Lesser General Public License along
 * with elPaSo-Core (COPYING.txt and COPYING.LESSER.txt). If not, see
 * <https://www.gnu.org/licenses/>. 
 */

#ifndef INFAM_MATRIX_H
#define INFAM_MATRIX_H

#include <cstring>
#include <iostream>

#include "vektor.h"


namespace infam {

// declaration of global functions before tMatrix class template may be
// necessary for some cases, read:
// https://stackoverflow.com/questions/30379397/template-in-gcc-doesnt-see-a-function-declared-after-it-but-before-instantiati
template <class TFLOAT1, class TFLOAT2>
bool approx_equal(const TFLOAT1 &a, const TFLOAT2 &b, const double &tol);

//! @brief simple matrix class
//! @author Dirk Clasen
//! @date 22.06.2005
//!
//! matrix template matrix. The values are stored in a vector row by row.
//! This is the same memory layout like C uses.
template <class TFLOAT, class TUINT>
class tMatrix {
 private:
  TUINT m_Rows;         //! number of rows
  TUINT m_Cols;         //! number of columns
  TFLOAT *m_Data;       //! matrix entries row by row
  TUINT m_InsertIndex;  //! stores next index to be accessed by '<<'

  //! compute determinant of matrix
  double computeDeterminant(double determinant, double detfaktor) const;

 public:
  tMatrix();
  tMatrix(TUINT Rows, TUINT Cols);
  tMatrix(const tMatrix<TFLOAT, TUINT> &other);
  ~tMatrix();

  //! return the number of rows of the matrix
  inline TUINT rows(void) const { return m_Rows; }

  //! return the number of columns of the matrix
  inline TUINT cols(void) const { return m_Cols; }

  //! return the current insert index
  inline TUINT getInsertIndex(void) const { return m_InsertIndex; }

  //! set new size of matrix
  void resize(TUINT Rows, TUINT Cols);

  //! return a single element of the matrix (read-only)
  inline TFLOAT operator()(TUINT row, TUINT col) const {
    return m_Data[col + row * m_Cols];
  }

  //! access a single element of the matrix (read/write)
  inline TFLOAT &operator()(TUINT row, TUINT col) {
    return m_Data[col + row * m_Cols];
  }

  //! allocation operator " = " for class
  tMatrix<TFLOAT, TUINT> operator=(const tMatrix<TFLOAT, TUINT> &other);

  //! assign value to m_Data[m_InsertIndex], then ++m_InsertIndex,
  //! exception if m_InsertIndex == m_Rows * m_Cols (true for copied matrices!)
  void insertNext(TFLOAT value);

  //! set the same value to all matrix cells
  void setValue(const TFLOAT &value);

  //! check for symmetry
  bool isSymmetric(const double &tol = 1.0e-6) const;

  //! returns determinant of matrix
  double getDeterminant(void) const;

  //! returns true if matrix is invertible
  bool isInvertible(const double &tol = 1.0e-15) const;

  //! returns true for a square matrix that is not invertible
  bool isSingular(const double &tol = 1.0e-15) const;

  //! returns true if the matrix is equal to zero within relative tolerance
  bool isZero(const double &tol = 1.0e-15) const;

  //! returns absolute maximum of matrix
  TFLOAT getAbsMax(void) const;

  //! returns inverted matrix
  tMatrix<TFLOAT, TUINT> inv(void) const;

  //! returns transposed matrix
  tMatrix<TFLOAT, TUINT> trans(void) const;

  //! returns element-wise negative of matrix
  tMatrix<TFLOAT, TUINT> operator-() const;

  //! write this matrix to a stream
  std::ostream &write(std::ostream &os) const;
};

template <class TFLOAT, class TUINT>
tMatrix<TFLOAT, TUINT>::tMatrix() {
  m_Rows = 0;
  m_Cols = 0;
  m_Data = new TFLOAT[m_Rows * m_Cols];
  m_InsertIndex = 0;
}

template <class TFLOAT, class TUINT>
tMatrix<TFLOAT, TUINT>::tMatrix(TUINT Rows, TUINT Cols) {
  m_Rows = Rows;
  m_Cols = Cols;
  m_Data = new TFLOAT[m_Rows * m_Cols];
  for (int k = 0; k < m_Rows * m_Cols; k++) m_Data[k] = 0.;
  m_InsertIndex = 0;
}

template <class TFLOAT, class TUINT>
tMatrix<TFLOAT, TUINT>::tMatrix(const tMatrix<TFLOAT, TUINT> &other) {
  m_Rows = other.rows();
  m_Cols = other.cols();
  m_Data = new TFLOAT[m_Rows * m_Cols];
  m_InsertIndex = m_Rows * m_Cols;

  for (TUINT r = 0; r < rows(); r++)
    for (TUINT c = 0; c < cols(); c++) m_Data[c + r * m_Cols] = other(r, c);
}

template <class TFLOAT, class TUINT>
tMatrix<TFLOAT, TUINT>::~tMatrix() {
  delete[] m_Data;
}

template <class TFLOAT, class TUINT>
inline void tMatrix<TFLOAT, TUINT>::resize(TUINT Rows, TUINT Cols) {
  if (m_Rows != Rows || m_Cols != Cols) {
    m_Rows = Rows;
    m_Cols = Cols;
    delete[] m_Data;
    m_Data = new TFLOAT[m_Rows * m_Cols];
    for (int k = 0; k < m_Rows * m_Cols; k++) m_Data[k] = 0.;
  }
}

template <class TFLOAT, class TUINT>
inline void tMatrix<TFLOAT, TUINT>::setValue(const TFLOAT &value) {
  for (TUINT k = 0; k < m_Rows * m_Cols; k++) m_Data[k] = value;
}

//! allocation operator " = " for class
template <class TFLOAT, class TUINT>
tMatrix<TFLOAT, TUINT> tMatrix<TFLOAT, TUINT>::operator=(
    const tMatrix<TFLOAT, TUINT> &other) {
  if (m_Rows != other.rows() || m_Cols != other.cols()) {
    m_Rows = other.rows();
    m_Cols = other.cols();
    delete[] m_Data;
    m_Data = new TFLOAT[m_Rows * m_Cols];
  }

  for (TUINT r = 0; r < rows(); r++)
    for (TUINT c = 0; c < cols(); c++) m_Data[c + r * m_Cols] = other(r, c);

  return *this;
}

// ----------------------------------------------------------------------------
//   '-'-operator to change sign of matrix elements, returns B = -A
// ----------------------------------------------------------------------------
template <class TFLOAT, class TUINT>
inline tMatrix<TFLOAT, TUINT> tMatrix<TFLOAT, TUINT>::operator-() const {
  tMatrix<TFLOAT, TUINT> B =
      tMatrix<TFLOAT, TUINT>((*this).rows(), (*this).cols());
  for (TUINT row = 0; row < rows(); row++) {
    for (TUINT col = 0; col < cols(); col++) {
      B(row, col) = (-1.0) * (*this)(row, col);
    }
  }
  return B;
}

template <class TFLOAT, class TUINT>
inline void tMatrix<TFLOAT, TUINT>::insertNext(TFLOAT value) {
  if (m_InsertIndex >= m_Rows * m_Cols) {
    throw cException(
        "illegal call tMatrix << value: tried to access out of bounds index",
        __FILE__, __LINE__);
  }
  m_Data[m_InsertIndex] = value;
  ++m_InsertIndex;
}

//! overloaded outputoperator
template <class TFLOAT, class TUINT>
inline std::ostream &operator<<(std::ostream &os,
                                const tMatrix<TFLOAT, TUINT> &M) {
  return M.write(os);
}

//! overloaded input operator for assigning values to the matrix
template <class TFLOAT1, class TFLOAT2, class TUINT>
inline tMatrix<TFLOAT1, TUINT> &operator<<(tMatrix<TFLOAT1, TUINT> &A,
                                           TFLOAT2 value) {
  A.insertNext(value);
  return A;
}

template <class TFLOAT, class TUINT>
bool tMatrix<TFLOAT, TUINT>::isSymmetric(const double &tol) const {
  if (rows() != cols()) return false;

  for (int r = 0; r < rows(); r++) {
    for (int c = r + 1; c < cols(); c++) {
      if (std::abs(operator()(r, c) - operator()(c, r)) > tol) return false;
    }
  }

  return true;
}

template <class TFLOAT, class TUINT>
double tMatrix<TFLOAT, TUINT>::getDeterminant() const {
  double determinant = 0.;
  double detfaktor = 1.;
  return computeDeterminant(determinant, detfaktor);
}

// ----------------------------------------------------------------------------
//   detA = |A|  // TODO: use Gauss elimination instead of Laplace's formular
// ----------------------------------------------------------------------------
template <class TFLOAT, class TUINT>
double tMatrix<TFLOAT, TUINT>::computeDeterminant(double determinant,
                                                  double detfaktor) const {
  if (rows() == 1 && cols() == 1) {
    return operator()(0, 0);
  }
  if (rows() == 2 && cols() == 2) {
    determinant +=
        detfaktor *
        (operator()(0, 0) * operator()(1, 1) - operator()(0, 1) * operator()(
                                                                      1, 0));
    // std::cout<<"det: "<<determinant<<" detf: "<<detfaktor<<"\n";
  } else {
    for (int k = 0; k < cols(); k++) {
      int pm = k % 2;
      double faktor = 0.0;
      if (pm != 0)
        faktor = -1. * operator()(0, k);
      else
        faktor = operator()(0, k);
      if (faktor != 0.0) {
        // std::cout<<"new matrix "<<rows()-1<<" "<<cols()-1<<"...";
        tMatrix<TFLOAT, TUINT> dummy(rows() - 1, cols() - 1);
        // std::cout<<"done\n";
        for (int j = 0; j < cols(); j++) {
          if (j < k) {
            for (int i = 1; i < rows(); i++) {
              dummy(i - 1, j) = operator()(i, j);
            }
          } else if (j > k) {
            for (int i = 1; i < cols(); i++) {
              dummy(i - 1, j - 1) = operator()(i, j);
            }
          }
        }
        // std::cout<<dummy;
        determinant = dummy.computeDeterminant(determinant, detfaktor * faktor);
      }
    }
  }

  return determinant;
}

// ----------------------------------------------------------------------------
//   returns the maximum absolute value of the matrix
// ----------------------------------------------------------------------------
template <class TFLOAT, class TUINT>
TFLOAT tMatrix<TFLOAT, TUINT>::getAbsMax(void) const {
  TFLOAT temp_val = 0.0;
  for (TUINT row = 0; row < this->rows(); row++) {
    for (TUINT col = 0; col < this->cols(); col++) {
      if (std::abs((*this)(0, 0)) > temp_val) {
        temp_val = std::abs((*this)(0, 0));
      }
    }
  }
  return temp_val;
}

template <class TFLOAT1, class TFLOAT2, class TUINT>
inline void scale(tMatrix<TFLOAT1, TUINT> &A, const TFLOAT2 &c) {
  for (TUINT z = 0; z < A.rows(); z++)
    for (TUINT s = 0; s < A.cols(); s++) A(z, s) *= c;
}

// ----------------------------------------------------------------------------
//   output of the matrix
// ----------------------------------------------------------------------------
template <class TFLOAT, class UINT>
std::ostream &tMatrix<TFLOAT, UINT>::write(std::ostream &os) const {
  UINT i, j;
  UINT sc, mc;    // first/last column
  UINT step = 4;  // only 4 columns

  int w = os.width();
  int p = os.precision();
  int myWidth = 14;

  os.setf(std::ios::scientific);
  sc = 0;
  while (sc < cols()) {
    mc = sc + step;
    if (mc > cols()) mc = cols();
    os << std::endl << "Spalten " << sc << " bis " << mc - 1 << std::endl;
    for (i = 0; i < rows(); i++) {
      for (j = sc; j < mc; j++) {
#ifdef PETSC_USE_COMPLEX
        os << "[";
        os << std::setw(myWidth - 1) << m_Data[i * cols() + j].real();
        os << std::setw(myWidth) << m_Data[i * cols() + j].imag();
        os << "]  ";
#else
        os << std::setw(myWidth) << m_Data[i * cols() + j];
#endif
      }
      os << std::endl;
    }
    sc += step;
  }
  os.unsetf(std::ios::scientific);

  os.width(w);
  os.precision(p);

  return os;
}

// ----------------------------------------------------------------------------
//   returns true if the matrix is equal to zero within relative tolerance
// ----------------------------------------------------------------------------
template <class TFLOAT, class TUINT>
bool tMatrix<TFLOAT, TUINT>::isZero(const double &tol) const {
  tMatrix<TFLOAT, TUINT> zeromatrix = tMatrix<TFLOAT, TUINT>(rows(), cols());
  return approx_equal((*this), zeromatrix, tol);
}

template <class TFLOAT, class TUINT>
bool tMatrix<TFLOAT, TUINT>::isInvertible(const double &tol) const {
  if (rows() != cols()) {
    return false;
  }
  double determinant = getDeterminant();
  if (approx_equal(determinant, 0.0, tol)) {  // determinant is zero
    return false;
  } else {
    return true;
  }
}

// ----------------------------------------------------------------------------
//   true for a square matrix that is not invertible
// ----------------------------------------------------------------------------
template <class TFLOAT, class TUINT>
bool tMatrix<TFLOAT, TUINT>::isSingular(const double &tol) const {
  if (rows() != cols()) {
    return false;
  }
  double determinant = getDeterminant();
  if (approx_equal(determinant, 0.0, tol)) {  // determinant is zero
    return true;
  } else {
    return false;
  }
}

template <class TFLOAT, class TUINT>
tMatrix<TFLOAT, TUINT> tMatrix<TFLOAT, TUINT>::inv(void) const {
  if (!isInvertible()) {
    throw cException("illegal call tMatrix.inv(): matrix is not invertible!",
                     __FILE__, __LINE__);
  }
  TUINT n_dim = cols();
  if (1 == n_dim) {
    tMatrix<TFLOAT, TUINT> tempmat = tMatrix<TFLOAT, TUINT>(1, 1);
    tempmat(0, 0) = 1.0 / operator()(0, 0);
    return tempmat;
  } else if (TUINT(2) == n_dim) {
    return invert2x2Matrix(*this);
  } else if (TUINT(3) == n_dim) {
    return invert3x3Matrix(*this);
  } else {
    throw cException(
        "illegal call tMatrix.inv(): matrix dimensions 4 or greater are "
        "currently not supported!",
        __FILE__, __LINE__);
  }
}

template <class TFLOAT, class TUINT>
inline tMatrix<TFLOAT, TUINT> invert2x2Matrix(const tMatrix<TFLOAT, TUINT> &A) {
  double det = A.getDeterminant();
  tMatrix<TFLOAT, TUINT> Ainv = tMatrix<TFLOAT, TUINT>(A.rows(), A.cols());
  Ainv(0, 0) = A(1, 1) / det;
  Ainv(0, 1) = -A(0, 1) / det;
  Ainv(1, 0) = -A(1, 0) / det;
  Ainv(1, 1) = A(0, 0) / det;
  return Ainv;
}

template <class TFLOAT, class TUINT>
inline tMatrix<TFLOAT, TUINT> invert3x3Matrix(const tMatrix<TFLOAT, TUINT> &A) {
  double det = A.getDeterminant();
  tMatrix<TFLOAT, TUINT> Ainv = tMatrix<TFLOAT, TUINT>(A.rows(), A.cols());
  // "adjugate / determinant" solution
  Ainv(0, 0) = (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) / det;
  Ainv(0, 1) = (-A(0, 1) * A(2, 2) + A(0, 2) * A(2, 1)) / det;
  Ainv(0, 2) = (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) / det;
  Ainv(1, 0) = (-A(1, 0) * A(2, 2) + A(1, 2) * A(2, 0)) / det;
  Ainv(1, 1) = (A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0)) / det;
  Ainv(1, 2) = (-A(0, 0) * A(1, 2) + A(0, 2) * A(1, 0)) / det;
  Ainv(2, 0) = (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0)) / det;
  Ainv(2, 1) = (-A(0, 0) * A(2, 1) + A(0, 1) * A(2, 0)) / det;
  Ainv(2, 2) = (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0)) / det;
  return Ainv;
}

// overloaded approx_equal for scalar values
// -> maybe move to new scalar.h class ?
// adapted from armadillo package (v. 8.300.0) "internal_approx_equal_rel_diff"
// function called from general arma::approx_equal function
template <class TFLOAT1, class TFLOAT2>
bool approx_equal(const TFLOAT1 &a, const TFLOAT2 &b, const double &tol) {
  if (a == b)  // handles exact zero and unlikely exact equality cases
  {            // therefore avoids zero division errors
    return true;
  }

  const double abs_a = std::abs(double(a));
  const double abs_b = std::abs(double(b));

  const double max_c = std::max(abs_a, abs_b);

  const double abs_d = std::abs(a - b);

  if (max_c >= 1.0) {
    if (abs_d > (tol * max_c)) {
      return false;
    }
  } else {
    if ((abs_d / max_c) > tol) {
      return false;
    }
  }

  return true;
}

template <class TFLOAT1, class TFLOAT2, class TUINT>
bool approx_equal(const tMatrix<TFLOAT1, TUINT> &A,
                  const tMatrix<TFLOAT2, TUINT> &B, const double &tol) {
  if (!isAddCompatible(A, B)) {
    return false;
  }

  for (TUINT row = 0; row < A.rows(); row++) {
    for (TUINT col = 0; col < A.cols(); col++) {
      const double a = A(row, col);
      const double b = B(row, col);

      if (false ==
          approx_equal(
              a, b,
              tol))  // call overloaded function for single float values a and b
      {
        return false;
      }
    }
  }

  return true;
}

// ----------------------------------------------------------------------------
//   A -> A^T
// ----------------------------------------------------------------------------
template <class TFLOAT, class TUINT>
inline void transpose(const tMatrix<TFLOAT, TUINT> &A,
                      tMatrix<TFLOAT, TUINT> &AT) {
  if (A.rows() == AT.cols() && A.cols() == AT.rows()) {
    for (TUINT zeile = 0; zeile < A.rows(); zeile++)
      for (TUINT spalte = 0; spalte < A.cols(); spalte++)
        AT(spalte, zeile) = A(zeile, spalte);
  } else {
    throw cException(
        "illegal call transpose A->A^T: improper size of matrices!", __FILE__,
        __LINE__);
  }
}

template <class TFLOAT, class TUINT>
inline tMatrix<TFLOAT, TUINT> tMatrix<TFLOAT, TUINT>::trans(void) const {
  tMatrix<TFLOAT, TUINT> AT = tMatrix<TFLOAT, TUINT>(cols(), rows());
  transpose((*this), AT);
  return AT;
}

// ----------------------------------------------------------------------------
//   Cholesky decomposition A=LL^T
// ----------------------------------------------------------------------------
template <class TFLOAT, class UINT>
inline void decompositionCholesky(const tMatrix<TFLOAT, UINT> &A,
                                  tMatrix<TFLOAT, UINT> &L) {
  // if(A.isSymmetric() == true)
  //{
  if (A.rows() == A.cols() && A.rows() == L.rows() && L.cols() == L.rows()) {
    //--- copy matrix A to L
    L = A;
    TFLOAT dummy;
    TFLOAT eps = 1.E-20;
    for (UINT row = 0; row < L.rows(); ++row) {
      for (UINT col = row; col < L.cols(); ++col) {
        dummy = L(row, col);
        for (UINT k = row - 1; k >= 0; --k) dummy -= L(row, k) * L(col, k);

        if (row == col)
          L(col, row) = sqrt(dummy);
        else
          L(col, row) = dummy / L(row, row);
      }
    }
  } else {
    throw cException(
        "illegal call Cholesky decomposition A=LL^T: improper size of "
        "matrices!",
        __FILE__, __LINE__);
  }
  //}
  // else
  //{
  //  throw cException("illegal call Cholesky decomposition A=LL^T: matrix is
  //  not symmetric!", __FILE__, __LINE__);
  //}
  //--- delete upper triangular
  for (UINT row = 0; row < L.rows(); ++row) {
    for (UINT col = row + 1; col < L.cols(); ++col) {
      L(row, col) = 0.0;
    }
  }
}

// ----------------------------------------------------------------------------
//   B += A
// ----------------------------------------------------------------------------
template <class TFLOAT1, class TFLOAT2, class UINT>
inline void add(const tMatrix<TFLOAT1, UINT> &A, tMatrix<TFLOAT2, UINT> &B) {
  if (B.rows() == A.rows() && B.cols() == A.cols()) {
    for (UINT zeile = 0; zeile < B.rows(); zeile++)
      for (UINT spalte = 0; spalte < B.cols(); spalte++)
        B(zeile, spalte) += A(zeile, spalte);
  } else {
    throw cException("illegal call add B += A : improper size of matrices!",
                     __FILE__, __LINE__);
  }
}

// ----------------------------------------------------------------------------
//   Check if A and B have the same dimensions
// ----------------------------------------------------------------------------
template <class TFLOAT1, class TFLOAT2, class UINT>
inline bool isAddCompatible(const tMatrix<TFLOAT1, UINT> &A,
                            const tMatrix<TFLOAT2, UINT> &B) {
  if (B.rows() == A.rows() && B.cols() == A.cols()) {
    return true;
  } else {
    return false;
  }
}

// ----------------------------------------------------------------------------
//   '+'-operator, returns C = A + B
// ----------------------------------------------------------------------------
template <class TFLOAT1, class TFLOAT2, class UINT>
inline tMatrix<TFLOAT1, UINT> operator+(const tMatrix<TFLOAT1, UINT> &A,
                                        const tMatrix<TFLOAT2, UINT> &B) {
  if (!isAddCompatible(A, B)) {
    throw cException("illegal call A + B operator : improper size of matrices!",
                     __FILE__, __LINE__);
  }
  tMatrix<TFLOAT1, UINT> C = tMatrix<TFLOAT1, UINT>(A.rows(), A.cols());
  add(A, C);
  add(B, C);
  return C;
}

// ----------------------------------------------------------------------------
//   '-'-operator, returns C = A - B
// ----------------------------------------------------------------------------
template <class TFLOAT1, class TFLOAT2, class UINT>
inline tMatrix<TFLOAT1, UINT> operator-(const tMatrix<TFLOAT1, UINT> &A,
                                        const tMatrix<TFLOAT2, UINT> &B) {
  if (!isAddCompatible(A, B)) {
    throw cException("illegal call A + B operator : improper size of matrices!",
                     __FILE__, __LINE__);
  }
  tMatrix<TFLOAT1, UINT> C = tMatrix<TFLOAT1, UINT>(A.rows(), A.cols());
  mult((-1.0), B, C);
  add(A, C);
  return C;
}

// ----------------------------------------------------------------------------
//   C += A*B   e.g. C(3,4)+=A(3,6)*B(6,4)
// ----------------------------------------------------------------------------
template <class TFLOAT1, class TFLOAT2, class TFLOAT3, class UINT>
inline void mult(const tMatrix<TFLOAT1, UINT> &A,
                 const tMatrix<TFLOAT2, UINT> &B, tMatrix<TFLOAT3, UINT> &C) {
  if (A.cols() == B.rows() && C.rows() == A.rows() && C.cols() == B.cols()) {
    for (UINT zeile = 0; zeile < C.rows(); zeile++)
      for (UINT spalte = 0; spalte < C.cols(); spalte++)
        for (UINT k = 0; k < A.cols(); k++)
          C(zeile, spalte) += A(zeile, k) * B(k, spalte);
  } else {
    throw cException("illegal call mult C += A*B : improper size of matrices!",
                     __FILE__, __LINE__);
  }
}

// ----------------------------------------------------------------------------
//   c += A*b   e.g. c(3)+=A(3,6)*b(6)
// ----------------------------------------------------------------------------
template <class TFLOAT1, class TFLOAT2, class TFLOAT3, class UINT>
inline void mult(const tMatrix<TFLOAT1, UINT> &A,
                 const tVector<TFLOAT2, UINT> &B, tVector<TFLOAT3, UINT> &C) {
  if (A.cols() == B.size() && C.size() == A.rows()) {
    for (UINT zeile = 0; zeile < C.size(); zeile++)
      for (UINT spalte = 0; spalte < A.cols(); spalte++)
        C[zeile] += A(zeile, spalte) * B[spalte];
  } else {
    throw cException(
        "illegal call mult c += A*b : improper size of matrix or vectors!",
        __FILE__, __LINE__);
  }
}

// ----------------------------------------------------------------------------
//   B += scalar*A
// ----------------------------------------------------------------------------
template <class TFLOAT1, class TFLOAT2, class TFLOAT3, class TUINT>
inline void mult(const TFLOAT1 &scalar, const tMatrix<TFLOAT2, TUINT> &A,
                 tMatrix<TFLOAT3, TUINT> &B) {
  if (!isAddCompatible(A, B)) {
    throw cException(
        "illegal call mult(scalar,A,B)-function : improper size of matrices A "
        "and B!",
        __FILE__, __LINE__);
  }
  for (TUINT row = 0; row < A.rows(); row++) {
    for (TUINT col = 0; col < A.cols(); col++) {
      B(row, col) += scalar * A(row, col);
    }
  }
}

// ----------------------------------------------------------------------------
//   Check if A.cols() == B.rows()
// ----------------------------------------------------------------------------
template <class TFLOAT1, class TFLOAT2, class TUINT>
inline bool isMultiplyCompatible(const tMatrix<TFLOAT1, TUINT> &A,
                                 const tMatrix<TFLOAT2, TUINT> &B) {
  if (A.cols() == B.rows()) {
    return true;
  } else {
    return false;
  }
}

// ----------------------------------------------------------------------------
//   '*'-operator, returns C = A * B
// ----------------------------------------------------------------------------
template <class TFLOAT1, class TFLOAT2, class TUINT>
inline tMatrix<TFLOAT1, TUINT> operator*(const tMatrix<TFLOAT1, TUINT> &A,
                                         const tMatrix<TFLOAT2, TUINT> &B) {
  if (!isMultiplyCompatible(A, B)) {
    throw cException("illegal call A * B operator : improper size of matrices!",
                     __FILE__, __LINE__);
  }
  tMatrix<TFLOAT1, TUINT> C = tMatrix<TFLOAT1, TUINT>(A.rows(), B.cols());
  mult(A, B, C);
  return C;
}

// ----------------------------------------------------------------------------
//   '*'-operator for scalar left multiplication, returns B = scalar * A
// ----------------------------------------------------------------------------
template <class TFLOAT1, class TFLOAT2, class TUINT>
inline tMatrix<TFLOAT1, TUINT> operator*(const TFLOAT2 &scalar,
                                         const tMatrix<TFLOAT1, TUINT> &A) {
  tMatrix<TFLOAT1, TUINT> B = tMatrix<TFLOAT1, TUINT>(A.rows(), A.cols());
  mult(scalar, A, B);
  return B;
}

// ----------------------------------------------------------------------------
//   '*'-operator for scalar right multiplication, returns B = A * scalar
// ----------------------------------------------------------------------------
template <class TFLOAT1, class TFLOAT2, class TUINT>
inline tMatrix<TFLOAT1, TUINT> operator*(const tMatrix<TFLOAT1, TUINT> &A,
                                         const TFLOAT2 &scalar) {
  tMatrix<TFLOAT1, TUINT> B = tMatrix<TFLOAT1, TUINT>(A.rows(), A.cols());
  mult(scalar, A, B);
  return B;
}

// ----------------------------------------------------------------------------
//   Matrix product of 2 vectors
//   C += b * d^T   e.g. C(3,4)+= b(3)*d(4)^T
// ----------------------------------------------------------------------------
template <class TFLOAT1, class TFLOAT2, class TFLOAT3, class UINT>
inline void mult(const tVector<TFLOAT1, UINT> &A,
                 const tVector<TFLOAT2, UINT> &B, tMatrix<TFLOAT3, UINT> &C) {
  if (A.size() == C.rows() && B.size() == C.cols()) {
    for (UINT zeile = 0; zeile < A.size(); zeile++)
      for (UINT spalte = 0; spalte < B.size(); spalte++)
        C(zeile, spalte) += A[zeile] * B[spalte];
  } else {
    throw cException(
        "illegal call mult C += a*b : improper size of matrix or vectors!",
        __FILE__, __LINE__);
  }
}

// ----------------------------------------------------------------------------
//   KK += B^T * C * B * wdJ  (useful for FEM)
// ----------------------------------------------------------------------------
template <class TFLOAT1, class TFLOAT2, class TFLOAT3, class TFLOAT4,
          class UINT>
inline void BT_C_B_wdJ(tMatrix<TFLOAT1, UINT> &KK,
                       const tMatrix<TFLOAT2, UINT> &B,
                       const tMatrix<TFLOAT3, UINT> &C, const TFLOAT4 &wdJ) {
  // --------------------------------------------------------------------------
  //   transpose B
  // --------------------------------------------------------------------------
  tMatrix<TFLOAT2, UINT> BT(B.cols(), B.rows());
  infam::transpose(B, BT);

  // --------------------------------------------------------------------------
  //   weight * detJ   einrechnen
  // --------------------------------------------------------------------------
  infam::scale(BT, wdJ);

  // --------------------------------------------------------------------------
  //   BT * C ( * wdJ )
  // --------------------------------------------------------------------------
  tMatrix<TFLOAT1, UINT> res1(BT.rows(), C.cols());
  infam::mult(BT, C, res1);

  // --------------------------------------------------------------------------
  //   BT * C * B * wdJ
  // --------------------------------------------------------------------------
  tMatrix<TFLOAT1, UINT> res2(KK.rows(), KK.cols());
  infam::mult(res1, B, res2);

  // --------------------------------------------------------------------------
  //   K += res
  // --------------------------------------------------------------------------
  infam::add(res2, KK);
}

// ----------------------------------------------------------------------------
//   KK += A^T * C * B * wdJ  (useful for SBFEM)
// ----------------------------------------------------------------------------
template <class TFLOAT1, class TFLOAT2, class TFLOAT3, class TFLOAT4,
          class TFLOAT5, class UINT>
inline void AT_C_B_wdJ(tMatrix<TFLOAT1, UINT> &KK,
                       const tMatrix<TFLOAT2, UINT> &A,
                       const tMatrix<TFLOAT3, UINT> &C,
                       const tMatrix<TFLOAT4, UINT> &B, const TFLOAT5 &wdJ) {
  // --------------------------------------------------------------------------
  //   transpose A
  // --------------------------------------------------------------------------
  tMatrix<TFLOAT2, UINT> AT(A.cols(), A.rows());
  infam::transpose(A, AT);

  // --------------------------------------------------------------------------
  //   weight * detJ   einrechnen
  // --------------------------------------------------------------------------
  infam::scale(AT, wdJ);

  // --------------------------------------------------------------------------
  //   AT * C ( * wdJ )
  // --------------------------------------------------------------------------
  tMatrix<TFLOAT1, UINT> res1(AT.rows(), C.cols());
  infam::mult(AT, C, res1);

  // --------------------------------------------------------------------------
  //   AT * C * B * wdJ
  // --------------------------------------------------------------------------
  tMatrix<TFLOAT1, UINT> res2(KK.rows(), KK.cols());
  infam::mult(res1, B, res2);

  // --------------------------------------------------------------------------
  //   K += res
  // --------------------------------------------------------------------------
  infam::add(res2, KK);
}

}  // namespace infam

#endif
