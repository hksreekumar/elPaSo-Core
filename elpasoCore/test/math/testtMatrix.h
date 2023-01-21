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

#include "../../source/math/matrix.h"
#include "../../source/misc/mytypes.h"

#ifdef HAVE_GTEST
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#endif  // HAVE_GTEST

//! @brief Class for testing cMatrix
//! @author Harikrishnan Sreekumar
//! @date 27.07.2022
class cTesttMatrix : public testing::Test {
 public:
  void SetUp() override {}

  void TearDown() override {}

  //! @brief Function to test determinant computation
  //! @author Harikrishnan Sreekumar
  //! @date 27.07.2022
  void TEST_3X3MATRIX_DETERMINANT() {
    // arrange
    // form a matrix [1 2 3; 4 5 6; 7 8 9] -> det is 0
    // form a matrix [1 2 3; 4 5 6; 7 8 0] -> det is 27
    cMatrix m_A(3, 3);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) m_A(i, j) = i * 3 + j + 1;

    // act
    PetscScalar det_case1 = m_A.getDeterminant();
    m_A(2, 2) = 0.;
    PetscScalar det_case2 = m_A.getDeterminant();

    // assert
    ASSERT_EQ(det_case1, 0.);
    ASSERT_EQ(det_case2, 27.);
  }

  //! @brief Function to test 3x3 inverse computation
  //! @author Harikrishnan Sreekumar
  //! @date 27.07.2022
  void TEST_3X3MATRIX_SOLVE() {
    // arrange
    // form a matrix [1 2 3; 4 5 6; 7 8 0] -> inv is [-1.77777777777778
    // 0.888888888888889-0.111111111111111; 1.55555555555556 - 0.777777777777778
    // 0.222222222222222; -0.111111111111111 0.222222222222222
    // -0.111111111111111]
    cMatrix m_A(3, 3);
    cMatrix m_Ainv(3, 3);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) m_A(i, j) = i * 3 + j + 1;

    m_A(2, 2) = 0.;
    std::vector<PetscScalar> referenceInverse{
        -1.77777777777778,  0.888888888888889,  -0.111111111111111,
        1.55555555555556,   -0.777777777777778, 0.222222222222222,
        -0.111111111111111, 0.222222222222222,  -0.111111111111111};

    // act
    m_Ainv = m_A.inv();

    // assert
    for (size_t i = 0; i < 3; i++)
      for (size_t j = 0; j < 3; j++)
        ASSERT_NEAR(m_Ainv(i, j), referenceInverse[i * 3 + j].real(), 1e-7);
  }

  //! @brief Function to test 2x2 inverse computation
  //! @author Harikrishnan Sreekumar
  //! @date 27.07.2022
  void TEST_2X2MATRIX_SOLVE() {
    // arrange
    // form a matrix [1 2; 3 4] -> inv is
    // [-2.00000000000000, 1.00000000000000, 1.50000000000000,
    // -0.500000000000000]
    cMatrix m_A(2, 2);
    cMatrix m_Ainv(2, 2);
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++) m_A(i, j) = i * 2 + j + 1;

    std::vector<PetscScalar> referenceInverse{
        -2.00000000000000, 1.00000000000000, 1.50000000000000,
        -0.500000000000000};

    // act
    m_Ainv = m_A.inv();

    // assert
    for (size_t i = 0; i < 2; i++)
      for (size_t j = 0; j < 2; j++)
        ASSERT_NEAR(m_Ainv(i, j), referenceInverse[i * 2 + j].real(), 1e-7);
  }

  //! @brief Function to test transpose
  //! @author Harikrishnan Sreekumar
  //! @date 27.07.2022
  void TEST_3X3MATRIX_TRANSPOSE() {
    // arrange
    // form a matrix [1 2 3; 4 5 6; 7 8 9] -> transpose is [1, 4, 7, 2, 5, 8, 3,
    // 6, 9]
    cMatrix m_A(3, 3);
    cMatrix m_Ainv(3, 3);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) m_A(i, j) = i * 3 + j + 1;

    std::vector<PetscScalar> referenceInverse{1, 4, 7, 2, 5, 8, 3, 6, 9};

    // act
    m_Ainv = m_A.trans();

    // assert
    for (size_t i = 0; i < 3; i++)
      for (size_t j = 0; j < 3; j++)
        ASSERT_EQ(m_Ainv(i, j), referenceInverse[i * 3 + j]);
  }

  //! @brief Function to test matrix-vector multiply
  //! @author Harikrishnan Sreekumar
  //! @date 27.07.2022
  void TEST_MATRIXVECTORMULTIPY() {
    // arrange
    // form a matrix [1 2; 3 4]
    cMatrix m_A(2, 2);
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++) m_A(i, j) = i * 2 + j + 1;
    // form a vector [5 6]
    cElementVector m_V(2);
    m_V[0] = 5;
    m_V[1] = 6;

    // act
    cElementVector pdt(2);
    infam::mult(m_A, m_V, pdt);

    // assert
    ASSERT_EQ(pdt[0], 17.);
    ASSERT_EQ(pdt[1], 39.);
  }

  //! @brief Function to test matrix-matrix multiply
  //! @author Harikrishnan Sreekumar
  //! @date 27.07.2022
  void TEST_MATRIXMATRIXMULTIPY() {
    // arrange
    // form a matrix [1 2; 3 4]
    cMatrix m_A(2, 2);
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++) m_A(i, j) = i * 2 + j + 1;
    // form a matrix [1 2 3; 4 5 6]
    cMatrix m_B(2, 3);
    m_B(0, 0) = 1;
    m_B(0, 1) = 2;
    m_B(0, 2) = 3;
    m_B(1, 0) = 4;
    m_B(1, 1) = 5;
    m_B(1, 2) = 6;

    // act
    cMatrix pdt(2, 3);
    infam::mult(m_A, m_B, pdt);

    // assert
    ASSERT_EQ(pdt(0, 0), 9.);
    ASSERT_EQ(pdt(0, 1), 12.);
    ASSERT_EQ(pdt(0, 2), 15.);
    ASSERT_EQ(pdt(1, 0), 19.);
    ASSERT_EQ(pdt(1, 1), 26.);
    ASSERT_EQ(pdt(1, 2), 33.);
  }

  //! @brief Function to test matrix-scalar multiply
  //! @author Harikrishnan Sreekumar
  //! @date 27.07.2022
  void TEST_MATRIXSCALARMULTIPY() {
    // arrange
    // form a matrix [1 2; 3 4]
    cMatrix m_A(2, 2);
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++) m_A(i, j) = i * 2 + j + 1;

    // act
    cMatrix pdt(2, 2);
    infam::mult(2., m_A, pdt);

    // assert
    ASSERT_EQ(pdt(0, 0), 1. * 2.);
    ASSERT_EQ(pdt(0, 1), 2. * 2.);
    ASSERT_EQ(pdt(1, 0), 3. * 2.);
    ASSERT_EQ(pdt(1, 1), 4. * 2.);
  }

  //! @brief Function to test matrix-matrix addition
  //! @author Harikrishnan Sreekumar
  //! @date 27.07.2022

  void TEST_MATRIXMATRIXADDITION() {
    // arrange
    // form a matrix [1 2; 3 4]
    cMatrix m_A(2, 2), m_B(2, 2);
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++) m_A(i, j) = i * 2 + j + 1;
    m_B = m_A;

    // act
    infam::add(m_A, m_B);

    // assert
    ASSERT_EQ(m_B(0, 0), 1. * 2.);
    ASSERT_EQ(m_B(0, 1), 2. * 2.);
    ASSERT_EQ(m_B(1, 0), 3. * 2.);
    ASSERT_EQ(m_B(1, 1), 4. * 2.);
  }

 private:
};

TEST_F(cTesttMatrix, when3x3MatrixPassedCorrectInverseComputation) {
  TEST_3X3MATRIX_SOLVE();
}

TEST_F(cTesttMatrix, when3x3MatrixPassedCorrectDeterminantComputation) {
  TEST_3X3MATRIX_DETERMINANT();
}

TEST_F(cTesttMatrix, when2x2MatrixPassedCorrectDeterminantComputation) {
  TEST_2X2MATRIX_SOLVE();
}

TEST_F(cTesttMatrix, whenMatrixPassedCorrectTransposeComputation) {
  TEST_3X3MATRIX_TRANSPOSE();
}

TEST_F(cTesttMatrix, correctMatrixVectorMultipy) { TEST_MATRIXVECTORMULTIPY(); }

TEST_F(cTesttMatrix, correctMatrixMatrixMultipy) { TEST_MATRIXMATRIXMULTIPY(); }

TEST_F(cTesttMatrix, correctMatrixScalarMultipy) { TEST_MATRIXSCALARMULTIPY(); }

TEST_F(cTesttMatrix, correctMatrixMatrixAddition) {
  TEST_MATRIXMATRIXADDITION();
}
