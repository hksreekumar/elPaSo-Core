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

#include "../../source/math/vektor.h"
#include "../../source/misc/mytypes.h"

#ifdef HAVE_GTEST
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#endif  // HAVE_GTEST

//! @brief Function to test vector L2 norm computation
//! @author Harikrishnan Sreekumar
//! @date 27.07.2022
TEST_F(tVector, correctVectorL2Norm) {
  // arrange
  // form a vector [1 11 -4 -1 7]
  infam::tVector<PetscReal, PetscInt> m_V(5);
  m_V[0] = 1;
  m_V[1] = 11;
  m_V[2] = -4;
  m_V[3] = -1;
  m_V[4] = 7;

  Vec m_VRef;
  VecCreate(PETSC_COMM_WORLD, &m_VRef);
  VecSetSizes(m_VRef, 5, PETSC_DECIDE);
  VecSetUp(m_VRef);

  VecSetValue(m_VRef, 0, 1., INSERT_VALUES);
  VecSetValue(m_VRef, 1, 11., INSERT_VALUES);
  VecSetValue(m_VRef, 2, -4., INSERT_VALUES);
  VecSetValue(m_VRef, 3, -1., INSERT_VALUES);
  VecSetValue(m_VRef, 4, 7., INSERT_VALUES);

  VecAssemblyBegin(m_VRef);
  VecAssemblyEnd(m_VRef);

  // act
  PetscScalar cmp_norm = m_V.abs();
  PetscScalar cmp_norm2 = m_V.abs2();

  // assert
  PetscReal ref_norm;
  VecNorm(m_VRef, NORM_2, &ref_norm);
  ASSERT_EQ(ref_norm, cmp_norm);
  ASSERT_EQ(ref_norm * ref_norm, cmp_norm2);
}

//! @brief Function to test vector intialization
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(tVector, correctInitialization) {
  // arrange
  // act
  infam::tVector<PetscReal, PetscInt> testVec(5);
  // assert
  ASSERT_EQ(testVec.size(), 5);
}

//! @brief Function to test vector resizing
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(tVector, correctVectorResize) {
  // arrange
  infam::tVector<PetscReal, PetscInt> testVec(5);
  // act
  testVec.resize(10);
  // assert
  ASSERT_EQ(testVec.size(), 10);
}

//! @brief Function to test vector indexing operation
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(tVector, correctVectorIndexingOperator) {
  // arrange
  infam::tVector<PetscReal, PetscInt> testVec(5);
  // act
  testVec[0] = 5.;
  // assert
  ASSERT_EQ(testVec[0], 5.);
}

//! @brief Function to test vector set value
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(tVector, correctVectorSetValue) {
  // arrange
  infam::tVector<PetscReal, PetscInt> testVec(5);
  // act
  testVec.setValue(500.1);
  // assert
  for (size_t i = 0; i < testVec.size(); i++) ASSERT_EQ(testVec[i], 500.1);
}

//! @brief Function to test vector dot product
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(tVector, correctVectorDotProduct) {
  // arrange
  infam::tVector<PetscReal, PetscInt> testVec1(5);
  testVec1.setValue(2.);
  infam::tVector<PetscReal, PetscInt> testVec2(5);
  testVec2.setValue(3.);
  // act
  PetscReal pdt = testVec1.dot(testVec2);
  // assert
  ASSERT_EQ(pdt, 30);
}
//! @brief Function to test vector cross product
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(tVector, correctVectorCrossProduct) {
  // arrange
  infam::tVector<PetscReal, PetscInt> testVec1(3);
  testVec1[0] = 2;
  testVec1[1] = 3;
  testVec1[2] = 4;
  infam::tVector<PetscReal, PetscInt> testVec2(3);
  testVec2[0] = 5;
  testVec2[1] = 6;
  testVec2[2] = 7;
  infam::tVector<PetscReal, PetscInt> testVec3(3);
  // act
  cross_product(testVec1, testVec2, testVec3);
  // assert
  ASSERT_EQ(testVec3[0], -3);
  ASSERT_EQ(testVec3[1], 6);
  ASSERT_EQ(testVec3[2], -3);
}
//! @brief Function to test vector scaling
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(tVector, correctVectorScale) {
  // arrange
  infam::tVector<PetscReal, PetscInt> testVec(5);
  testVec.setValue(5.);
  // act
  scale(testVec, 2.);
  // assert
  for (size_t i = 0; i < testVec.size(); i++) ASSERT_EQ(testVec[i], 10.);
}
//! @brief Function to test vector sorting
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(tVector, correctVectorSort) {
  // arrange
  infam::tVector<PetscReal, PetscInt> testVec(3);
  testVec[0] = 3.;
  testVec[1] = 2.;
  testVec[2] = 1.;
  // act
  testVec.sort();
  // assert
  ASSERT_EQ(testVec[0], 1.);
  ASSERT_EQ(testVec[1], 2.);
  ASSERT_EQ(testVec[2], 3.);
}
//! @brief Function to test vector minimum and maximum
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(tVector, correctVectorMinimumMaximum) {
  // arrange
  infam::tVector<PetscReal, PetscInt> testVec(3);
  testVec[0] = 3.;
  testVec[1] = 2.;
  testVec[2] = 1.;
  // act
  // assert
  ASSERT_EQ(testVec.getMin(), 1.);
  ASSERT_EQ(testVec.getMax(), 3.);
}