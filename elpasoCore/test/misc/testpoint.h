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

#include "../../source/misc/point.h"

#ifdef HAVE_GTEST
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#endif  // HAVE_GTEST

//! @brief Test for cPoint correct initialization
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(cPoint, properConstructorIntialization) {
  // arrange
  // act
  cPoint testObject;
  // assert
  ASSERT_EQ(testObject.getComponent(0), 0.);
  ASSERT_EQ(testObject.getComponent(1), 0.);
  ASSERT_EQ(testObject.getComponent(2), 0.);
}

//! @brief Test for cPoint correct copy constructor
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(cPoint, properCopyConstructorIntialization) {
  // arrange
  cPoint testObject1;
  testObject1[0] = 1.;
  testObject1[1] = 2.;
  testObject1[2] = 3.;
  // act
  cPoint testObject2(testObject1);
  // assert
  ASSERT_EQ(testObject2.getComponent(0), 1.);
  ASSERT_EQ(testObject2.getComponent(1), 2.);
  ASSERT_EQ(testObject2.getComponent(2), 3.);
}

//! @brief Test for cPoint correct distance computation
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(cPoint, properDistanceComputation) {
  // arrange
  cPoint testObject1;
  testObject1[0] = 1.;
  testObject1[1] = 2.;
  testObject1[2] = 3.;
  cPoint testObject2;
  testObject1[0] = 4.;
  testObject1[1] = 5.;
  testObject1[2] = 6.;
  // act
  PetscReal dist = testObject1.distance(testObject2);
  // assert
  ASSERT_NEAR(dist, 8.77496, 1e-4);
}

//! @brief Test for cPoint correct distance2 computation
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(cPoint, properDistance2Computation) {
  // arrange
  cPoint testObject1;
  testObject1[0] = 1.;
  testObject1[1] = 2.;
  testObject1[2] = 3.;
  cPoint testObject2;
  testObject1[0] = 4.;
  testObject1[1] = 5.;
  testObject1[2] = 6.;
  // act
  PetscReal dist = testObject1.distance2(testObject2);
  // assert
  ASSERT_NEAR(dist, 8.77496 * 8.77496, 1e-4);
}