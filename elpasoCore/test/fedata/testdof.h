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

#include "../../source/fedata/dof.h"

#ifdef HAVE_GTEST
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#endif  // HAVE_GTEST

//! @brief Tests for cDegreesOfFreedom - correct initialization
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(cDegreesOfFreedom, correctConstructorInitialization) {
  // arrange
  // act
  cDegreesOfFreedom testObject;
  // assert
  ASSERT_EQ(testObject.getAllGlobalRows().size(), cstNumberOfKnownDofs);
  ASSERT_EQ(testObject.getAllDofFlags().size(), cstNumberOfKnownDofs);
}

//! @brief Tests for cDegreesOfFreedom - correct copy constructor
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(cDegreesOfFreedom, correctCopyConstructorInitialization) {
  // arrange
  cDegreesOfFreedom testObject1;
  // act
  cDegreesOfFreedom testObject2(testObject1);
  // assert
  ASSERT_EQ(testObject2.getAllGlobalRows().size(),
            testObject1.getAllGlobalRows().size());
  ASSERT_EQ(testObject2.getAllDofFlags().size(),
            testObject1.getAllDofFlags().size());
}

//! @brief Tests for cDegreesOfFreedom - correct activation of dofs
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(cDegreesOfFreedom, correctActivateDofFunction) {
  // arrange
  cDegreesOfFreedom testObject;
  // act
  testObject.activateDof(10);
  // assert
  ASSERT_EQ(testObject.checkIfActive(10), 1);
}

//! @brief Tests for cDegreesOfFreedom - correct deactivation of dofs
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(cDegreesOfFreedom, correctDeactivateDofFunction) {
  // arrange
  cDegreesOfFreedom testObject;
  // act
  testObject.deactivateDof(10);
  // assert
  ASSERT_EQ(testObject.checkIfActive(10), 0);
}

//! @brief Tests for cDegreesOfFreedom - correct setting of global row
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(cDegreesOfFreedom, correctSetGlobalRowFunction) {
  // arrange
  cDegreesOfFreedom testObject;
  // act
  testObject.setGlobalRow(1, 0);
  testObject.setGlobalRow(10, 1);
  testObject.setGlobalRow(2, 2);
  // assert
  ASSERT_EQ(testObject.getGlobalRow(1), 0);
  ASSERT_EQ(testObject.getGlobalRow(10), 1);
  ASSERT_EQ(testObject.getGlobalRow(2), 2);
}