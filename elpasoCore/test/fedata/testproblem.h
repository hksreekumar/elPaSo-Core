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

#include "../../source/analysis/frequency/analysisfrequencybasic.h"
#include "../../source/fedata/problem.h"

#ifdef HAVE_GTEST
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#endif  // HAVE_GTEST

//! @brief Tests for cProblem - correct initialization
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(cProblem, correctInitialization) {
  // arrange
  // act
  cProblem testObject;
  // assert
  ASSERT_TRUE(testObject.getAnalysis() == NULL);
  ASSERT_TRUE(testObject.getMesh() != NULL);
}

//! @brief Tests for cProblem - correct copy constructor
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(cProblem, correctCopyConstructor) {
  // arrange
  cProblem testObject;
  // testObject.setAnalysis(new cAnalysisFrequencyBasic);
  // act
  cProblem testObjectCopied(testObject);
  // assert
  ASSERT_TRUE(testObjectCopied.getAnalysis() == NULL);
  ASSERT_TRUE(testObjectCopied.getMesh() != NULL);
}

//! @brief Tests for cProblem - correct setting of analysis object
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(cProblem, correctSetAnalysis) {
  // arrange
  cProblem testObject;
  // act
  testObject.setAnalysis(new cAnalysisFrequencyBasic);
  // assert
  ASSERT_TRUE(testObject.getAnalysis() != NULL);
}
//! @brief Tests for cProblem - correct setting of filename
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
TEST(cProblem, correctSetFilename) {
  // arrange
  cProblem testObject;
  // act
  testObject.setFilename("test.xyz");
  // assert
  ASSERT_EQ(testObject.getFilename(), "test.xyz");
}