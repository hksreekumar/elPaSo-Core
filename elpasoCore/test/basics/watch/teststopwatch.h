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

#include "../../../source/basics/watch/stopwatch.h"

#ifdef HAVE_GTEST
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#endif  // HAVE_GTEST

//! @brief Mock class for cStopWatch
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
class cTestStopWatch : protected cStopwatch, public testing::Test {
 public:
  void testCorrectIntialization() {
    // arrange
    // act
    // assert
    ASSERT_EQ(m_StartWall.size(), cstNumTimingTasks);
    ASSERT_EQ(m_DurationWall.size(), cstNumTimingTasks);
  }

  void testCorrectStartStop() {
    // arrange
    // act
    startClock(Overall);
    stopClock(Overall);
    // assert
    char buff[258];
    convertToString(buff, Overall);
    ASSERT_TRUE(m_DurationWall[(int)Overall] >= 0);
  }
};

TEST_F(cTestStopWatch, correctInitalization) { testCorrectIntialization(); }

TEST_F(cTestStopWatch, correctStartStop) { testCorrectStartStop(); }