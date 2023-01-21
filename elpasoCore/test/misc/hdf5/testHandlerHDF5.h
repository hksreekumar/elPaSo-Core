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

#include "../../../source/misc/hdf5/handlerhdf5.h"

#ifdef HAVE_GTEST
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#endif  // HAVE_GTEST

//! @brief A misc function to test as a handle for HDF5::Handles
//! @author Harikrishnan Sreekumar
//! @date 21.10.2020
void miscTestFunction1(void *) {
  // do nothing
}

//! @brief Another misc function to test as a handle for HDF5::Handles
//! @author Harikrishnan Sreekumar
//! @date 21.10.2020
void miscTestFunction2(void *) {
  // do nothing
}

//! @class cTestIOHDF5
//! @brief Test class for HDF5::Handles
//! @author Harikrishnan Sreekumar
//! @date 21.10.2020
class cTestIOHDF5 : public testing::Test {
 public:
  //! @brief Misc class member function to test as a handle for HDF5::Handles
  //! @author Harikrishnan Sreekumar
  //! @date 21.10.2020
  static void miscTestMemberFunction1(void *) {
    // do nothing
  }
};

namespace HDF5 {
//! @brief GTEST: Test to ensusre proper working of Handlers as list
//! @author Harikrishnan Sreekumar
//! @date 31.10.2020
TEST(Handler, checkHandlerHDF5WorksForExternalFunctions) {
  // arrange and act
  Handler handles[] = {HDF5::Handler("HandleFunction1", miscTestFunction1),
                       HDF5::Handler("HandleFunction2", miscTestFunction2)};
  // assert
  ASSERT_EQ(handles[0].functionpair.second, miscTestFunction1);
  ASSERT_NE(handles[0].functionpair.second, miscTestFunction2);
  ASSERT_EQ(handles[1].functionpair.second, miscTestFunction2);
  ASSERT_NE(handles[1].functionpair.second, miscTestFunction1);
}

//! @brief GTEST: Test to ensusre proper working of Handlers - where functions
//! are static class members
//! @author Harikrishnan Sreekumar
//! @date 31.10.2020
TEST(Handler, checkHandlerHDF5WorksForMemberFunctions) {
  // arrange and act
  HDF5::Handler handles[] = {HDF5::Handler(
      "HandleMemberFunction1", cTestIOHDF5::miscTestMemberFunction1)};
  // assert
  ASSERT_EQ(handles[0].functionpair.second,
            cTestIOHDF5::miscTestMemberFunction1);
}
}  // namespace HDF5