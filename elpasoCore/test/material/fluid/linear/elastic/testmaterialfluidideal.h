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

#include "../../../../../source/misc/parser/femparserhdf5.h"
#include "../../../../resources/properties/testproperties.h"

#ifdef HAVE_GTEST
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#endif  // HAVE_GTEST

//! @brief Mock class for cMaterialFluidIdeal
//! @author Akash Doddamane Lingaraja
//! @date 01.03.2022
class cTestMaterialFluidIdeal : private cFemParserHDF5, public testing::Test {
 private:
  cProblem myTestProblem;

 public:
  //! @brief GTEST Setup
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void SetUp() {
    // do nothing
  }

  //! @brief GTEST Teardown
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TearDown() {
    // close HDF5 File
#ifdef HAVE_HDF5
    cInputSingletonH5::getInstance()->closeContainer();
#endif
  }

  //! @brief Function to arrange - do parse the HDF5 file
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ARRANGE() {
    std::string elpasoTestResourceDir = ELPASO_TEST_RESOURCE_DIR;
    std::string testFileHDF5 = elpasoTestResourceDir +
                               "unitHDF5/material/fluid/linear/elastic/" +
                               HDF5_AF_LIN_UAF_ISO_DIR;

    // assign hdf5 to problem
    myTestProblem.setFilename(testFileHDF5);

    // open HDF5 File
#ifdef HAVE_HDF5
    cInputSingletonH5::getInstance()->setGlobalFileName(
        testFileHDF5);  // Set file name
    cInputSingletonH5::getInstance()->openContainer(ELPASO_H5_READONLY);
#endif
  }

  //! @brief Function to act - do get Material properties
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ACT() {
    cFemParserInterface::parseNodeEntity(myTestProblem);
    cFemParserInterface::parseMaterialsEntity(myTestProblem);
    cFemParserInterface::parseElementsEntity(myTestProblem);

    for (ItMapElements it = myTestProblem.getMesh()->getFirstElement();
         it != myTestProblem.getMesh()->getLastElement(); it++) {
      // set omega to 1.0 as poroelastic materials otherwise don't compute
      // matrices
      it->second->getMaterial()->setOmega(1.0);
      it->second->getMaterial()->updateMaterial();
    }
  }

  //! @brief Function to assert - do get frequency dependent wave number.
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ASSERT_GETCFOMEGA() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialFluidIdeal* mat;
    mat = dynamic_cast<cMaterialFluidIdeal*>(ptrMaterial);

    ASSERT_EQ(mat->getCfOmega(), 343.0);
  }

  //! @brief Function to assert - do get frequency dependent density
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ASSERT_GETRHOOMEGA() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialFluidIdeal* mat;
    mat = dynamic_cast<cMaterialFluidIdeal*>(ptrMaterial);

    ASSERT_EQ(mat->getRhoOmega(), 1.21);
  }
};

TEST_F(cTestMaterialFluidIdeal,
       checkForCorrectLinearViscoElasticFluidRealFreqDependentWaveNumber) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_GETCFOMEGA();
}

TEST_F(cTestMaterialFluidIdeal,
       checkForCorrectLinearViscoElasticFluidRealFreqDependentDensity) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_GETRHOOMEGA();
}