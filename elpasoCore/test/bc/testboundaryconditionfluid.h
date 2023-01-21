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

#include "../../source/misc/parser/femparserhdf5.h"
#include "../resources/properties/testproperties.h"

#ifdef HAVE_GTEST
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#endif  // HAVE_GTEST

//! @brief Mock class for cBoundaryConditionFluid
//! @author Akash Doddamane Lingaraja
//! @date 08.01.2023
class cTestBoundaryConditionFluid : private cFemParserHDF5,
                                    public testing::Test {
 private:
  cProblem myTestProblem;

 public:
  //! @brief GTEST Setup
  //! @author Akash Doddamane Lingaraja
  //! @date 21.03.2022
  void SetUp() {
    // do nothing
  }

  //! @brief GTEST Teardown
  //! @author Akash Doddamane Lingaraja
  //! @date 21.03.2022
  void TearDown() {
    // close HDF5 File
#ifdef HAVE_HDF5
    cInputSingletonH5::getInstance()->closeContainer();
#endif
  }

  //! @brief Function to arrange - do parse the HDF5 file
  //! @author Akash Doddamane Lingaraja
  //! @date 21.03.2022
  void TEST_ARRANGE() {
    std::string elpasoTestResourceDir = ELPASO_TEST_RESOURCE_DIR;
    std::string testFileHDF5 =
        elpasoTestResourceDir + "unitHDF5/bc/" + HDF5_AF_LIN_UAF_ISO_DIR;

    // assign hdf5 to problem
    myTestProblem.setFilename(testFileHDF5);

    // open HDF5 File
#ifdef HAVE_HDF5
    cInputSingletonH5::getInstance()->setGlobalFileName(
        testFileHDF5);  // Set file name
    cInputSingletonH5::getInstance()->openContainer(ELPASO_H5_READONLY);
#endif
  }

  //! @brief Function to act - do get Node constraints from HDF5 file
  //! @author Akash Doddamane Lingaraja
  //! @date 21.03.2022
  void TEST_ACT() {
    cFemParserInterface::parseNodeEntity(myTestProblem);
    cFemParserInterface::parseNodeConstraintsEntity(myTestProblem);
  }

  //! @brief Function to assert - do check if dof is fixed or not
  //! @author Akash Doddamane Lingaraja
  //! @date 21.03.2022
  void TEST_ASSERT_CHECKIFFIXED() {
    cBoundaryCondition* bc = myTestProblem.getMesh()->getBoundaryCondition(17);
    cBoundaryConditionFluid* bcFluid;
    bcFluid = dynamic_cast<cBoundaryConditionFluid*>(bc);
    for (int i = 0; i < 21; i++) {
      if (i == 16) {
        ASSERT_EQ(bcFluid->checkIfFixed(i), true);
      } else {
        ASSERT_EQ(bcFluid->checkIfFixed(i), false);
      }
    }
  }

  //! @brief Function to assert - do check if dof is fixed or not using enum
  //! @author Akash Doddamane Lingaraja
  //! @date 21.03.2022
  void TEST_ASSERT_CHECKIFFIXED_ENUM() {
    cBoundaryCondition* bc = myTestProblem.getMesh()->getBoundaryCondition(17);
    cBoundaryConditionFluid* bcFluid;
    bcFluid = dynamic_cast<cBoundaryConditionFluid*>(bc);

    ASSERT_EQ(bcFluid->checkIfFixed(disp_x1), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(disp_x2), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(disp_x3), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(disp_w1), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(disp_w2), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(disp_w3), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(pore0), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(pore1), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(pore2), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(pore3), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(disp_xd3), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(disp_wd1), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(disp_wd2), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(disp_dwdx), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(disp_dwdy), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(disp_dwdxy), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(fluid), 1.);
    ASSERT_EQ(bcFluid->checkIfFixed(disp_z_1), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(disp_z_3), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(disp_x1_2), 0.);
    ASSERT_EQ(bcFluid->checkIfFixed(disp_x2_2), 0.);
  }

  //! @brief Function to assert - do get prescribed value of a specific degree
  //! of freedom
  //! @author Akash Doddamane Lingaraja
  //! @date 21.03.2022
  void TEST_ASSERT_GETPRESCRIBEDVALUE() {
    cBoundaryCondition* bc = myTestProblem.getMesh()->getBoundaryCondition(17);
    cBoundaryConditionFluid* bcFluid;
    bcFluid = dynamic_cast<cBoundaryConditionFluid*>(bc);

    ASSERT_EQ(bcFluid->getPrescribedValue(disp_x1, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(disp_x2, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(disp_x3, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(disp_w1, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(disp_w2, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(disp_w3, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(pore0, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(pore1, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(pore2, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(pore3, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(disp_xd3, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(disp_wd1, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(disp_wd2, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(disp_dwdx, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(disp_dwdy, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(disp_dwdxy, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(fluid, NULL, NULL), 1.0);
    ASSERT_EQ(bcFluid->getPrescribedValue(disp_z_1, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(disp_z_3, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(disp_x1_2, NULL, NULL), 0.);
    ASSERT_EQ(bcFluid->getPrescribedValue(disp_x2_2, NULL, NULL), 0.);
  }

  //! @brief Function to assert - do get pressure value
  //! @author Akash Doddamane Lingaraja
  //! @date 21.03.2022
  void TEST_ASSERT_SET_AND_GET_PRESSURE() {
    cBoundaryCondition* bc = myTestProblem.getMesh()->getBoundaryCondition(17);
    cBoundaryConditionFluid* bcFluid;
    bcFluid = dynamic_cast<cBoundaryConditionFluid*>(bc);
    bcFluid->setPressure(2.5);

    ASSERT_EQ(bcFluid->getPressure(), 2.5);
  }

  //! @brief Function to assert - do get Identifier of the boundary condition
  //! @author Akash Doddamane Lingaraja
  //! @date 20.04.2022
  void TEST_ASSERT_GETIDENTIFIER() {
    cBoundaryCondition* bc = myTestProblem.getMesh()->getBoundaryCondition(17);
    cBoundaryConditionFluid* bcFluid;
    bcFluid = dynamic_cast<cBoundaryConditionFluid*>(bc);
    bcFluid->setPressure(2.5);

    ASSERT_EQ(bcFluid->getIdentifier(), "acoustic");
  }
};

TEST_F(cTestBoundaryConditionFluid, checkIfDOFIsFixedOrNot) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_CHECKIFFIXED();
}

TEST_F(cTestBoundaryConditionFluid, checkIfDOFIsFixedOrNotUsingEnumerator) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_CHECKIFFIXED_ENUM();
}

TEST_F(cTestBoundaryConditionFluid,
       checkForCorrectPrescribedValueOfASpecificDOF) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_GETPRESCRIBEDVALUE();
}

TEST_F(cTestBoundaryConditionFluid, checkIfCorrectPressureValueIsSetAndGet) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_SET_AND_GET_PRESSURE();
}

TEST_F(cTestBoundaryConditionFluid, checkIfCorrectIdentifierIsReturned) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_GETIDENTIFIER();
}