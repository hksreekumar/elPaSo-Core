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

//! @brief Mock class for cBoundaryConditionStructure
//! @author Akash Doddamane Lingaraja
//! @date 21.03.2022
class cTestBoundaryConditionStructure : private cFemParserHDF5,
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
        elpasoTestResourceDir + "unitHDF5/bc/" + HDF5_BC_STRUCTURE;

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

  //! @brief Function to assert - do check if dof is fixed or not for a specific
  //! node
  //! @author Akash Doddamane Lingaraja
  //! @date 21.03.2022
  void TEST_ASSERT_CHECKIFFIXED() {
    cBoundaryCondition* bc = myTestProblem.getMesh()->getBoundaryCondition(11);
    cBoundaryConditionStructure* bcStruct;
    bcStruct = dynamic_cast<cBoundaryConditionStructure*>(bc);

    for (int i = 0; i < 21; i++) {
      if (i < 6) {
        ASSERT_EQ(bcStruct->checkIfFixed(i), true);
      } else {
        ASSERT_EQ(bcStruct->checkIfFixed(i), false);
      }
    }
  }

  //! @brief Function to assert - do check if dof is fixed or not for a specific
  //! node using enum
  //! @author Akash Doddamane Lingaraja
  //! @date 21.03.2022
  void TEST_ASSERT_CHECKIFFIXED_ENUM() {
    cBoundaryCondition* bc = myTestProblem.getMesh()->getBoundaryCondition(11);
    cBoundaryConditionStructure* bcStruct;
    bcStruct = dynamic_cast<cBoundaryConditionStructure*>(bc);

    ASSERT_EQ(bcStruct->checkIfFixed(disp_x1), 1.);
    ASSERT_EQ(bcStruct->checkIfFixed(disp_x2), 1.);
    ASSERT_EQ(bcStruct->checkIfFixed(disp_x3), 1.);
    ASSERT_EQ(bcStruct->checkIfFixed(disp_w1), 1.);
    ASSERT_EQ(bcStruct->checkIfFixed(disp_w2), 1.);
    ASSERT_EQ(bcStruct->checkIfFixed(disp_w3), 1.);
    ASSERT_EQ(bcStruct->checkIfFixed(pore0), 0.);
    ASSERT_EQ(bcStruct->checkIfFixed(pore1), 0.);
    ASSERT_EQ(bcStruct->checkIfFixed(pore2), 0.);
    ASSERT_EQ(bcStruct->checkIfFixed(pore3), 0.);
    ASSERT_EQ(bcStruct->checkIfFixed(disp_xd3), 0.);
    ASSERT_EQ(bcStruct->checkIfFixed(disp_wd1), 0.);
    ASSERT_EQ(bcStruct->checkIfFixed(disp_wd2), 0.);
    ASSERT_EQ(bcStruct->checkIfFixed(disp_dwdx), 0.);
    ASSERT_EQ(bcStruct->checkIfFixed(disp_dwdy), 0.);
    ASSERT_EQ(bcStruct->checkIfFixed(disp_dwdxy), 0.);
    ASSERT_EQ(bcStruct->checkIfFixed(fluid), 0.);
    ASSERT_EQ(bcStruct->checkIfFixed(disp_z_1), 0.);
    ASSERT_EQ(bcStruct->checkIfFixed(disp_z_3), 0.);
    ASSERT_EQ(bcStruct->checkIfFixed(disp_x1_2), 0.);
    ASSERT_EQ(bcStruct->checkIfFixed(disp_x2_2), 0.);
  }

  //! @brief Function to assert - do get prescribed value of a specific degree
  //! of freedom
  //! @author Akash Doddamane Lingaraja
  //! @date 21.03.2022
  void TEST_ASSERT_GETPRESCRIBEDVALUE() {
    cBoundaryCondition* bc = myTestProblem.getMesh()->getBoundaryCondition(11);
    cBoundaryConditionStructure* bcStruct;
    bcStruct = dynamic_cast<cBoundaryConditionStructure*>(bc);

    ASSERT_EQ(bcStruct->getPrescribedValue(disp_x1, NULL, NULL), 0.);
    ASSERT_EQ(bcStruct->getPrescribedValue(disp_x2, NULL, NULL), 1.);
    ASSERT_EQ(bcStruct->getPrescribedValue(disp_x3, NULL, NULL), 2.);
    ASSERT_EQ(bcStruct->getPrescribedValue(disp_w1, NULL, NULL), 3.);
    ASSERT_EQ(bcStruct->getPrescribedValue(disp_w2, NULL, NULL), 4.);
    ASSERT_EQ(bcStruct->getPrescribedValue(disp_w3, NULL, NULL), 5.);
    ASSERT_EQ(bcStruct->getPrescribedValue(pore0, NULL, NULL), 0.);
    ASSERT_EQ(bcStruct->getPrescribedValue(pore1, NULL, NULL), 0.);
    ASSERT_EQ(bcStruct->getPrescribedValue(pore2, NULL, NULL), 0.);
    ASSERT_EQ(bcStruct->getPrescribedValue(pore3, NULL, NULL), 0.);
    ASSERT_EQ(bcStruct->getPrescribedValue(disp_xd3, NULL, NULL), 0.);
    ASSERT_EQ(bcStruct->getPrescribedValue(disp_wd1, NULL, NULL), 0.);
    ASSERT_EQ(bcStruct->getPrescribedValue(disp_wd2, NULL, NULL), 0.);
    ASSERT_EQ(bcStruct->getPrescribedValue(disp_dwdx, NULL, NULL), 0.);
    ASSERT_EQ(bcStruct->getPrescribedValue(disp_dwdy, NULL, NULL), 0.);
    ASSERT_EQ(bcStruct->getPrescribedValue(disp_dwdxy, NULL, NULL), 0.);
    ASSERT_EQ(bcStruct->getPrescribedValue(fluid, NULL, NULL), 0.);
    ASSERT_EQ(bcStruct->getPrescribedValue(disp_z_1, NULL, NULL), 0.);
    ASSERT_EQ(bcStruct->getPrescribedValue(disp_z_3, NULL, NULL), 0.);
    ASSERT_EQ(bcStruct->getPrescribedValue(disp_x1_2, NULL, NULL), 0.);
    ASSERT_EQ(bcStruct->getPrescribedValue(disp_x2_2, NULL, NULL), 0.);
  }

  //! @brief Function to assert - do set a specific degree of freedom to be
  //! fixed
  //! @author Akash Doddamane Lingaraja
  //! @date 21.03.2022
  void TEST_ASSERT_FIXDOF() {
    cBoundaryCondition* bc = myTestProblem.getMesh()->getBoundaryCondition(11);
    cBoundaryConditionStructure* bcStruct;
    bcStruct = dynamic_cast<cBoundaryConditionStructure*>(bc);
    bcStruct->fixDof(6);

    ASSERT_EQ(bcStruct->checkIfFixed(6), 1);
  }

  //! @brief Function to assert - do set a value of specific degree of freedom
  //! @author Akash Doddamane Lingaraja
  //! @date 21.03.2022
  void TEST_ASSERT_SETPRESCRIBEDVALUE() {
    cBoundaryCondition* bc = myTestProblem.getMesh()->getBoundaryCondition(11);
    cBoundaryConditionStructure* bcStruct;
    bcStruct = dynamic_cast<cBoundaryConditionStructure*>(bc);
    bcStruct->setPrescribedValue(1, 25);

    ASSERT_EQ(bcStruct->getPrescribedValue(disp_x2, NULL, NULL), 25.);
  }

  //! @brief Function to assert - do get array of fixed DOF
  //! @author Akash Doddamane Lingaraja
  //! @date 21.03.2022
  void TEST_ASSERT_GETFIXEDARRAY() {
    cBoundaryCondition* bc = myTestProblem.getMesh()->getBoundaryCondition(11);
    cBoundaryConditionStructure* bcStruct;
    bcStruct = dynamic_cast<cBoundaryConditionStructure*>(bc);

    std::bitset<21> refVal{"000000000000000111111"};
    std::bitset<21> val = bcStruct->getFixedArray();
    // std::cout << val << std::endl;

    ASSERT_EQ(refVal, bcStruct->getFixedArray());
  }

  //! @brief Function to assert - do get array of value of DOF
  //! @author Akash Doddamane Lingaraja
  //! @date 21.03.2022
  void TEST_ASSERT_GETVALUEARRAY() {
    cBoundaryCondition* bc = myTestProblem.getMesh()->getBoundaryCondition(11);
    cBoundaryConditionStructure* bcStruct;
    bcStruct = dynamic_cast<cBoundaryConditionStructure*>(bc);

    std::vector<PetscReal> val = bcStruct->getValueArray();
    std::vector<PetscReal> refVal = {0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    ASSERT_EQ((refVal == val), true);
  }

  //! @brief Function to assert - do get Identifier of the boundary condition
  //! @author Akash Doddamane Lingaraja
  //! @date 20.04.2022
  void TEST_ASSERT_GETIDENTIFIER() {
    cBoundaryCondition* bc = myTestProblem.getMesh()->getBoundaryCondition(11);
    cBoundaryConditionStructure* bcStruct;
    bcStruct = dynamic_cast<cBoundaryConditionStructure*>(bc);

    ASSERT_EQ(bcStruct->getIdentifier(), "structure");
  }
};

TEST_F(cTestBoundaryConditionStructure, checkIfDOFIsFixedOrNot) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_CHECKIFFIXED();
}

TEST_F(cTestBoundaryConditionStructure, checkIfDOFIsFixedOrNotUsingEnumerator) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_CHECKIFFIXED_ENUM();
}

TEST_F(cTestBoundaryConditionStructure,
       checkForCorrectPrescribedValueOfASpecificDOF) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_GETPRESCRIBEDVALUE();
}

TEST_F(cTestBoundaryConditionStructure, checkIfDOFIsFixedOrNotUsingFixDof) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_FIXDOF();
}

TEST_F(cTestBoundaryConditionStructure,
       checkIfCorrectPrescribedValueisSetForDOF) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_SETPRESCRIBEDVALUE();
}

TEST_F(cTestBoundaryConditionStructure, checkIfCorrectFixedArrayIsReturned) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_GETFIXEDARRAY();
}

TEST_F(cTestBoundaryConditionStructure, checkIfCorrectValueArrayIsReturned) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_GETVALUEARRAY();
}

TEST_F(cTestBoundaryConditionStructure, checkIfCorrectIdentifierIsReturned) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_GETIDENTIFIER();
}