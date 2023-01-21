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

#include "../../../source/fedata/problem.h"
#include "../../../source/misc/hdf5/inputh5.h"
#include "../../../source/misc/hdf5/outputh5.h"
#include "../../../source/misc/parser/femparserhdf5.h"
#include "../../../source/misc/parser/femparserinterface.h"
#include "../../resources/properties/testproperties.h"


#ifdef HAVE_GTEST
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#endif  // HAVE_GTEST

//! @class cTestFemParserHDF5
//! @brief Testing frame for cFemParserHDF5
//! @author Harikrishnan Sreekumar
//! @date 31.10.2020
class cTestFemParserHDF5 : public cFemParserHDF5, public testing::Test {
 public:
  //! @brief GTEST Setup
  //! @author Harikrishnan Sreekumar
  //! @date 31.10.2020
  void SetUp() override {
    // do nothing
  }

  //! @brief GTEST Teardown
  //! @author Harikrishnan Sreekumar
  //! @date 31.10.2020
  void TearDown() override {
    /// close HDF5 File
#ifdef HAVE_HDF5
    cInputSingletonH5::getInstance()->closeContainer();
#endif
  }

  //! @brief Function to arrange
  //! @author Harikrishnan Sreekumar
  //! @date 31.10.2020
  void TEST_ARRANGE(std::string _filename) {
    // Arrange
    /// test hdf5 file
    std::string elpasoTestResourceDir = ELPASO_TEST_RESOURCE_DIR;
    std::string testFileHDF5 = elpasoTestResourceDir + "unitHDF5/" + _filename;

    /// assign hdf5 to problem
    myTestProblem.setFilename(testFileHDF5);

    /// open HDF5 File
#ifdef HAVE_HDF5
    cInputSingletonH5::getInstance()->setGlobalFileName(
        testFileHDF5);  // Set file name
    cInputSingletonH5::getInstance()->openContainer(ELPASO_H5_READONLY);
#endif
  }

  //! @brief Function to act - frequency analysis
  //! @author Harikrishnan Sreekumar
  //! @date 31.10.2020
  void TEST_FREQUENCYANALYSIS_ACT() {
    // act
    /// perform HDF5 parse for analysis -> frequency
    cFemParserInterface::parseAnalysisEntity(myTestProblem);
  }

  //! @brief Function to assert - frequency analysis
  //! @author Harikrishnan Sreekumar
  //! @date 31.10.2020
  void TEST_FREQUENCYANALYSIS_ASSERT() {
    // assert
    cAnalysisFrequencyBasic* testAnalysisFrequency =
        dynamic_cast<cAnalysisFrequencyBasic*>(myTestProblem.getAnalysis());
    ASSERT_EQ(testAnalysisFrequency->getStartFrequency(),
              10);  // Assert for starting frequency
    ASSERT_EQ(testAnalysisFrequency->getIncrement(), 10);  // Assert for delta
    ASSERT_EQ(testAnalysisFrequency->getNumberOfSteps(),
              1);  // Assert for number of steps
  }

  //! @brief Function to act - nodes
  //! @author Harikrishnan Sreekumar
  //! @date 31.10.2020
  void TEST_NODES_ACT() {
    // act
    cFemParserInterface::parseNodeEntity(myTestProblem);
  }

  //! @brief Function to assert - nodes
  //! @author Harikrishnan Sreekumar
  //! @date 31.10.2020
  void TEST_NODES_ASSERT() {
    // assert
    /// if parsed data is correct
    ASSERT_EQ(mEstimatedNumberOfEntities.Nodes, 6);
    ASSERT_EQ(myTestProblem.getMesh()->getNumberOfNodes(),
              mEstimatedNumberOfEntities.Nodes);
    for (size_t i = 1; i <= 6; i++)
      ASSERT_EQ(myTestProblem.getMesh()->nodeInMesh(i), true);
    ASSERT_NE(myTestProblem.getMesh()->nodeInMesh(100), true);
  }

  //! @brief Function to act - elements
  //! @author Harikrishnan Sreekumar
  //! @date 31.10.2020
  void TEST_ELEMENTS_ACT() {
    // act
    /// perform HDF5 parse
    cFemParserInterface::parseNodeEntity(myTestProblem);
    TEST_MATERIALS_ACT();
    cFemParserInterface::parseElementsEntity(myTestProblem);
  }

  //! @brief Function to assert - elements
  //! @author Harikrishnan Sreekumar
  //! @date 31.10.2020
  void TEST_ELEMENTS_ASSERT() {
    // assert
    /// assert if parsed data is correct
    ASSERT_EQ(mEstimatedNumberOfEntities.AttributesExistInInputFile, true);
    ASSERT_EQ(mEstimatedNumberOfEntities.Elements, 2);
    ASSERT_EQ(m_SectionsOrientation, "global");
    ASSERT_EQ(m_SectionsOrientationFilename, "");
    ASSERT_EQ(m_SectionsMaterial, 1);
    ASSERT_EQ(myTestProblem.getMesh()->getNumberOfElements(),
              mEstimatedNumberOfEntities.Elements);
  }

  //! @brief Function to act - materials
  //! @author Harikrishnan Sreekumar
  //! @date 26.04.2021
  void TEST_MATERIALS_ACT() {
    // act
    cFemParserInterface::parseMaterialsEntity(myTestProblem);
  }

  //! @brief Function to assert - materials
  //! @author Harikrishnan Sreekumar
  //! @date 26.04.2021
  void TEST_MATERIALS_ASSERT() {
    // assert
    /// if parsed data is correct
    ASSERT_EQ(mEstimatedNumberOfEntities.Materials, 1);
    ASSERT_EQ(myTestProblem.getMesh()->getMaterial(1)->getIdentifier(), "name");
    ASSERT_EQ(myTestProblem.getMesh()->getMaterial(1)->getRho(), 2700.);
  }

  //! @brief Function to act - NodalLoads
  //! @author Harikrishnan Sreekumar
  //! @date 11.05.2021
  void TEST_NODALLOADS_ACT() {
    // extendended arrange
    cFemParserInterface::parseNodeEntity(myTestProblem);
    cFemParserInterface::parseMaterialsEntity(myTestProblem);
    cFemParserInterface::parseElementsEntity(myTestProblem);

    // act
    cFemParserInterface::parseNodeLoadsEntity(myTestProblem);
  }

  //! @brief Function to assert - NodalLoads
  //! @author Harikrishnan Sreekumar
  //! @date 11.05.2021
  void TEST_NODALLOADS_ASSERT() {
    // assert
    ASSERT_EQ(mEstimatedNumberOfEntities.NodalForces, 1);
    ASSERT_EQ(myTestProblem.getMesh()->getNumberOfNodalForces(), 1);
    ASSERT_EQ(mEstimatedNumberOfEntities.LoadedNodes, 1);
  }

  //! @brief Function to act - NodeConstraints
  //! @author Harikrishnan Sreekumar
  //! @date 11.05.2021
  void TEST_NODECONSTRAINTS_ACT() {
    // extendended arrange
    cFemParserInterface::parseNodeEntity(myTestProblem);

    // act
    cFemParserInterface::parseNodeConstraintsEntity(myTestProblem);
  }

  //! @brief Function to assert - NodeConstraints
  //! @author Harikrishnan Sreekumar
  //! @date 11.05.2021
  void TEST_NODECONSTRAINTS_ASSERT() {
    // assert
    ASSERT_EQ(mEstimatedNumberOfEntities.NodeBCs, 4);
    ASSERT_EQ(mEstimatedNumberOfEntities.FixedNodes, 4);
  }

 private:
  /// Object for the problem
  cProblem myTestProblem;
};

//! @brief GTEST: Test to ensusre proper reading of Frequency-Anaylsis HDF5
//! group
//! @author Harikrishnan Sreekumar
//! @date 31.10.2020
TEST_F(cTestFemParserHDF5,
       whenHdf5InputStateAreFineExpectProperFrequencyInformation) {
  TEST_ARRANGE(HDF5_PARSER_BASIC_FILE);
  TEST_FREQUENCYANALYSIS_ACT();
  TEST_FREQUENCYANALYSIS_ASSERT();
}

//! @brief GTEST: Test to ensusre proper reading of nodes HDF5 group
//! @author Harikrishnan Sreekumar
//! @date 31.10.2020
TEST_F(cTestFemParserHDF5,
       whenHdf5InputStateAreFineExpectProperNodesInformation) {
  TEST_ARRANGE(HDF5_PARSER_BASIC_FILE);
  TEST_NODES_ACT();
  TEST_NODES_ASSERT();
}

//! @brief GTEST: Test to ensusre proper reading of materials HDF5 group
//! @author Harikrishnan Sreekumar
//! @date 26.04.2021
TEST_F(cTestFemParserHDF5,
       whenHdf5InputStateAreFineExpectProperMaterialsInformation) {
  TEST_ARRANGE(HDF5_PARSER_BASIC_FILE);
  TEST_MATERIALS_ACT();
  TEST_MATERIALS_ASSERT();
}

//! @brief GTEST: Test to ensusre proper reading of elements HDF5 group
//! @author Harikrishnan Sreekumar
//! @date 31.10.2020
TEST_F(cTestFemParserHDF5,
       whenHdf5InputStateAreFineExpectProperElementsInformation) {
  TEST_ARRANGE(HDF5_PARSER_BASIC_FILE);
  TEST_ELEMENTS_ACT();
  TEST_ELEMENTS_ASSERT();
}

//! @brief GTEST: Test to ensusre proper reading of Nodal Loads HDF5 group
//! @author Harikrishnan Sreekumar
//! @date 11.05.2021
TEST_F(cTestFemParserHDF5,
       whenHdf5InputStateAreFineExpectProperNodalLoadsInformation) {
  TEST_ARRANGE(HDF5_PARSER_BASIC_FILE);
  TEST_NODALLOADS_ACT();
  TEST_NODALLOADS_ASSERT();
}

//! @brief GTEST: Test to ensusre proper reading of NodeConstraints HDF5 group
//! @author Harikrishnan Sreekumar
//! @date 11.05.2021
TEST_F(cTestFemParserHDF5,
       whenHdf5InputStateAreFineExpectProperNodeConstraintsInformation) {
  TEST_ARRANGE(HDF5_PARSER_NODECONSTRAINTS_FILE);
  TEST_NODECONSTRAINTS_ACT();
  TEST_NODECONSTRAINTS_ASSERT();
}
