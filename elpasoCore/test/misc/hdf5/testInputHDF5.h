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

#include "../../../source/misc/hdf5/inputh5.h"
#include "../../resources/properties/testproperties.h"

#ifdef HAVE_GTEST
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#endif  // HAVE_GTEST

//! @brief Testing frame for cInputSingletonH5
//! @author Harikrishnan Sreekumar
//! @date 14.05.2021

class cTestInputSingletonH5 : protected cInputSingletonH5,
                              public testing::Test {
 public:
  //! @brief GTEST Setup
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void SetUp() override {
    // do nothing
  }

  //! @brief GTEST Teardown
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void TearDown() override {
    /// close HDF5 File
#ifdef HAVE_HDF5
    cInputSingletonH5::getInstance()->closeContainer();
#endif
  }

  //! @brief Function to arrange
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void TEST_ARRANGE(std::string _filename) {
    // Arrange
    /// test hdf5 file
    std::string elpasoTestResourceDir = ELPASO_TEST_RESOURCE_DIR;
    std::string testFileHDF5 = elpasoTestResourceDir + "unitHDF5/" + _filename;

    /// open HDF5 File
#ifdef HAVE_HDF5
    cInputSingletonH5::getInstance()->setGlobalFileName(
        testFileHDF5);  // Set file name
    cInputSingletonH5::getInstance()->openContainer(ELPASO_H5_READONLY);
#endif
  }

  //! @brief Function to act and assert integer vector read from dataset
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void TEST_ACT_ARRANGE_INTEGERVECTOR_DATASET() {
    std::vector<int> intvector;
    cInputSingletonH5::getInstance()->readIntegerVector(intvector, "/Nodesets",
                                                        "", "/vecNodeset1");

    ASSERT_EQ(intvector.size(), 1);
    ASSERT_EQ(intvector[0], 1);
  }

  //! @brief Function to act and assert integer vector read from attribute
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void TEST_ACT_ARRANGE_INTEGERVECTOR_ATTRIBUTE() {
    std::vector<int> intvector;
    cInputSingletonH5::getInstance()->readIntegerVectorFromAttributeDataset(
        intvector, "/Vector", "/IntegerVector_Attribute", "NS");

    ASSERT_EQ(intvector.size(), 9);
    ASSERT_EQ(intvector[0], 76);
    ASSERT_EQ(intvector[1], 75);
    ASSERT_EQ(intvector[2], 90);
    ASSERT_EQ(intvector[3], 89);
  }

  //! @brief Function to act and assert integer vector read from attribute
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void TEST_ACT_ARRANGE_DENSEDOUBLE_MATRIX() {
    int num_row, num_col;
    std::vector<double> doublevector;
    cInputSingletonH5::getInstance()->readDenseDoubleMatrix(
        doublevector, num_row, num_col, "/Matrices", "", "/DoubleMatrix");

    ASSERT_EQ(doublevector.size(), 4);
    ASSERT_EQ(num_row, 2);
    ASSERT_EQ(num_col, 2);
    ASSERT_EQ(doublevector[0], 10);
    ASSERT_EQ(doublevector[1], 0.01);
    ASSERT_EQ(doublevector[2], 100);
    ASSERT_EQ(doublevector[3], 0.001);
  }

  //! @brief Function to act and assert double vector read from dataset
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void TEST_ACT_ARRANGE_DOUBLEVECTOR_DATASET() {
    std::vector<double> doublevector;
    cInputSingletonH5::getInstance()->readDoubleVector(doublevector, "/Vector",
                                                       "", "/DoubleVector64");

    ASSERT_EQ(doublevector.size(), 4);
    ASSERT_EQ(doublevector[0], 4.5);
    ASSERT_EQ(doublevector[1], 0.0001);
    ASSERT_EQ(doublevector[2], 500.01);
    ASSERT_EQ(doublevector[3], 4.2e-4);
  }

  //! @brief Function to act and assert proper reading of node compound data
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void TEST_ACT_ARRANGE_NODEDATA_COMPOUND() {
    std::vector<int> ids;
    std::vector<double> coord;
    cInputSingletonH5::getInstance()->readNodeCompoundData(ids, coord, "/Nodes",
                                                           "", "/mtxFemNodes");

    ASSERT_EQ(ids.size(), 6);
    ASSERT_EQ(coord.size(), 6 * 3);

    ASSERT_EQ(ids[0], 1);
    ASSERT_EQ(ids[1], 2);
    ASSERT_EQ(ids[2], 4);
    ASSERT_EQ(ids[3], 5);

    ASSERT_EQ(coord[0], 0.8);
    ASSERT_EQ(coord[1], 0.4);
    ASSERT_EQ(coord[2], 0.0);
  }

  //! @brief Function to act and assert proper reading of complex compound data
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void TEST_ACT_ARRANGE_COMPLEXDATA_COMPOUND() {
    std::vector<PetscComplex> complexvector;
    cInputSingletonH5::getInstance()->readComplexVector(
        complexvector, "/Compound", "", "/vecFemStep1");

    ASSERT_EQ(complexvector.size(), 918);
    ASSERT_EQ(complexvector[3].real(), 2.0054695597675435E-7);
    ASSERT_EQ(complexvector[3].imag(), 0.);
  }

  //! @brief Function to act and assert proper reading of integer attribute from
  //! group
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void TEST_ACT_ARRANGE_INTEGERATTRIBUTE_GROUP() {
    int attribute =
        cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(
            "steps", "/Analysis");

    ASSERT_EQ(attribute, 1);
  }

  //! @brief Function to act and assert proper reading of double attribute from
  //! group
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void TEST_ACT_ARRANGE_DOUBLEATTRIBUTE_GROUP() {
    double attribute =
        cInputSingletonH5::getInstance()->readDoubleAttributeFromGroup(
            "start", "/Analysis");

    ASSERT_EQ(attribute, 10);
  }

  //! @brief Function to act and assert proper reading of double attribute from
  //! dataset
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void TEST_ACT_ARRANGE_DOUBLEATTRIBUTE_DATASET() {
    double attribute =
        cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(
            "valu3", "/NodeConstraints", "/nodeConstraint1_1");

    ASSERT_EQ(attribute, 0);
  }

  //! @brief Function to act and assert proper reading of non-existing double
  //! attribute from dataset
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void TEST_ACT_ARRANGE_NOEXIST_DOUBLEATTRIBUTE_DATASET() {
    double attribute =
        cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(
            "valp", "/NodeConstraints", "/nodeConstraint1_1");

    ASSERT_EQ(attribute, 0);
  }

  //! @brief Function to act and assert proper reading of string attribute from
  //! group
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void TEST_ACT_ARRANGE_STRINGATTRIBUTE_GROUP() {
    std::string attribute =
        cInputSingletonH5::getInstance()->readStringAttributeFromGroup(
            "type", "/Analysis");

    ASSERT_EQ(attribute, "frequency-basic");
  }

  //! @brief Function to act and assert proper reading of string attribute from
  //! dataset
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void TEST_ACT_ARRANGE_STRINGATTRIBUTE_DATASET() {
    std::string attribute =
        cInputSingletonH5::getInstance()->readStringAttributeFromDataset(
            "MaterialType", "/Materials", "/material1");

    ASSERT_EQ(attribute, "STR_LIN_ELA_ISO_DIR");
  }

  //! @brief Function to act and assert proper reading of integer attribute from
  //! dataset
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void TEST_ACT_ARRANGE_INTEGERATTRIBUTE_DATASET() {
    int attribute =
        cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(
            "MaterialId", "/Elements", "/mtxFemElemGroup1");

    ASSERT_EQ(attribute, 1);
  }

  //! @brief Function to act and assert proper reading of number of members in
  //! group
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void TEST_ACT_ARRANGE_NUMBERMEMBERS() {
    int members = cInputSingletonH5::getInstance()->getNumberOfMembersInGroup(
        "/ElemLoads");

    ASSERT_EQ(members, 0);
  }

  //! @brief Function to act and assert proper reading of member name
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void TEST_ACT_ARRANGE_MEMBERNAME() {
    std::string member =
        cInputSingletonH5::getInstance()->getMemberName("/NodeLoads", 0);

    ASSERT_EQ(member, "/mtxFemNodeLoad1_1");
  }

  //! @brief Function to act and assert proper reading of complex dense matrix
  //! @author Harikrishnan Sreekumar
  //! @date 16.05.2021
  void TEST_ACT_ARRANGE_COMPLEXDENSEMAT() {
    std::vector<PetscComplex> complex_mat;
    int num_row, num_col;
    cInputSingletonH5::getInstance()->readDenseComplexMatrix(
        complex_mat, num_row, num_col, "/Compound", "", "/mtxStiffnessRom");

    ASSERT_EQ(num_row * num_col, 2025);
    ASSERT_EQ(complex_mat.size(), 2025);
    ASSERT_NEAR(complex_mat[0].real(), 36621.04312976473, 1e-4);
    ASSERT_EQ(complex_mat[0].imag(), 0);
    ASSERT_NEAR(complex_mat[1].real(), -1937.9915546127222, 1e-4);
  }
};

//! @brief GTEST: Test to ensusre proper reading of integer vector
//! @author Harikrishnan Sreekumar
//! @date 14.05.2021
TEST_F(cTestInputSingletonH5,
       whenHd5InputStatesAreFineExpectProperReadingOfIntegerVectorFromDataSet) {
  TEST_ARRANGE(HDF5_PARSER_BASIC_FILE);
  TEST_ACT_ARRANGE_INTEGERVECTOR_DATASET();
}

//! @brief GTEST: Test to ensusre proper reading of integer vector from
//! attribute
//! @author Harikrishnan Sreekumar
//! @date 14.05.2021
TEST_F(
    cTestInputSingletonH5,
    whenHd5InputStatesAreFineExpectProperReadingOfIntegerVectorFromAttribute) {
  TEST_ARRANGE(HDF5_PARSER_TEST_FILE);
  TEST_ACT_ARRANGE_INTEGERVECTOR_ATTRIBUTE();
}

//! @brief GTEST: Test to ensusre proper reading of dense double matrix
//! @author Harikrishnan Sreekumar
//! @date 14.05.2021
TEST_F(cTestInputSingletonH5,
       whenHd5InputStatesAreFineExpectProperReadingOfDenseDoubleMatrix) {
  TEST_ARRANGE(HDF5_PARSER_TEST_FILE);
  TEST_ACT_ARRANGE_DENSEDOUBLE_MATRIX();
}

//! @brief GTEST: Test to ensusre proper reading of double vector
//! @author Harikrishnan Sreekumar
//! @date 14.05.2021
TEST_F(cTestInputSingletonH5,
       whenHd5InputStatesAreFineExpectProperReadingOfDoubleVectorFromDataSet) {
  TEST_ARRANGE(HDF5_PARSER_TEST_FILE);
  TEST_ACT_ARRANGE_DOUBLEVECTOR_DATASET();
}

//! @brief GTEST: Test to ensusre proper reading of complex compound data
//! @author Harikrishnan Sreekumar
//! @date 14.05.2021
TEST_F(cTestInputSingletonH5,
       whenHd5InputStatesAreFineExpectProperReadingOfCompoundComplexData) {
  TEST_ARRANGE(HDF5_PARSER_TEST_FILE);
  TEST_ACT_ARRANGE_COMPLEXDATA_COMPOUND();
}

//! @brief GTEST: Test to ensusre proper reading of node compound data
//! @author Harikrishnan Sreekumar
//! @date 14.05.2021
TEST_F(cTestInputSingletonH5,
       whenHd5InputStatesAreFineExpectProperReadingOfCompoundNodeData) {
  TEST_ARRANGE(HDF5_PARSER_BASIC_FILE);
  TEST_ACT_ARRANGE_NODEDATA_COMPOUND();
}

//! @brief GTEST: Test to ensusre proper reading of integer attribute from group
//! @author Harikrishnan Sreekumar
//! @date 14.05.2021
TEST_F(
    cTestInputSingletonH5,
    whenHd5InputStatesAreFineExpectProperReadingOfIntegerAttributeFromGroup) {
  TEST_ARRANGE(HDF5_PARSER_BASIC_FILE);
  TEST_ACT_ARRANGE_INTEGERATTRIBUTE_GROUP();
}

//! @brief GTEST: Test to ensusre proper reading of double attribute from group
//! @author Harikrishnan Sreekumar
//! @date 14.05.2021
TEST_F(cTestInputSingletonH5,
       whenHd5InputStatesAreFineExpectProperReadingOfDoubleAttributeFromGroup) {
  TEST_ARRANGE(HDF5_PARSER_BASIC_FILE);
  TEST_ACT_ARRANGE_DOUBLEATTRIBUTE_GROUP();
}

//! @brief GTEST: Test to ensusre proper reading of double attribute from
//! dataset
//! @author Harikrishnan Sreekumar
//! @date 14.05.2021
TEST_F(
    cTestInputSingletonH5,
    whenHd5InputStatesAreFineExpectProperReadingOfDoubleAttributeFromDataset) {
  TEST_ARRANGE(HDF5_PARSER_NODECONSTRAINTS_FILE);
  TEST_ACT_ARRANGE_DOUBLEATTRIBUTE_DATASET();
}

//! @brief GTEST: Test to ensusre no error is thrown
//! @author Harikrishnan Sreekumar
//! @date 14.05.2021
TEST_F(
    cTestInputSingletonH5,
    whenHd5InputStatesAreNotFineExpectReadingOfDoubleAttributeFromDatasetThrowsNoError) {
  TEST_ARRANGE(HDF5_PARSER_NODECONSTRAINTS_FILE);
  TEST_ACT_ARRANGE_NOEXIST_DOUBLEATTRIBUTE_DATASET();
}

//! @brief GTEST: Test to ensusre proper reading of string attribute from group
//! @author Harikrishnan Sreekumar
//! @date 14.05.2021
TEST_F(cTestInputSingletonH5,
       whenHd5InputStatesAreFineExpectProperReadingOfStringAttributeFromGroup) {
  TEST_ARRANGE(HDF5_PARSER_BASIC_FILE);
  TEST_ACT_ARRANGE_STRINGATTRIBUTE_GROUP();
}

//! @brief GTEST: Test to ensusre proper reading of string attribute from
//! dataset
//! @author Harikrishnan Sreekumar
//! @date 14.05.2021
TEST_F(
    cTestInputSingletonH5,
    whenHd5InputStatesAreFineExpectProperReadingOfStringAttributeFromDataset) {
  TEST_ARRANGE(HDF5_PARSER_BASIC_FILE);
  TEST_ACT_ARRANGE_STRINGATTRIBUTE_DATASET();
}

//! @brief GTEST: Test to ensusre proper reading of integer attribute from
//! dataset
//! @author Harikrishnan Sreekumar
//! @date 14.05.2021
TEST_F(
    cTestInputSingletonH5,
    whenHd5InputStatesAreFineExpectProperReadingOfIntegerAttributeFromDataset) {
  TEST_ARRANGE(HDF5_PARSER_NODECONSTRAINTS_FILE);
  TEST_ACT_ARRANGE_INTEGERATTRIBUTE_DATASET();
}

//! @brief GTEST: Test to ensusre proper reading of number of members in group
//! @author Harikrishnan Sreekumar
//! @date 14.05.2021
TEST_F(cTestInputSingletonH5,
       whenHd5InputStatesAreFineExpectProperReadingOfNumberOfMembersInGroup) {
  TEST_ARRANGE(HDF5_PARSER_BASIC_FILE);
  TEST_ACT_ARRANGE_NUMBERMEMBERS();
}

//! @brief GTEST: Test to ensusre proper reading of member name by index
//! @author Harikrishnan Sreekumar
//! @date 14.05.2021
TEST_F(cTestInputSingletonH5,
       whenHd5InputStatesAreFineExpectProperReadingMemberName) {
  TEST_ARRANGE(HDF5_PARSER_BASIC_FILE);
  TEST_ACT_ARRANGE_MEMBERNAME();
}

//! @brief GTEST: Test to ensusre proper reading of complex dense matrix
//! @author Harikrishnan Sreekumar
//! @date 16.05.2021
TEST_F(cTestInputSingletonH5,
       whenHd5InputStatesAreFineExpectProperReadingComplexDenseMatrix) {
  TEST_ARRANGE(HDF5_PARSER_TEST_FILE);
  TEST_ACT_ARRANGE_COMPLEXDENSEMAT();
}