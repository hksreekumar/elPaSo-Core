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

#include "../../../source/misc/hdf5/outputh5.h"
#include "../../../source/misc/hdf5/readerh5.h"
#include "../../resources/properties/testproperties.h"

#ifdef HAVE_GTEST
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#endif  // HAVE_GTEST

//! @brief Testing frame for cOutputH5
//! @author Harikrishnan Sreekumar
//! @date 14.05.2021
class cTestOutputH5 : protected cOutputH5, public testing::Test {
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
  void TearDown() override {}

  //! @brief Function to arrange
  //! @author Harikrishnan Sreekumar
  //! @date 14.05.2021
  void TEST_ARRANGE(std::string _filename) {
    // Arrange
    /// test hdf5 file
    std::string elpasoTestResourceDir = ELPASO_TEST_RESOURCE_DIR;
    std::string testFileHDF5 = elpasoTestResourceDir + "unitHDF5/" + _filename;

    /// open HDF5 File
    myOutputWriter.setGlobalFileName(testFileHDF5);
    myOutputWriter.openContainer(ELPASO_H5_READWRITE_FORCE);
  }

  //! @brief Function to act - proper group creation
  //! @author Harikrishnan Sreekumar
  //! @date 15.05.2021
  void TEST_ACT_CREATEGROUPS() {
    myOutputWriter.createGroup("/ACT_GROUP_MAIN");
    myOutputWriter.createSecoGroup("/ACT_GROUP_MAIN", "/ACT_GROUP_SECO");
    myOutputWriter.closeContainer();
  }

  //! @brief Function to assert - proper group creation
  //! @author Harikrishnan Sreekumar
  //! @date 15.05.2021
  void TEST_ASSERT_CREATEGROUPS() {
    cReaderH5 myReader;
    myReader.setGlobalFileName(myOutputWriter.getGlobalFileName());
    myReader.openContainer(ELPASO_H5_READONLY);

    std::string seco_group_name = myReader.getMemberName("/ACT_GROUP_MAIN", 0);
    ASSERT_EQ(seco_group_name, "/ACT_GROUP_SECO");

    myReader.closeContainer();
  }

  //! @brief Function to act - proper int vector creation
  //! @author Harikrishnan Sreekumar
  //! @date 15.05.2021
  void TEST_ACT_INTEGERVECTOR() {
    std::vector<int> exportvec;
    exportvec.reserve(100);
    for (size_t i = 0; i < 100; i++) exportvec.push_back(i);

    myOutputWriter.createGroup("/ACT_GROUP_MAIN");
    myOutputWriter.appendIntegerVector(exportvec, "/ACT_GROUP_MAIN", "",
                                       "/testvec");
    myOutputWriter.closeContainer();
  }

  //! @brief Function to assert - proper int vector creation
  //! @author Harikrishnan Sreekumar
  //! @date 15.05.2021
  void TEST_ASSERT_INTEGERVECTOR() {
    cReaderH5 myReader;
    myReader.setGlobalFileName(myOutputWriter.getGlobalFileName());
    myReader.openContainer(ELPASO_H5_READONLY);

    std::vector<int> readvector;
    myReader.readIntegerVector(readvector, "/ACT_GROUP_MAIN", "", "/testvec");

    ASSERT_EQ(readvector.size(), 100);
    for (size_t i = 0; i < 100; i++) ASSERT_EQ(readvector[i], i);

    myReader.closeContainer();
  }

#ifdef PETSC_USE_COMPLEX

  //! @brief Function to act - proper petsc complex vector creation
  //! @author Harikrishnan Sreekumar
  //! @date 15.05.2021
  void TEST_ACT_PETSCCOMPLEXVECTOR() {
    std::vector<PetscComplex> exportvec;
    exportvec.reserve(100);
    for (size_t i = 0; i < 100; i++)
      exportvec.push_back((double)i + (double)i * PETSC_i);

    myOutputWriter.createGroup("/ACT_GROUP_MAIN");
    myOutputWriter.appendComplexVector(exportvec, "/ACT_GROUP_MAIN", "",
                                       "/testvec");
    myOutputWriter.closeContainer();
  }

  //! @brief Function to act - proper petsc complex vector creation - with
  //! hyperslab
  //! @author Harikrishnan Sreekumar
  //! @date 22.04.2022
  void TEST_ACT_PETSCCOMPLEXVECTORHYPERSLABS() {
    std::vector<PetscComplex> exportvec;
    exportvec.reserve(100);
    for (size_t i = 0; i < 100; i++)
      exportvec.push_back((double)i + (double)i * PETSC_i);

    myOutputWriter.createGroup("/ACT_GROUP_MAIN");
    myOutputWriter.appendComplexVector(exportvec, "/ACT_GROUP_MAIN", "",
                                       "/testvec", 0, 100);
    myOutputWriter.closeContainer();
  }

  //! @brief Function to assert - proper petsc complex vector creation
  //! @author Harikrishnan Sreekumar
  //! @date 15.05.2021
  void TEST_ASSERT_PETSCCOMPLEXVECTOR() {
    cReaderH5 myReader;
    myReader.setGlobalFileName(myOutputWriter.getGlobalFileName());
    myReader.openContainer(ELPASO_H5_READONLY);

    std::vector<PetscComplex> readvector;
    myReader.readComplexVector(readvector, "/ACT_GROUP_MAIN", "", "/testvec");

    ASSERT_EQ(readvector.size(), 100);
    for (size_t i = 0; i < 100; i++) {
      ASSERT_EQ(readvector[i].real(), i);
      ASSERT_EQ(readvector[i].imag(), i);
    }

    myReader.closeContainer();
  }
#endif

#ifdef PETSC_USE_COMPLEX

  //! @brief Function to act - proper petsc complex matrix creation
  //! @author Harikrishnan Sreekumar
  //! @date 15.05.2021
  void TEST_ACT_PETSCCOMPLEXMATRIX() {
    std::vector<PetscComplex> exportmat;
    exportmat.reserve(100);
    for (size_t i = 0; i < 100; i++)
      exportmat.push_back((double)i + (double)i * PETSC_i);

    myOutputWriter.createGroup("/ACT_GROUP_MAIN");
    myOutputWriter.appendComplexMatrix(exportmat, 10, 10, "/ACT_GROUP_MAIN", "",
                                       "/testvec", PETSC_TRUE);
    myOutputWriter.closeContainer();
  }

  //! @brief Function to assert - proper petsc complex matrix creation
  //! @author Harikrishnan Sreekumar
  //! @date 15.05.2021
  void TEST_ASSERT_PETSCCOMPLEXMATRIX() {
    cReaderH5 myReader;
    myReader.setGlobalFileName(myOutputWriter.getGlobalFileName());
    myReader.openContainer(ELPASO_H5_READONLY);

    std::vector<PetscComplex> readmatrix;
    int num_row, num_col;
    myReader.readDenseComplexMatrix(readmatrix, num_row, num_col,
                                    "/ACT_GROUP_MAIN", "", "/testvec");

    ASSERT_EQ(num_row * num_col, 100);
    for (size_t i = 0; i < 100; i++) {
      ASSERT_EQ(readmatrix[i].real(), i);
      ASSERT_EQ(readmatrix[i].imag(), i);
    }

    myReader.closeContainer();
  }
#endif

  //! @brief Function to act - proper string attribute creation
  //! @author Harikrishnan Sreekumar
  //! @date 15.05.2021
  void TEST_ACT_STRINGATTRIBUTE() {
    myOutputWriter.createGroup("/ACT_GROUP_MAIN");
    myOutputWriter.writeStringAttributeToGroup("elPaSo", "/ACT_GROUP_MAIN",
                                               "2021");
    myOutputWriter.closeContainer();
  }

  //! @brief Function to assert - proper string attribute creation
  //! @author Harikrishnan Sreekumar
  //! @date 15.05.2021
  void TEST_ASSERT_STRINGATTRIBUTE() {
    cReaderH5 myReader;
    myReader.setGlobalFileName(myOutputWriter.getGlobalFileName());
    myReader.openContainer(ELPASO_H5_READONLY);

    std::string readstring =
        myReader.readStringAttributeFromGroup("elPaSo", "/ACT_GROUP_MAIN");

    ASSERT_EQ(readstring, "2021");

    myReader.closeContainer();
  }

  //! @brief Function to act - proper integer attribute creation
  //! @author Harikrishnan Sreekumar
  //! @date 25.05.2021
  void TEST_ACT_INTEGERATTRIBUTE() {
    myOutputWriter.createGroup("/ACT_GROUP_MAIN");
    myOutputWriter.writeIntegerAttributeToGroup("INTEGER", "/ACT_GROUP_MAIN",
                                                2021);
    myOutputWriter.closeContainer();
  }

  //! @brief Function to assert - proper integer attribute creation
  //! @author Harikrishnan Sreekumar
  //! @date 25.05.2021
  void TEST_ASSERT_INTEGERATTRIBUTE() {
    cReaderH5 myReader;
    myReader.setGlobalFileName(myOutputWriter.getGlobalFileName());
    myReader.openContainer(ELPASO_H5_READONLY);

    int readint =
        myReader.readIntegerAttributeFromGroup("INTEGER", "/ACT_GROUP_MAIN");

    ASSERT_EQ(readint, 2021);

    myReader.closeContainer();
  }

  //! @brief Function to act - proper petsc real matrix creation
  //! @author Harikrishnan Sreekumar
  //! @date 15.05.2021
  void TEST_ACT_PETSCREALMATRIX() {
    std::vector<PetscReal> exportmat;
    exportmat.reserve(100);
    for (size_t i = 0; i < 100; i++) exportmat.push_back((double)i);

    myOutputWriter.createGroup("/ACT_GROUP_MAIN");
    myOutputWriter.appendRealDenseMatrix(exportmat, 10, 10, "/ACT_GROUP_MAIN",
                                         "", "/testmat");
    myOutputWriter.closeContainer();
  }

  //! @brief Function to assert - proper petsc real matrix creation
  //! @author Harikrishnan Sreekumar
  //! @date 15.05.2021
  void TEST_ASSERT_PETSCREALMATRIX() {
    cReaderH5 myReader;
    myReader.setGlobalFileName(myOutputWriter.getGlobalFileName());
    myReader.openContainer(ELPASO_H5_READONLY);

    std::vector<PetscReal> readmatrix;
    int num_row, num_col;
    myReader.readDenseDoubleMatrix(readmatrix, num_row, num_col,
                                   "/ACT_GROUP_MAIN", "", "/testmat");

    ASSERT_EQ(num_row * num_col, 100);
    for (size_t i = 0; i < 100; i++) {
      ASSERT_EQ(readmatrix[i], i);
    }

    myReader.closeContainer();
  }

  //! @brief Function to act - proper petsc integer matrix creation
  //! @author Harikrishnan Sreekumar
  //! @date 14.12.2021
  void TEST_ACT_PETSCINTEGERMATRIX() {
    std::vector<PetscInt> exportmat;
    exportmat.reserve(100);
    for (size_t i = 0; i < 100; i++) exportmat.push_back((int)i);

    myOutputWriter.createGroup("/ACT_GROUP_MAIN");
    myOutputWriter.appendIntegerDenseMatrix(exportmat, 10, 10,
                                            "/ACT_GROUP_MAIN", "", "/testmat");
    myOutputWriter.closeContainer();
  }

  //! @brief Function to assert - proper petsc integer matrix creation
  //! @author Harikrishnan Sreekumar
  //! @date 14.12.2021
  void TEST_ASSERT_PETSCINTEGERMATRIX() {
    cReaderH5 myReader;
    myReader.setGlobalFileName(myOutputWriter.getGlobalFileName());
    myReader.openContainer(ELPASO_H5_READONLY);

    std::vector<PetscInt> readmatrix;
    int num_row, num_col;
    myReader.readDenseIntegerMatrix(readmatrix, num_row, num_col,
                                    "/ACT_GROUP_MAIN", "", "/testmat");

    ASSERT_EQ(num_row * num_col, 100);
    for (size_t i = 0; i < 100; i++) {
      ASSERT_EQ(readmatrix[i], i);
    }

    myReader.closeContainer();
  }

  //! @brief Function to act and assert - proper petsc integer vector attribute
  //! write
  //! @author Harikrishnan Sreekumar
  //! @date 28.01.2022
  void TEST_ACT_ASSERT_PETSCINTEGERVECTORATTRIBUTE() {
    // act
    std::vector<int> exportvec;
    exportvec.reserve(10);
    for (size_t i = 1; i <= 10; i++) exportvec.push_back(i);

    myOutputWriter.createGroup("/ACT_GROUP2_MAIN");
    myOutputWriter.writeIntegerVectorAttributeToGroup(
        exportvec, "/ACT_GROUP2_MAIN", "att_vec");
    myOutputWriter.closeContainer();

    // assert
    cReaderH5 myReader;
    myReader.setGlobalFileName(myOutputWriter.getGlobalFileName());
    myReader.openContainer(ELPASO_H5_READONLY);

    std::vector<int> readvector;
    myReader.readIntegerVectorFromAttributeGroup(readvector, "/ACT_GROUP2_MAIN",
                                                 "att_vec");

    ASSERT_EQ(readvector.size(), 10);
    for (size_t i = 1; i <= 10; i++) ASSERT_EQ(readvector[i - 1], i);

    myReader.closeContainer();
  }

 private:
  cOutputH5 myOutputWriter;
};

//! @brief GTEST: Test to ensure proper creation of group
//! @author Harikrishnan Sreekumar
//! @date 16.05.2021
TEST_F(cTestOutputH5, checkProperOutputOfGroupsInHDF5) {
  TEST_ARRANGE(HDF5_PARSER_OUTPUT_EXAMPLE_FILE);
  TEST_ACT_CREATEGROUPS();
  TEST_ASSERT_CREATEGROUPS();
}

//! @brief GTEST: Test to ensure proper creation of integer vector
//! @author Harikrishnan Sreekumar
//! @date 16.05.2021
TEST_F(cTestOutputH5, checkProperOutputOfIntegerVectorInHDF5) {
  TEST_ARRANGE(HDF5_PARSER_OUTPUT_EXAMPLE_FILE);
  TEST_ACT_INTEGERVECTOR();
  TEST_ASSERT_INTEGERVECTOR();
}

#ifdef PETSC_USE_COMPLEX

//! @brief GTEST: Test to ensure proper creation of petsc complex vector
//! @author Harikrishnan Sreekumar
//! @date 16.05.2021
TEST_F(cTestOutputH5, checkProperOutputOfPetscComplexVectorInHDF5) {
  TEST_ARRANGE(HDF5_PARSER_OUTPUT_EXAMPLE_FILE);
  TEST_ACT_PETSCCOMPLEXVECTOR();
  TEST_ASSERT_PETSCCOMPLEXVECTOR();
}
#endif

#ifdef PETSC_USE_COMPLEX

//! @brief GTEST: Test to ensusre proper creation of petsc complex vector - with
//! hyperslabs
//! @author Harikrishnan Sreekumar
//! @date 22.04.2022
TEST_F(cTestOutputH5,
       checkProperOutputOfPetscComplexVectorForHyperSlabsInHDF5) {
  TEST_ARRANGE(HDF5_PARSER_OUTPUT_EXAMPLE_FILE);
  TEST_ACT_PETSCCOMPLEXVECTORHYPERSLABS();
  TEST_ASSERT_PETSCCOMPLEXVECTOR();
}
#endif

#ifdef PETSC_USE_COMPLEX

//! @brief GTEST: Test to ensusre proper creation of petsc complex matrix
//! @author Harikrishnan Sreekumar
//! @date 16.05.2021
TEST_F(cTestOutputH5, checkProperOutputOfPetscComplexMatrixInHDF5) {
  TEST_ARRANGE(HDF5_PARSER_OUTPUT_EXAMPLE_FILE);
  TEST_ACT_PETSCCOMPLEXMATRIX();
  TEST_ASSERT_PETSCCOMPLEXMATRIX();
}
#endif

//! @brief GTEST: Test to ensure proper creation of string attribute
//! @author Harikrishnan Sreekumar
//! @date 16.05.2021
TEST_F(cTestOutputH5, checkProperOutputOfStringAttributeInHDF5) {
  TEST_ARRANGE(HDF5_PARSER_OUTPUT_EXAMPLE_FILE);
  TEST_ACT_STRINGATTRIBUTE();
  TEST_ASSERT_STRINGATTRIBUTE();
}

//! @brief GTEST: Test to ensure proper creation of integer attribute
//! @author Harikrishnan Sreekumar
//! @date 25.05.2021
TEST_F(cTestOutputH5, checkProperOutputOfIntegerAttributeInHDF5) {
  TEST_ARRANGE(HDF5_PARSER_OUTPUT_EXAMPLE_FILE);
  TEST_ACT_INTEGERATTRIBUTE();
  TEST_ASSERT_INTEGERATTRIBUTE();
}

//! @brief GTEST: Test to ensure proper creation of petsc real matrix
//! @author Harikrishnan Sreekumar
//! @date 16.05.2021
TEST_F(cTestOutputH5, checkProperOutputOfPetscRealMatrixInHDF5) {
  TEST_ARRANGE(HDF5_PARSER_OUTPUT_EXAMPLE_FILE);
  TEST_ACT_PETSCREALMATRIX();
  TEST_ASSERT_PETSCREALMATRIX();
}

//! @brief GTEST: Test to ensusre proper creation of petsc integer matrix
//! @author Harikrishnan Sreekumar
//! @date 14.12.2021
TEST_F(cTestOutputH5, checkProperOutputOfPetscIntegerMatrixInHDF5) {
  TEST_ARRANGE(HDF5_PARSER_OUTPUT_EXAMPLE_FILE);
  TEST_ACT_PETSCINTEGERMATRIX();
  TEST_ASSERT_PETSCINTEGERMATRIX();
}

//! @brief GTEST: Test to ensure proper writing of petsc integer vector as
//! attribute
//! @author Harikrishnan Sreekumar
//! @date 28.01.2022
TEST_F(cTestOutputH5,
       checkProperOutputOfPetscIntegerVectorAttributeGroupInHDF5) {
  TEST_ARRANGE(HDF5_PARSER_OUTPUT_EXAMPLE_FILE);
  TEST_ACT_ASSERT_PETSCINTEGERVECTORATTRIBUTE();
}