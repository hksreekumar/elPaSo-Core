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

#include "../../../../../../source/misc/parser/femparserhdf5.h"
#include "../../../../../resources/properties/testproperties.h"

#ifdef HAVE_GTEST
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#endif  // HAVE_GTEST

//! @brief Mock class for cMaterialStructureVisco
//! @author Akash Doddamane Lingaraja
//! @date 07.01.2022
class cTestMaterialStructureVisco : private cFemParserHDF5,
                                    public testing::Test {
 public:
  //! @brief GTEST Setup
  //! @author Akash Doddamane Lingaraja
  //! @date 07.01.2022
  void SetUp() {
    // do nothing
  }

  //! @brief GTEST Teardown
  //! @author Akash Doddamane Lingaraja
  //! @date 07.01.2022
  void TearDown() {
    // close HDF5 File
#ifdef HAVE_HDF5
    cInputSingletonH5::getInstance()->closeContainer();
#endif
  }

  //! @brief Function to arrange - do parse the HDF5 file
  //! @author Akash Doddamane Lingaraja
  //! @date 07.01.2022
  void TEST_ARRANGE(std::string _filename) {
    std::string filename = _filename;
    std::string elpasoTestResourceDir = ELPASO_TEST_RESOURCE_DIR;
    std::string testFileHDF5 =
        elpasoTestResourceDir +
        "unitHDF5/material/structure/isotrop/linear/viscoelastic/" + filename;

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
  //! @date 07.01.2022
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

  //! @brief Function to assert - do get frequency dependent Young's modulus
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_EOMEGA() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    // std::cout << (mat->getEOmega()) << "\n";
    PetscScalar eVis = 7e10 + 7e7 * PETSC_i;
    ASSERT_EQ((mat->getEOmega()), eVis);
  }

  //! @brief Function to assert - do get Young's modulus
  //! @author Akash Doddamane Lingaraja
  //! @date 07.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_GETE() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    // std::cout << (mat->getE()) << "\n";
    ASSERT_EQ((mat->getE()), 7e10);
  }

  //! @brief Function to assert - do get Eta
  //! @author Akash Doddamane Lingaraja
  //! @date 17.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_GETETA() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    // std::cout << (mat->getEta()) << "\n";
    ASSERT_EQ((mat->getEta()), 0.001);
  }

  //! @brief Function to assert - do get Poisson's ratio
  //! @author Akash Doddamane Lingaraja
  //! @date 07.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_GETNU() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    // std::cout << (mat->getNu()) << "\n";
    ASSERT_EQ((mat->getNu()), 0.3);
  }

  //! @brief Function to assert - do get Visco-elastic type
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_GETTYPE() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    // std::cout << (mat->getType()) << "\n";
    ASSERT_EQ((mat->getType()), 1);
  }

  //! @brief Function to assert - do set Youngs Modulus
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_SETE() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    mat->setE(5e5);
    // std::cout << (mat->getE()) << "\n";
    ASSERT_EQ((mat->getE()), 5e5);
  }

  //! @brief Function to assert - do set Eta
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_SETETA() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    mat->setEta(5);
    // std::cout << (mat->getEta()) << "\n";
    ASSERT_EQ((mat->getEta()), 5);
  }

  //! @brief Function to assert - do set Poisson's ratio
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_SETNU() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    mat->setNu(5);
    // std::cout << (mat->getNu()) << "\n";
    ASSERT_EQ((mat->getNu()), 5);
  }

  //! @brief Function to assert - do set Visco-elastic type
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_SETTYPE() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    // std::cout << (mat->getType()) << "\n";
    mat->setType(1);
    ASSERT_EQ((mat->getType()), 1);
  }

  //! @brief Function to assert - do get name of the file storing the eta values
  //! @author Akash Doddamane Lingaraja
  //! @date 28.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_GETFILENAMEETA() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    // std::cout << (mat->getFilenameEta()) << "\n";
    std::string refVal = "/Parameters/material1_eta";
    ASSERT_EQ(refVal, mat->getFilenameEta());
  }

  //! @brief Function to assert - do get name of the file storing the youngs
  //! modulus values
  //! @author Akash Doddamane Lingaraja
  //! @date 28.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_GETFILENAMEE() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    // std::cout << (mat->getFilenameE()) << "\n";
    std::string refVal = "/Parameters/material1_E";
    ASSERT_EQ(refVal, mat->getFilenameE());
  }

  //! @brief Function to assert - do get Spring constant
  //! @author Akash Doddamane Lingaraja
  //! @date 17.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_SPRINGCONST() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    // std::cout << (mat->getSpringConst()) << "\n";
    PetscScalar eVis = 7e10 + 7e7 * PETSC_i;
    ASSERT_EQ((mat->getSpringConst()), eVis);
  }

  //! @brief Function to assert - do get Thickness
  //! @author Akash Doddamane Lingaraja
  //! @date 07.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_THICKNESS() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    // std::cout << (mat->getT()) << "\n";
    ASSERT_EQ((mat->getT()), 0.003);
  }

  //! @brief Function to assert - do get centroid of the element
  //! @author Akash Doddamane Lingaraja
  //! @date 28.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_CALCMEANCOORDINATE() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);

    cElement* elem1 = myTestProblem.getMesh()->getElement(1);
    // std::cout << mat->calculateMeanCoordinate(elem1, 0) << "	" <<
    // mat->calculateMeanCoordinate(elem1, 1) << "	" <<
    // mat->calculateMeanCoordinate(elem1, 2) << "\n";
    ASSERT_NEAR(0.6, mat->calculateMeanCoordinate(elem1, 0), 1e-6);
    ASSERT_NEAR(0.2, mat->calculateMeanCoordinate(elem1, 1), 1e-6);
    ASSERT_NEAR(0.0, mat->calculateMeanCoordinate(elem1, 2), 1e-6);

    cElement* elem2 = myTestProblem.getMesh()->getElement(2);
    // std::cout << mat->calculateMeanCoordinate(elem2, 0) << "	" <<
    // mat->calculateMeanCoordinate(elem2, 1) << "	" <<
    // mat->calculateMeanCoordinate(elem2, 2) << "\n";
    ASSERT_NEAR(0.2, mat->calculateMeanCoordinate(elem2, 0), 1e-6);
    ASSERT_NEAR(0.2, mat->calculateMeanCoordinate(elem2, 1), 1e-6);
    ASSERT_NEAR(0.0, mat->calculateMeanCoordinate(elem2, 2), 1e-6);
  }

  //! @brief Function to assert - do get distance between two points
  //! @author Akash Doddamane Lingaraja
  //! @date 28.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_CALCDIST() {
    PetscReal x1 = 2.0, y1 = 1.0;
    PetscReal x2 = 3.0, y2 = 4.0;
    PetscReal refVal = 3.1622776601684;

    cMaterialStructureVisco* mat;
    PetscReal val = mat->calculateDistance(x1, y1, x2, y2);
    // std::cout << val << "\n";
    ASSERT_NEAR(val, refVal, 1e-6);
  }

  //! @brief Function to assert - do get shear part of the elasticity matrix of
  //! Mindlin plate
  //! @author Akash Doddamane Lingaraja
  //! @date 17.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_SETUPCS() {
#ifdef PETSC_USE_COMPLEX
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    cElementMatrix Cs(2, 2);
    mat->setupCs(Cs);

    // std::cout << Cs(0,0) << "\n";
    PetscScalar refVal =
        6.730769230769230e+07 + 6.730769230769233e+04 * PETSC_i;
    ASSERT_EQ(Cs(0, 0), refVal);

    // std::cout << Cs(1, 1) << "\n";
    ASSERT_EQ(Cs(1, 1), refVal);
#else
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    cElementMatrix Cs(2, 2);

    testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH(mat->setupCs(Cs), ".*");
    // FAIL() << "Recompile the code with support for complex numbers.";
#endif
  }

  //! @brief Function to assert - do get bending part of the elasticity matrix
  //! of Mindlin plate
  //! @author Akash Doddamane Lingaraja
  //! @date 17.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_SETUPCB() {
#ifdef PETSC_USE_COMPLEX
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    cElementMatrix Cb(3, 3);
    mat->setupCb(Cb);

    PetscScalar refVal1 =
        1.730769230769231e+02 + 1.730769230769231e-01 * PETSC_i;
    PetscScalar val1 = Cb(0, 0);
    // std::cout << val1 << "	" << refVal1 << "\n";
    ASSERT_NEAR(val1.real(), refVal1.real(), 1e-06);
    ASSERT_NEAR(val1.imag(), refVal1.imag(), 1e-06);

    PetscScalar val2 = Cb(1, 1);
    // std::cout << val2 << "	" << refVal1 << "\n";
    ASSERT_NEAR(val2.real(), refVal1.real(), 1e-06);
    ASSERT_NEAR(val2.imag(), refVal1.imag(), 1e-06);

    PetscScalar refVal2 = 60.576923076923070 + 0.060576923076923 * PETSC_i;
    PetscScalar val3 = Cb(2, 2);
    // std::cout << val3 << "	" << refVal2 << "\n";
    ASSERT_NEAR(val3.real(), refVal2.real(), 1e-06);
    ASSERT_NEAR(val3.imag(), refVal2.imag(), 1e-06);

    PetscScalar refVal3 = 51.923076923076920 + 0.051923076923077 * PETSC_i;
    PetscScalar val4 = Cb(0, 1);
    // std::cout << val4 << "	" << refVal3 << "\n";
    ASSERT_NEAR(val4.real(), refVal3.real(), 1e-06);
    ASSERT_NEAR(val4.imag(), refVal3.imag(), 1e-06);

    PetscScalar val5 = Cb(1, 0);
    // std::cout << val5 << refVal3 << "\n";
    ASSERT_NEAR(val5.real(), refVal3.real(), 1e-06);
    ASSERT_NEAR(val5.imag(), refVal3.imag(), 1e-06);
#else
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    cElementMatrix Cb(3, 3);

    testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH(mat->setupCb(Cb), ".*");
    // FAIL() << "Recompile the code with support for complex numbers.";
#endif
  }

  //! @brief Function to assert - do get bending part of the elasticity matrix
  //! of Mindlin plate
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_SETUPCM() {
#ifdef PETSC_USE_COMPLEX
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    cElementMatrix Cm(3, 3);
    mat->setupCm(Cm);

    PetscScalar refVal1 =
        7.692307692307692e+10 + 7.692307692307693e+07 * PETSC_i;
    PetscScalar val1 = Cm(0, 0);
    // std::cout << val1 << "	" << refVal1 << "\n";
    ASSERT_NEAR(val1.real(), refVal1.real(), 1e-06);
    ASSERT_NEAR(val1.imag(), refVal1.imag(), 1e-06);

    PetscScalar val2 = Cm(1, 1);
    // std::cout << val2 << "	" << refVal1 << "\n";
    ASSERT_NEAR(val2.real(), refVal1.real(), 1e-06);
    ASSERT_NEAR(val2.imag(), refVal1.imag(), 1e-06);

    PetscScalar refVal2 =
        2.692307692307692e+10 + 2.692307692307692e+07 * PETSC_i;
    PetscScalar val3 = Cm(2, 2);
    // std::cout << val3 << "	" << refVal2 << "\n";
    ASSERT_NEAR(val3.real(), refVal2.real(), 1e-06);
    ASSERT_NEAR(val3.imag(), refVal2.imag(), 1e-06);

    PetscScalar refVal3 =
        2.307692307692308e+10 + 2.307692307692308e+07 * PETSC_i;
    PetscScalar val4 = Cm(0, 1);
    // std::cout << val4 << "	" << refVal3 << "\n";
    ASSERT_NEAR(val4.real(), refVal3.real(), 1e-05);
    ASSERT_NEAR(val4.imag(), refVal3.imag(), 1e-06);

    PetscScalar val5 = Cm(1, 0);
    // std::cout << val5 << refVal3 << "\n";
    ASSERT_NEAR(val5.real(), refVal3.real(), 1e-05);
    ASSERT_NEAR(val5.imag(), refVal3.imag(), 1e-06);
#else
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    cElementMatrix Cm(3, 3);

    testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH(mat->setupCm(Cm), ".*");
    // FAIL() << "Recompile the code with support for complex numbers.";
#endif
  }

  //! @brief Function to assert - do get Rho
  //! @author Akash Doddamane Lingaraja
  //! @date 17.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_RHO() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    // std::cout << (mat->getRho()) << "\n";
    ASSERT_EQ((mat->getRho()), 2700.0);
  }

  //! @brief Function to assert - do get Cross-section area
  //! @author Akash Doddamane Lingaraja
  //! @date 17.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_CROSSSECTIONAREA() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    // std::cout << (mat->getA()) << "\n";
    ASSERT_EQ((mat->getA()), 0);
  }

  //! @brief Function to assert - do get shear correction factor
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_ISOLINVISCOELASTIC_KS() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureVisco* mat;
    mat = dynamic_cast<cMaterialStructureVisco*>(ptrMaterial);
    // std::cout << (mat->getKs()) << "\n";
    ASSERT_EQ((mat->getKs()), 5.0 / 6.0);
  }

 private:
  cProblem myTestProblem;
};

// Material type "STR_LIN_VIS_ISO_DIR"

TEST_F(cTestMaterialStructureVisco,
       CheckForCorrectIsoLinViscoElasticDefaultValues) {
  cMaterialStructureVisco mat;
  ASSERT_EQ(mat.getE(), 0.0);
  ASSERT_EQ(mat.getNu(), 0.0);
  ASSERT_EQ(mat.getSpringConst(), 0.0);
  ASSERT_NEAR(mat.getKs(), 5. / 6., 1e-6);
  ASSERT_EQ(mat.getEta(), 0.0);
  ASSERT_EQ(mat.getType(), -1);
}

TEST_F(cTestMaterialStructureVisco,
       checkForCorrectIsoLinViscoElasticGetEOmegaMethod) {
  TEST_ARRANGE(HDF5_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_EOMEGA();
}

TEST_F(cTestMaterialStructureVisco,
       checkForCorrectIsoLinViscoElasticGetEMethod) {
  TEST_ARRANGE(HDF5_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_GETE();
}

TEST_F(cTestMaterialStructureVisco,
       checkForCorrectIsoLinViscoElasticGetEtaMethod) {
  TEST_ARRANGE(HDF5_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_GETETA();
}

TEST_F(cTestMaterialStructureVisco,
       checkForCorrectIsoLinViscoElasticGetNuMethod) {
  TEST_ARRANGE(HDF5_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_GETNU();
}

TEST_F(cTestMaterialStructureVisco,
       checkForCorrectIsoLinViscoElasticGetTypeMethod) {
  TEST_ARRANGE(HDF5_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_GETTYPE();
}

TEST_F(cTestMaterialStructureVisco,
       checkForCorrectIsoLinViscoElasticSetEMethod) {
  TEST_ARRANGE(HDF5_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_SETE();
}

TEST_F(cTestMaterialStructureVisco,
       checkForCorrectIsoLinViscoElasticSetEtaMethod) {
  TEST_ARRANGE(HDF5_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_SETETA();
}

TEST_F(cTestMaterialStructureVisco,
       checkForCorrectIsoLinViscoElasticSetNuMethod) {
  TEST_ARRANGE(HDF5_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_SETNU();
}

TEST_F(cTestMaterialStructureVisco,
       checkForCorrectIsoLinViscoElasticSetTypeMethod) {
  TEST_ARRANGE(HDF5_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_SETTYPE();
}

TEST_F(cTestMaterialStructureVisco,
       checkForCorrectIsoLinViscoElasticGetFilenameEtaMethod) {
  TEST_ARRANGE(HDF5_FREQ_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_GETFILENAMEETA();
}

TEST_F(cTestMaterialStructureVisco,
       checkForCorrectIsoLinViscoElasticGetFilenameEMethod) {
  TEST_ARRANGE(HDF5_FREQ_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_GETFILENAMEE();
}

TEST_F(cTestMaterialStructureVisco,
       checkForCorrectIsoLinViscoElasticSpringConstantValue) {
  TEST_ARRANGE(HDF5_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_SPRINGCONST();
}

TEST_F(cTestMaterialStructureVisco,
       checkForCorrectIsoLinViscoElasticThicknessValue) {
  TEST_ARRANGE(HDF5_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_THICKNESS();
}

TEST_F(cTestMaterialStructureVisco,
       checkForCorrectIsoLinViscoElasticMeanCoordinateOfElement) {
  TEST_ARRANGE(HDF5_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_CALCMEANCOORDINATE();
}

TEST_F(cTestMaterialStructureVisco,
       checkForCorrectIsoLinViscoElasticDistanceBetweenTwoPoints) {
  TEST_ASSERT_ISOLINVISCOELASTIC_CALCDIST();
}

TEST_F(cTestMaterialStructureVisco, checkForCorrectIsoLinViscoElasticCsValue) {
  TEST_ARRANGE(HDF5_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_SETUPCS();
}

TEST_F(cTestMaterialStructureVisco, checkForCorrectIsoLinViscoElasticCbValue) {
  TEST_ARRANGE(HDF5_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_SETUPCB();
}

TEST_F(cTestMaterialStructureVisco, checkForCorrectIsoLinViscoElasticCmValue) {
  TEST_ARRANGE(HDF5_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_SETUPCM();
}

TEST_F(cTestMaterialStructureVisco, checkForCorrectIsoLinViscoElasticRhoValue) {
  TEST_ARRANGE(HDF5_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_RHO();
}

TEST_F(cTestMaterialStructureVisco,
       checkForCorrectIsoLinViscoElasticCrossSectionAreaValue) {
  TEST_ARRANGE(HDF5_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_CROSSSECTIONAREA();
}

TEST_F(cTestMaterialStructureVisco, checkForCorrectIsoLinViscoElasticKsValue) {
  TEST_ARRANGE(HDF5_STR_LIN_VIS_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINVISCOELASTIC_KS();
}
