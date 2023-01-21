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

//! @brief Mock class for cMaterialSpring
//! @author Akash Doddamane Lingaraja
//! @date 01.03.2022
class cTestMaterialSpring : private cFemParserHDF5, public testing::Test {
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
    std::string testFileHDF5 =
        elpasoTestResourceDir +
        "unitHDF5/material/structure/isotrop/linear/elastic/" +
        HDF5_STR_LIN_SPR_ORT_DIR;

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

  //! @brief Function to assert - do get stiffness trans x
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ASSERT_ISOLINELASTICSPRING_GETCX() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialSpring* mat;
    mat = dynamic_cast<cMaterialSpring*>(ptrMaterial);

#ifdef PETSC_USE_COMPLEX
    PetscScalar refVal = 100.0 * (1 + PETSC_i * 0.003);
    ASSERT_EQ((mat->getCx()), refVal);
#else
    ASSERT_EQ((mat->getCx()), 100.0);
#endif
  }

  //! @brief Function to assert - do get stiffness trans y
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ASSERT_ISOLINELASTICSPRING_GETCY() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialSpring* mat;
    mat = dynamic_cast<cMaterialSpring*>(ptrMaterial);

#ifdef PETSC_USE_COMPLEX
    PetscScalar refVal = 200.0 * (1 + PETSC_i * 0.003);
    ASSERT_EQ((mat->getCy()), refVal);
#else
    ASSERT_EQ((mat->getCy()), 100.0);
#endif
  }

  //! @brief Function to assert - do get stiffness trans z
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ASSERT_ISOLINELASTICSPRING_GETCZ() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialSpring* mat;
    mat = dynamic_cast<cMaterialSpring*>(ptrMaterial);

#ifdef PETSC_USE_COMPLEX
    PetscScalar refVal = 300.0 * (1 + PETSC_i * 0.003);
    ASSERT_EQ((mat->getCz()), refVal);
#else
    ASSERT_EQ((mat->getCz()), 300.0);
#endif
  }

  //! @brief Function to assert - do get stiffness rot x
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ASSERT_ISOLINELASTICSPRING_GETCRX() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialSpring* mat;
    mat = dynamic_cast<cMaterialSpring*>(ptrMaterial);

#ifdef PETSC_USE_COMPLEX
    PetscScalar refVal = 400.0 * (1 + PETSC_i * 0.003);
    ASSERT_EQ((mat->getCrx()), refVal);
#else
    ASSERT_EQ((mat->getCrx()), 400.0);
#endif
  }

  //! @brief Function to assert - do get stiffness rot y
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ASSERT_ISOLINELASTICSPRING_GETCRY() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialSpring* mat;
    mat = dynamic_cast<cMaterialSpring*>(ptrMaterial);

#ifdef PETSC_USE_COMPLEX
    PetscScalar refVal = 500.0 * (1 + PETSC_i * 0.003);
    ASSERT_EQ((mat->getCry()), refVal);
#else
    ASSERT_EQ((mat->getCry()), 500.0);
#endif
  }

  //! @brief Function to assert - do get stiffness rot z
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ASSERT_ISOLINELASTICSPRING_GETCRZ() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialSpring* mat;
    mat = dynamic_cast<cMaterialSpring*>(ptrMaterial);

#ifdef PETSC_USE_COMPLEX
    PetscScalar refVal = 600.0 * (1 + PETSC_i * 0.003);
    ASSERT_EQ((mat->getCrz()), refVal);
#else
    ASSERT_EQ((mat->getCrz()), 600.0);
#endif
  }

  //! @brief Function to assert - do get Eta
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ASSERT_ISOLINELASTICSPRING_GETETASPRING() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialSpring* mat;
    mat = dynamic_cast<cMaterialSpring*>(ptrMaterial);

    ASSERT_EQ((mat->getEtaSpring()), 0.003);
  }

  //! @brief Function to assert - do get Type of material
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ASSERT_ISOLINELASTICSPRING_GETTYPE() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialSpring* mat;
    mat = dynamic_cast<cMaterialSpring*>(ptrMaterial);

    ASSERT_EQ((mat->getType()), 1);
  }

  //! @brief Function to assert - do set Eta
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ASSERT_ISOLINELASTICSPRING_SETTYPE() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialSpring* mat;
    mat = dynamic_cast<cMaterialSpring*>(ptrMaterial);

    mat->setType(3);
    ASSERT_EQ((mat->getType()), 3);
  }

  //! @brief Function to assert - do get Young's modulus
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ASSERT_ISOLINELASTICSPRING_GETE() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialSpring* mat;
    mat = dynamic_cast<cMaterialSpring*>(ptrMaterial);

    ASSERT_EQ((mat->getE()), 0.0);
  }

  //! @brief Function to assert - do get frequency dependent Young's modulus
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ASSERT_ISOLINELASTICSPRING_GETEOMEGA() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialSpring* mat;
    mat = dynamic_cast<cMaterialSpring*>(ptrMaterial);

    ASSERT_EQ((mat->getEOmega()), 0.0);
  }

  //! @brief Function to assert - do get possoins ratio
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ASSERT_ISOLINELASTICSPRING_GETNU() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialSpring* mat;
    mat = dynamic_cast<cMaterialSpring*>(ptrMaterial);

    ASSERT_EQ((mat->getNu()), 0.0);
  }

  //! @brief Function to assert - do get spring constant for structural
  //! interface elements
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ASSERT_ISOLINELASTICSPRING_GETSPRINGCONST() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialSpring* mat;
    mat = dynamic_cast<cMaterialSpring*>(ptrMaterial);

    ASSERT_EQ((mat->getSpringConst()), 0.0);
  }

  //! @brief Function to assert - do get thickness
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ASSERT_ISOLINELASTICSPRING_GETT() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialSpring* mat;
    mat = dynamic_cast<cMaterialSpring*>(ptrMaterial);

    ASSERT_DEATH(mat->getT(), "");
  }
};

using cDeathTestMaterialSpring = cTestMaterialSpring;

TEST_F(cTestMaterialSpring, checkForCorectLinSpringDefaultValues) {
  cMaterialSpring mat;
  ASSERT_EQ(mat.getCx(), 0.0);
  ASSERT_EQ(mat.getCy(), 0.0);
  ASSERT_EQ(mat.getCz(), 0.0);
  ASSERT_EQ(mat.getCrx(), 0.0);
  ASSERT_EQ(mat.getCry(), 0.0);
  ASSERT_EQ(mat.getCrz(), 0.0);
  ASSERT_EQ(mat.getEtaSpring(), 0.0);
  ASSERT_EQ(mat.getType(), -1);
}

TEST_F(cTestMaterialSpring, checkForCorectLinSpringStiffnessTranslationInXDir) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTICSPRING_GETCX();
}

TEST_F(cTestMaterialSpring, checkForCorectLinSpringStiffnessTranslationInYDir) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTICSPRING_GETCY();
}

TEST_F(cTestMaterialSpring, checkForCorectLinSpringStiffnessTranslationInZDir) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTICSPRING_GETCZ();
}

TEST_F(cTestMaterialSpring, checkForCorectLinSpringStiffnessRotationInXDir) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTICSPRING_GETCRX();
}

TEST_F(cTestMaterialSpring, checkForCorectLinSpringStiffnessRotationInYDir) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTICSPRING_GETCRY();
}

TEST_F(cTestMaterialSpring, checkForCorectLinSpringStiffnessRotationInZDir) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTICSPRING_GETCRZ();
}

TEST_F(cTestMaterialSpring, checkForCorectLinSpringEtaValue) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTICSPRING_GETETASPRING();
}

TEST_F(cTestMaterialSpring, checkForCorectLinSpringTypeValue) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTICSPRING_GETTYPE();
}

TEST_F(cTestMaterialSpring, checkForCorectLinSpringSetTypeValue) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTICSPRING_SETTYPE();
}

TEST_F(cTestMaterialSpring, checkForCorectLinSpringYoungsModulusValue) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTICSPRING_GETE();
}

TEST_F(cTestMaterialSpring,
       checkForCorectLinSpringFreqDependentYoungsModulusValue) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTICSPRING_GETEOMEGA();
}

TEST_F(cTestMaterialSpring, checkForCorectLinSpringPoissonsRatioValue) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTICSPRING_GETNU();
}

TEST_F(cTestMaterialSpring, checkForCorectLinSpringSpringConstantValue) {
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTICSPRING_GETSPRINGCONST();
}

TEST_F(cDeathTestMaterialSpring, checkForCorectLinSpringThickness) {
  testing::FLAGS_gtest_death_test_style = "threadsafe";
  TEST_ARRANGE();
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTICSPRING_GETT();
}