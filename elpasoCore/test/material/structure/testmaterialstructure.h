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

#include "../../../source/misc/parser/femparserhdf5.h"
#include "../../resources/properties/testproperties.h"

#ifdef HAVE_GTEST
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#endif  // HAVE_GTEST

//! @brief Mock class for cMaterialStructure
//! @author Akash Doddamane Lingaraja
//! @date 21.04.2022
class cTestMaterialStructure : private cFemParserHDF5, public testing::Test {
 private:
  cProblem myTestProblem;

 public:
  //! @brief GTEST Setup
  //! @author Akash Doddamane Lingaraja
  //! @date 21.04.2022
  void SetUp() {
    // do nothing
  }

  //! @brief GTEST Teardown
  //! @author Akash Doddamane Lingaraja
  //! @date 21.04.2022
  void TearDown() {
    // close HDF5 File
#ifdef HAVE_HDF5
    cInputSingletonH5::getInstance()->closeContainer();
#endif
  }

  //! @brief Function to arrange - do parse the HDF5 file
  //! @author Akash Doddamane Lingaraja
  //! @date 21.04.2022
  void TEST_ARRANGE(std::string filename) {
    std::string elpasoTestResourceDir = ELPASO_TEST_RESOURCE_DIR;
    std::string testFileHDF5 =
        elpasoTestResourceDir +
        "unitHDF5/material/structure/isotrop/linear/elastic/" + filename;

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
  //! @date 21.04.2022
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

  //! @brief Function to assert - do get mass of structural elements
  //! @author Akash Doddamane Lingaraja
  //! @date 21.04.2022
  void TEST_ASSERT_GETMASS() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructure* mat;
    mat = dynamic_cast<cMaterialStructure*>(ptrMaterial);

    ASSERT_EQ((mat->getMass()), -1.);
  }

  //! @brief Function to assert - do get stiffness trans x
  //! @author Akash Doddamane Lingaraja
  //! @date 21.04.2022
  void TEST_ASSERT_GETCX() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructure* mat;
    mat = dynamic_cast<cMaterialStructure*>(ptrMaterial);

    ASSERT_EQ((mat->getCx()), -1.);
  }

  //! @brief Function to assert - do get stiffness trans y
  //! @author Akash Doddamane Lingaraja
  //! @date 21.04.2022
  void TEST_ASSERT_GETCY() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructure* mat;
    mat = dynamic_cast<cMaterialStructure*>(ptrMaterial);

    ASSERT_EQ((mat->getCy()), -1.);
  }

  //! @brief Function to assert - do get stiffness trans z
  //! @author Akash Doddamane Lingaraja
  //! @date 21.04.2022
  void TEST_ASSERT_GETCZ() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructure* mat;
    mat = dynamic_cast<cMaterialStructure*>(ptrMaterial);

    ASSERT_EQ((mat->getCz()), -1.);
  }

  //! @brief Function to assert - do get stiffness rot x
  //! @author Akash Doddamane Lingaraja
  //! @date 21.04.2022
  void TEST_ASSERT_GETCRX() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructure* mat;
    mat = dynamic_cast<cMaterialStructure*>(ptrMaterial);

    ASSERT_EQ((mat->getCrx()), -1.);
  }

  //! @brief Function to assert - do get stiffness rot y
  //! @author Akash Doddamane Lingaraja
  //! @date 21.04.2022
  void TEST_ASSERT_GETCRY() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructure* mat;
    mat = dynamic_cast<cMaterialStructure*>(ptrMaterial);

    ASSERT_EQ((mat->getCry()), -1.);
  }

  //! @brief Function to assert - do get stiffness rot z
  //! @author Akash Doddamane Lingaraja
  //! @date 21.04.2022
  void TEST_ASSERT_GETCRZ() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructure* mat;
    mat = dynamic_cast<cMaterialStructure*>(ptrMaterial);

    ASSERT_EQ((mat->getCrz()), -1.);
  }

  //! @brief Function to assert - do get Eta
  //! @author Akash Doddamane Lingaraja
  //! @date 21.04.2022
  void TEST_ASSERT_GETETASPRING() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructure* mat;
    mat = dynamic_cast<cMaterialStructure*>(ptrMaterial);

    ASSERT_EQ((mat->getEtaSpring()), -1.);
  }

  //! @brief Function to assert - do get moment of inertia
  //! @author Akash Doddamane Lingaraja
  //! @date 21.04.2022
  void TEST_ASSERT_GETI() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructure* mat;
    mat = dynamic_cast<cMaterialStructure*>(ptrMaterial);

    ASSERT_EQ((mat->getI()), 0.0);
  }

  //! @brief Function to assert - do get moment of inertia in X
  //! @author Akash Doddamane Lingaraja
  //! @date 21.04.2022
  void TEST_ASSERT_GETIX() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructure* mat;
    mat = dynamic_cast<cMaterialStructure*>(ptrMaterial);

    ASSERT_EQ((mat->getIx()), 0.0);
  }

  //! @brief Function to assert - do get moment of inertia in Y
  //! @author Akash Doddamane Lingaraja
  //! @date 21.04.2022
  void TEST_ASSERT_GETIY() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructure* mat;
    mat = dynamic_cast<cMaterialStructure*>(ptrMaterial);

    ASSERT_EQ((mat->getIy()), 0.0);
  }

  //! @brief Function to assert - do get moment of inertia in Z
  //! @author Akash Doddamane Lingaraja
  //! @date 21.04.2022
  void TEST_ASSERT_GETIZ() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructure* mat;
    mat = dynamic_cast<cMaterialStructure*>(ptrMaterial);

    ASSERT_EQ((mat->getIz()), 0.0);
  }

  //! @brief Function to assert - do get thickness
  //! @author Akash Doddamane Lingaraja
  //! @date 21.04.2022
  void TEST_ASSERT_GETFI() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructure* mat;
    mat = dynamic_cast<cMaterialStructure*>(ptrMaterial);

    ASSERT_EQ((mat->getFi()), 0.0);
  }

  //! @brief Function to assert - do get bulkmodulus
  //! @author Akash Doddamane Lingaraja
  //! @date 21.04.2022
  void TEST_ASSERT_GETRHOLAMBDA() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructure* mat;
    mat = dynamic_cast<cMaterialStructure*>(ptrMaterial);

    cPoint point;
    PetscReal rho = 2.0;
    ASSERT_EQ((mat->getRhoLambda(point, rho)), 0.0);
  }
};

TEST_F(cTestMaterialStructure, checkForCorrectMassOfStructuralElements) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_GETMASS();
}

TEST_F(cTestMaterialStructure,
       checkForCorrectStiffnessTranslationInXOfStructuralElements) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_GETCX();
}

TEST_F(cTestMaterialStructure,
       checkForCorrectStiffnessTranslationInYOfStructuralElements) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_GETCY();
}

TEST_F(cTestMaterialStructure,
       checkForCorrectStiffnessTranslationInZOfStructuralElements) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_GETCZ();
}

TEST_F(cTestMaterialStructure,
       checkForCorrectStiffnessRotationInXOfStructuralElements) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_GETCRX();
}

TEST_F(cTestMaterialStructure,
       checkForCorrectStiffnessRotationInYOfStructuralElements) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_GETCRY();
}

TEST_F(cTestMaterialStructure,
       checkForCorrectStiffnessRotationInZOfStructuralElements) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_GETCRZ();
}

TEST_F(cTestMaterialStructure, checkForCorrectETAOfStructuralElements) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_GETETASPRING();
}

TEST_F(cTestMaterialStructure,
       checkForCorrectMomentOfInertiaOfStructuralElements) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_GETI();
}

TEST_F(cTestMaterialStructure,
       checkForCorrectMomentOfInertiaInXOfStructuralElements) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_GETIX();
}

TEST_F(cTestMaterialStructure,
       checkForCorrectMomentOfInertiaInYOfStructuralElements) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_GETIY();
}

TEST_F(cTestMaterialStructure,
       checkForCorrectMomentOfInertiaInZOfStructuralElements) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_GETIZ();
}

TEST_F(cTestMaterialStructure, checkForCorrectThicknessOfStructuralElements) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_GETFI();
}

TEST_F(cTestMaterialStructure, checkForCorrectRhoLambdaOfStructuralElements) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_GETRHOLAMBDA();
}