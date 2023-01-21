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

//! @brief Mock class for cMaterialStructureIsotrop
//! @author Akash Doddamane Lingaraja
//! @date 07.01.2022
class cTestMaterialStructureIsotrop : private cFemParserHDF5,
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

  //! @brief Function to assert - do get Young's modulus
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_ISOLINELASTIC_YOUNGSMODULUS() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureIsotrop* mat;
    mat = dynamic_cast<cMaterialStructureIsotrop*>(ptrMaterial);
    // std::cout << (mat->getE()) << "\n";
    ASSERT_EQ((mat->getE()), 7e10);
  }

  //! @brief Function to assert - do get frequency dependent Young's modulus
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_ISOLINELASTIC_EOMEGA() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureIsotrop* mat;
    mat = dynamic_cast<cMaterialStructureIsotrop*>(ptrMaterial);
    // std::cout << (mat->getEOmega()) << "\n";
    ASSERT_EQ((mat->getEOmega()), 7e10);
  }

  //! @brief Function to assert - do get Poisson's ratio
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_ISOLINELASTIC_NU() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureIsotrop* mat;
    mat = dynamic_cast<cMaterialStructureIsotrop*>(ptrMaterial);
    // std::cout << (mat->getNu()) << "\n";
    ASSERT_EQ((mat->getNu()), 0.3);
  }

  //! @brief Function to assert - do get Spring constant
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_ISOLINELASTIC_SPRINGCONST() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureIsotrop* mat;
    mat = dynamic_cast<cMaterialStructureIsotrop*>(ptrMaterial);
    // std::cout << (mat->getSpringConst()) << "\n";
    ASSERT_EQ((mat->getSpringConst()), 7e10);
  }

  //! @brief Function to assert - do get Thickness
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_ISOLINELASTIC_THICKNESS() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureIsotrop* mat;
    mat = dynamic_cast<cMaterialStructureIsotrop*>(ptrMaterial);
    // std::cout << (mat->getT()) << "\n";
    ASSERT_EQ((mat->getT()), 0.003);
  }

  //! @brief Function to assert - do get Thickness of beam element
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ASSERT_BEAM_ISOLINELASTIC_THICKNESS() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureIsotrop* mat;
    mat = dynamic_cast<cMaterialStructureIsotrop*>(ptrMaterial);
    // std::cout << (mat->getT()) << "\n";
    ASSERT_EQ((mat->getT()), 0.0);
  }

  //! @brief Function to assert - do get shear part of the elasticity matrix of
  //! Mindlin plate
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_ISOLINELASTIC_SETUPCS() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureIsotrop* mat;
    mat = dynamic_cast<cMaterialStructureIsotrop*>(ptrMaterial);
    cElementMatrix Cs(2, 2);
    mat->setupCs(Cs);

    // std::cout << Cs(0,0) << "\n";
    double refVal = 6.730769230769230e+07;
    ASSERT_EQ(Cs(0, 0).real(), refVal);

    // std::cout << Cs(1, 1) << "\n";
    ASSERT_EQ(Cs(1, 1).real(), refVal);
  }

  //! @brief Function to assert - do get shear part of the elasticity matrix of
  //! Mindlin plate for beam element
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022

  void TEST_ASSERT_BEAM_ISOLINELASTIC_SETUPCS() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureIsotrop* mat;
    mat = dynamic_cast<cMaterialStructureIsotrop*>(ptrMaterial);
    cElementMatrix Cs(2, 2);
    mat->setupCs(Cs);

    // std::cout << Cs(0,0) << "\n";
    ASSERT_EQ(Cs(0, 0).real(), 0.0);
    ASSERT_EQ(Cs(1, 1).real(), 0.0);
  }

  //! @brief Function to assert - do get bending part of the elasticity matrix
  //! of Mindlin plate
  //! @author Akash Doddamane Lingaraja
  //! @date 17.01.2022
  void TEST_ASSERT_ISOLINELASTIC_SETUPCB() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureIsotrop* mat;
    mat = dynamic_cast<cMaterialStructureIsotrop*>(ptrMaterial);
    cElementMatrix Cb(3, 3);
    mat->setupCb(Cb);

    double refVal1 = 1.730769230769231e+02;
    // std::cout << Cb(0, 0) << "\n";
    ASSERT_NEAR(Cb(0, 0).real(), refVal1, 1e-06);

    // std::cout << Cb(1, 1) << "\n";
    ASSERT_NEAR(Cb(1, 1).real(), refVal1, 1e-06);

    double refVal2 = 60.576923076923070;
    // std::cout << Cb(2, 2) << "\n";
    ASSERT_NEAR(Cb(2, 2).real(), refVal2, 1e-06);

    double refVal3 = 51.923076923076920;
    // std::cout << Cb(0, 1) << "\n";
    ASSERT_NEAR(Cb(0, 1).real(), refVal3, 1e-06);

    // std::cout << Cb(1, 0) << "\n";
    ASSERT_NEAR(Cb(1, 0).real(), refVal3, 1e-06);
  }

  //! @brief Function to assert - do get bending part of the elasticity matrix
  //! of Mindlin plate for beam element
  //! @author Akash Doddamane Lingaraja
  //! @date 01.03.2022
  void TEST_ASSERT_BEAM_ISOLINELASTIC_SETUPCB() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureIsotrop* mat;
    mat = dynamic_cast<cMaterialStructureIsotrop*>(ptrMaterial);
    cElementMatrix Cb(3, 3);
    mat->setupCb(Cb);

    ASSERT_EQ(Cb(0, 0).real(), 0.0);
    ASSERT_EQ(Cb(1, 1).real(), 0.0);
    ASSERT_EQ(Cb(2, 2).real(), 0.0);
    ASSERT_EQ(Cb(0, 1).real(), 0.0);
    ASSERT_EQ(Cb(1, 0).real(), 0.0);
  }

  //! @brief Function to assert - do get bending part of the elasticity matrix
  //! of Mindlin plate
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_ISOLINELASTIC_SETUPCM() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureIsotrop* mat;
    mat = dynamic_cast<cMaterialStructureIsotrop*>(ptrMaterial);
    cElementMatrix Cm(3, 3);
    mat->setupCm(Cm);

    double refVal1 = 7.692307692307692e+10;
    // std::cout << Cm(0, 0) << "\n";
    ASSERT_NEAR(Cm(0, 0).real(), refVal1, 1e-06);

    // std::cout << Cm(1, 1) << "\n";
    ASSERT_NEAR(Cm(1, 1).real(), refVal1, 1e-06);

    double refVal2 = 2.692307692307692e+10;
    // std::cout << Cm(2, 2) << "\n";
    ASSERT_NEAR(Cm(2, 2).real(), refVal2, 1e-06);

    double refVal3 = 2.307692307692308e+10;
    // std::cout << Cm(0, 1) << "\n";
    ASSERT_NEAR(Cm(0, 1).real(), refVal3, 1e-05);

    // std::cout <<  Cm(1, 0) << "\n";
    ASSERT_NEAR(Cm(1, 0).real(), refVal3, 1e-05);
  }

  //! @brief Function to assert - do get Rho
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_ISOLINELASTIC_RHO() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureIsotrop* mat;
    mat = dynamic_cast<cMaterialStructureIsotrop*>(ptrMaterial);
    // std::cout << (mat->getRho()) << "\n";
    ASSERT_EQ((mat->getRho()), 2700.0);
  }

  //! @brief Function to assert - do get Cross-section area
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_ISOLINELASTIC_CROSSSECTIONAREA() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureIsotrop* mat;
    mat = dynamic_cast<cMaterialStructureIsotrop*>(ptrMaterial);
    // std::cout << (mat->getA()) << "\n";
    ASSERT_EQ((mat->getA()), 0);
  }

  //! @brief Function to assert - do get Cross-section area for beam elememt
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_BEAM_ISOLINELASTIC_CROSSSECTIONAREA() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureIsotrop* mat;
    mat = dynamic_cast<cMaterialStructureIsotrop*>(ptrMaterial);

    ASSERT_EQ((mat->getA()), 0.0001);
  }

  //! @brief Function to assert - do get shear correction factor
  //! @author Akash Doddamane Lingaraja
  //! @date 19.01.2022
  void TEST_ASSERT_ISOLINELASTIC_KS() {
    cMaterial* ptrMaterial = myTestProblem.getMesh()->getMaterial(1);
    cMaterialStructureIsotrop* mat;
    mat = dynamic_cast<cMaterialStructureIsotrop*>(ptrMaterial);
    // std::cout << (mat->getKs()) << "\n";
    ASSERT_EQ((mat->getKs()), 5.0 / 6.0);
  }

 private:
  cProblem myTestProblem;
};

TEST_F(cTestMaterialStructureIsotrop,
       CheckForCorrectIsoLinElasticDefaultValues) {
  cMaterialStructureIsotrop mat;
  ASSERT_EQ(mat.getE(), 0.0);
  ASSERT_EQ(mat.getNu(), 0.0);
  ASSERT_EQ(mat.getSpringConst(), 0.0);
  ASSERT_NEAR(mat.getKs(), 5. / 6., 1e-6);
}

TEST_F(cTestMaterialStructureIsotrop, checkForCorrectIsoLinElasticGetEMethod) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTIC_YOUNGSMODULUS();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectIsoLinElasticGetEOmegaMethod) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTIC_EOMEGA();
}

TEST_F(cTestMaterialStructureIsotrop, checkForCorrectIsoLinElasticGetNuMethod) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTIC_NU();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectIsoLinElasticSpringConstantValue) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTIC_SPRINGCONST();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectIsoLinElasticThicknessValue) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTIC_THICKNESS();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectIsoLinViscoElasticCsValue) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTIC_SETUPCS();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectIsoLinViscoElasticCbValue) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTIC_SETUPCB();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectIsoLinViscoElasticCmValue) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTIC_SETUPCM();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectIsoLinViscoElasticRhoValue) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTIC_RHO();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectIsoLinViscoElasticCrossSectionAreaValue) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTIC_CROSSSECTIONAREA();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectIsoLinViscoElasticKsValue) {
  TEST_ARRANGE(HDF5_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTIC_KS();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectBeamIsoLinElasticGetEMethod) {
  TEST_ARRANGE(HDF5_BEAM_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTIC_YOUNGSMODULUS();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectBeamIsoLinElasticGetEOmegaMethod) {
  TEST_ARRANGE(HDF5_BEAM_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTIC_EOMEGA();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectIsoBeamLinElasticGetNuMethod) {
  TEST_ARRANGE(HDF5_BEAM_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTIC_NU();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectBeamIsoLinElasticSpringConstantValue) {
  TEST_ARRANGE(HDF5_BEAM_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTIC_SPRINGCONST();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectBeamIsoLinElasticThicknessValue) {
  TEST_ARRANGE(HDF5_BEAM_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_BEAM_ISOLINELASTIC_THICKNESS();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectBeamIsoLinViscoElasticCsValue) {
  TEST_ARRANGE(HDF5_BEAM_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_BEAM_ISOLINELASTIC_SETUPCS();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectBeamIsoLinViscoElasticCbValue) {
  TEST_ARRANGE(HDF5_BEAM_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_BEAM_ISOLINELASTIC_SETUPCB();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectBeamIsoLinViscoElasticCmValue) {
  TEST_ARRANGE(HDF5_BEAM_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTIC_SETUPCM();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectBeamIsoLinViscoElasticRhoValue) {
  TEST_ARRANGE(HDF5_BEAM_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTIC_RHO();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectBeamIsoLinViscoElasticCrossSectionAreaValue) {
  TEST_ARRANGE(HDF5_BEAM_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_BEAM_ISOLINELASTIC_CROSSSECTIONAREA();
}

TEST_F(cTestMaterialStructureIsotrop,
       checkForCorrectBeamIsoLinViscoElasticKsValue) {
  TEST_ARRANGE(HDF5_BEAM_STR_LIN_ELA_ISO_DIR);
  TEST_ACT();
  TEST_ASSERT_ISOLINELASTIC_KS();
}