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

#include "../../source/shape/shapequad4.h"

#ifdef HAVE_GTEST
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#endif // HAVE_GTEST

//! @brief Mock class for cShapeQuad4
//! @author Akash Doddamane Lingaraja
//! @date 04.05.2022
class cTestShapeQuad4 : public testing::Test, public cShapeQuad4 {
private:

public:
	
	//! @brief GTEST Setup
	//! @author Akash Doddamane Lingaraja
	//! @date 04.05.2022	
	void SetUp() {
		// do nothing
	}

	
	//! @brief GTEST Teardown
	//! @author Akash Doddamane Lingaraja
	//! @date 04.05.2022	
	void TearDown() {
		// do nothing
	}

	
	//! @brief Function to arrange, act, assert - do get Shape function
	//! @author Akash Doddamane Lingaraja
	//! @date 04.05.2022	
	void TEST_GETSHAPEFUNCTIONS() {
		// Arrange
		cShapeQuad4 element;
		cArray3d N;
		int gaussianPoints = 48;

		// Act
		N = element.getShapeFunctions();
		PetscReal* data = N.getBase();
		
		/*
		std::cout.precision(8);
		for (int i = 0; i < gaussianPoints; i++)
			std::cout << data[i] << "\n";
		std::cout << std::endl;
		*/

		// Assert
		double refData[48] = { 0.044658199,	0.16666667,	0.62200847,	0.16666667,	-0.10566243,  0.10566243,  0.39433757,	-0.39433757,
										  -0.10566243,	-0.39433757,  0.39433757,	0.10566243,	0.16666667,	0.62200847,	0.16666667,	0.044658199,
										  -0.39433757,	0.39433757,	0.10566243,	-0.10566243,	-0.10566243,	-0.39433757,	0.39433757,
										   0.10566243,	0.16666667,	0.044658199,	0.16666667,	0.62200847,	-0.10566243,	0.10566243,	0.39433757,
										  -0.39433757,	-0.39433757,	-0.10566243,	0.10566243,	0.39433757,	0.62200847,	0.16666667,	0.044658199,
										   0.16666667,	-0.39433757,	0.39433757,	0.10566243,	-0.10566243,	-0.39433757,	-0.10566243,
										   0.10566243,	0.39433757 };
		for (size_t i = 0; i < gaussianPoints; i++)
			ASSERT_NEAR(data[i], refData[i], 1e-6);			
	}

	
	//! @brief Function to arrange, act, assert - do get Shape of the element
	//! @author Akash Doddamane Lingaraja
	//! @date 04.05.2022	
	void TEST_GETELEMENTSHAPE() {
		// Arrange and Act 
		int shape = cShapeQuad4::getElementShape();

		// Assert
		ASSERT_EQ(shape, Quadrilateral);
	}

	
	//! @brief Function to arrange, act, assert - do get number of faces of Tria3 elements
	//! @author Akash Doddamane Lingaraja
	//! @date 04.05.2022	
	void TEST_GETNUMBEROFFACES() {
		// Arrange and Act 
		short faces = cShapeQuad4::getNumberOfFaces();

		// Assert
		ASSERT_EQ(faces, 4);
	}

	
	//! @brief Function to arrange, act, assert - do get Numbering of one element's face
	//! @author Akash Doddamane Lingaraja
	//! @date 04.05.2022	
	void TEST_GETINDICESOFFACENODES() {
		// Arrange
		cShapeQuad4 element;
		std::vector<short> nodesFace0 { 1, 0 };
		std::vector<short> nodesFace1 { 2, 1 };
		std::vector<short> nodesFace2 { 3, 2 };
		std::vector<short> nodesFace3 { 0, 3 };

		// Assert
		ASSERT_EQ(element.getIndicesOfFaceNodes(0), nodesFace0);
		ASSERT_EQ(element.getIndicesOfFaceNodes(1), nodesFace1);
		ASSERT_EQ(element.getIndicesOfFaceNodes(2), nodesFace2);
		ASSERT_EQ(element.getIndicesOfFaceNodes(3), nodesFace3);
		ASSERT_ANY_THROW(element.getIndicesOfFaceNodes(4));
	}

	
	//! @brief Function to arrange, act, assert - do get number of nodes that describe one face of the element
	//! @author Akash Doddamane Lingaraja
	//! @date 04.05.2022	
	void TEST_GETNUMBEROFNODESPERFACE() {
		// Arrange
		short nodes = 2;
		// Assert
		ASSERT_EQ(getNumberOfNodesPerFace(), nodes);
	}

	
	//! @brief Function to arrange, act, assert - do Initialize shape functions for PoroPlateKienzler
	//! @author Akash Doddamane Lingaraja
	//! @date 04.05.2022	
	void TEST_INITIALIZESHAPEFUNCTIONS_N_MAP() {
		// Arrange
		cShapeQuad4* element = new cShapeQuad4;
		int gaussianPoints = 48;

		// Act
		initializeShapeFunctions_N_map();
		//cArray3d N = cShapeQuad4::N_map;
		PetscReal* data = cShapeQuad4::N_map.getBase();

		// Assert
		for (int i = 0; i < gaussianPoints; i++)
			ASSERT_EQ(data[i], 0.);
	}
};

TEST_F(cTestShapeQuad4, checkForCorrectShapeFunctionValues) {
	TEST_GETSHAPEFUNCTIONS();
}

TEST_F(cTestShapeQuad4, checkForCorrectShapeOfTheElement) {
	TEST_GETELEMENTSHAPE();
}

TEST_F(cTestShapeQuad4, checkForCorrectNumberOfFacesOfTheQuad4Element) {
	TEST_GETNUMBEROFFACES();
}

TEST_F(cTestShapeQuad4, checkForCorrectNumberingOfFacesOfTheQuad4Element) {
	TEST_GETINDICESOFFACENODES();
}

TEST_F(cTestShapeQuad4, checkForCorrectNumberOfNodesPerFaceOfTheQuad4Element) {
	TEST_GETNUMBEROFNODESPERFACE();
}

TEST_F(cTestShapeQuad4, checkForCorrectIntializationOfNMapOfTheQuad4Element) {
	TEST_INITIALIZESHAPEFUNCTIONS_N_MAP();
}