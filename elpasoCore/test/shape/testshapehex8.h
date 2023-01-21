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

#include "../../source/shape/shapehex8.h"

#ifdef HAVE_GTEST
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#endif // HAVE_GTEST

//! @brief Mock class for cShapeHex8
//! @author Akash Doddamane Lingaraja
//! @date 10.05.2022
class cTestShapeHex8 : public testing::Test, public cShapeHex8 {
private:

public:
	
	//! @brief GTEST Setup
	//! @author Akash Doddamane Lingaraja
	//! @date 10.05.2022
	void SetUp() {
		// do nothing
	}

	
	//! @brief GTEST Teardown
	//! @author Akash Doddamane Lingaraja
	//! @date 10.05.2022	
	void TearDown() {
		// do nothing
	}

	
	//! @brief Function to arrange, act, assert - do get Shape function
	//! @author Akash Doddamane Lingaraja
	//! @date 10.05.2022	
	void TEST_GETSHAPEFUNCTIONS() {
		// Arrange
		cShapeHex8 element;
		cArray3d N;
		int gaussianPoints = 256;

		// Act
		N = element.getShapeFunctions();
		PetscReal* data = N.getBase();
		/*
		std::cout.precision(8);
		for (int i = 0; i < gaussianPoints; i++)
			std::cout << data[i] << ",\n";
		std::cout << std::endl;
		*/
		
		// Assert
		double refData[256] = { 0.0094373878, 0.035220811, 0.13144586, 0.035220811, 0.035220811, 0.13144586, 0.49056261, 0.13144586,
								-0.022329099, 0.022329099, 0.083333333, -0.083333333, -0.083333333, 0.083333333, 0.31100423, -0.31100423,
								-0.022329099, -0.083333333, 0.083333333, 0.022329099, -0.083333333, -0.31100423, 0.31100423, 0.083333333,
								-0.022329099, -0.083333333, -0.31100423, -0.083333333, 0.022329099, 0.083333333, 0.31100423, 0.083333333,
								0.035220811, 0.0094373878, 0.035220811, 0.13144586, 0.13144586, 0.035220811, 0.13144586, 0.49056261,
								-0.022329099, 0.022329099, 0.083333333, -0.083333333, -0.083333333, 0.083333333, 0.31100423, -0.31100423,
								-0.083333333, -0.022329099, 0.022329099, 0.083333333, -0.31100423, -0.083333333, 0.083333333, 0.31100423,
								-0.083333333, -0.022329099, -0.083333333, -0.31100423, 0.083333333, 0.022329099, 0.083333333, 0.31100423,
								0.035220811, 0.13144586, 0.035220811, 0.0094373878, 0.13144586, 0.49056261, 0.13144586, 0.035220811,
								-0.083333333, 0.083333333, 0.022329099, -0.022329099, -0.31100423, 0.31100423, 0.083333333, -0.083333333,
								-0.022329099, -0.083333333, 0.083333333, 0.022329099, -0.083333333, -0.31100423, 0.31100423, 0.083333333,
								-0.083333333, -0.31100423, -0.083333333, -0.022329099, 0.083333333, 0.31100423, 0.083333333, 0.022329099,
								0.13144586, 0.035220811, 0.0094373878, 0.035220811, 0.49056261, 0.13144586, 0.035220811, 0.13144586,
								-0.083333333, 0.083333333, 0.022329099, -0.022329099, -0.31100423, 0.31100423, 0.083333333, -0.083333333,
								-0.083333333, -0.022329099, 0.022329099, 0.083333333, -0.31100423, -0.083333333, 0.083333333, 0.31100423,
								-0.31100423, -0.083333333, -0.022329099, -0.083333333, 0.31100423, 0.083333333, 0.022329099, 0.083333333,
								0.035220811, 0.13144586, 0.49056261, 0.13144586, 0.0094373878, 0.035220811, 0.13144586, 0.035220811,
								-0.083333333, 0.083333333, 0.31100423, -0.31100423, -0.022329099, 0.022329099, 0.083333333, -0.083333333,
								-0.083333333, -0.31100423, 0.31100423, 0.083333333, -0.022329099, -0.083333333, 0.083333333, 0.022329099,
								-0.022329099, -0.083333333, -0.31100423, -0.083333333, 0.022329099, 0.083333333, 0.31100423, 0.083333333,
								0.13144586, 0.035220811, 0.13144586, 0.49056261, 0.035220811, 0.0094373878, 0.035220811, 0.13144586, -0.083333333,
								0.083333333, 0.31100423, -0.31100423, -0.022329099, 0.022329099, 0.083333333, -0.083333333, -0.31100423,
								-0.083333333, 0.083333333, 0.31100423, -0.083333333, -0.022329099, 0.022329099, 0.083333333, -0.083333333,
								-0.022329099, -0.083333333, -0.31100423, 0.083333333, 0.022329099, 0.083333333, 0.31100423, 0.13144586,
								0.49056261, 0.13144586, 0.035220811, 0.035220811, 0.13144586, 0.035220811, 0.0094373878, -0.31100423,
								0.31100423, 0.083333333, -0.083333333, -0.083333333, 0.083333333, 0.022329099, -0.022329099, -0.083333333,
								-0.31100423, 0.31100423, 0.083333333, -0.022329099, -0.083333333, 0.083333333, 0.022329099, -0.083333333,
								-0.31100423, -0.083333333, -0.022329099, 0.083333333, 0.31100423, 0.083333333, 0.022329099, 0.49056261,
								0.13144586, 0.035220811, 0.13144586, 0.13144586, 0.035220811, 0.0094373878, 0.035220811, -0.31100423,
								0.31100423, 0.083333333, -0.083333333, -0.083333333, 0.083333333, 0.022329099, -0.022329099, -0.31100423,
								-0.083333333, 0.083333333, 0.31100423, -0.083333333, -0.022329099, 0.022329099, 0.083333333, -0.31100423,
								-0.083333333, -0.022329099, -0.083333333, 0.31100423, 0.083333333, 0.022329099, 0.083333333 };
		for (size_t i = 0; i < gaussianPoints; i++) {
			//std::cout << i << "\t" << data[i] << "\t" << refData[i] << "\n";
			ASSERT_NEAR(data[i], refData[i], 1e-6);
		}
		
	}

	
	//! @brief Function to arrange, act, assert - do get Shape of the element
	//! @author Akash Doddamane Lingaraja
	//! @date 10.05.2022	
	void TEST_GETELEMENTSHAPE() {
		// Arrange and Act 
		int shape = cShapeHex8::getElementShape();

		// Assert
		ASSERT_EQ(shape, Hexahedron);
	}

	
	//! @brief Function to arrange, act, assert - do get number of faces of Tria3 elements
	//! @author Akash Doddamane Lingaraja
	//! @date 10.05.2022	
	void TEST_GETNUMBEROFFACES() {
		// Arrange and Act 
		short faces = cShapeHex8::getNumberOfFaces();

		// Assert
		ASSERT_EQ(faces, 6);
	}

	
	//! @brief Function to arrange, act, assert - do get Numbering of one element's face
	//! @author Akash Doddamane Lingaraja
	//! @date 10.05.2022	
	void TEST_GETINDICESOFFACENODES() {
		// Arrange
		cShapeHex8 element;
		std::vector<short> nodesFace0{ 0, 1, 5, 4 };
		std::vector<short> nodesFace1{ 2, 3, 7, 6 };
		std::vector<short> nodesFace2{ 0, 3, 2, 1 };
		std::vector<short> nodesFace3{ 1, 2, 6, 5 };
		std::vector<short> nodesFace4{ 4, 5, 6, 7 };
		std::vector<short> nodesFace5{ 0, 4, 7, 3 };

		// Assert
		ASSERT_EQ(element.getIndicesOfFaceNodes(0), nodesFace0);
		ASSERT_EQ(element.getIndicesOfFaceNodes(1), nodesFace1);
		ASSERT_EQ(element.getIndicesOfFaceNodes(2), nodesFace2);
		ASSERT_EQ(element.getIndicesOfFaceNodes(3), nodesFace3);
		ASSERT_EQ(element.getIndicesOfFaceNodes(4), nodesFace4);
		ASSERT_EQ(element.getIndicesOfFaceNodes(5), nodesFace5);
	}

	
	//! @brief Function to arrange, act, assert - do get number of nodes that describe one face of the element
	//! @author Akash Doddamane Lingaraja
	//! @date 10.05.2022	
	void TEST_GETNUMBEROFNODESPERFACE() {
		// Arrange
		short nodes = 4;
		// Assert
		ASSERT_EQ(getNumberOfNodesPerFace(), nodes);
	}

	
	//! @brief Function to arrange, act, assert - do get Shape function of a face
	//! @author Akash Doddamane Lingaraja
	//! @date 10.05.2022	
	void TEST_GETSHAPEFUNCTIONSFACE() {
		// Arrange
		cArray3d Nface;
		int gaussianPoints = 48;

		// Act
		Nface = getShapeFunctionsFace();
		PetscReal* data = Nface.getBase();
		/*
		std::cout.precision(8);
		int count = 0;
		for (int i = 0; i < gaussianPoints; i++) {
			count++;
			if (count == 5) {
				std::cout << "\n";
				count = 0;
			}
			std::cout << data[i] << ",\t";
		}
		std::cout << std::endl;
		*/
		// Assert
		double refData[48] = { 0.044658199,	0.16666667,	0.62200847,	0.16666667,	-0.10566243, 0.10566243, 0.39433757,
								-0.39433757, -0.10566243, -0.39433757,	0.39433757,	0.10566243,	0.16666667,	0.62200847,
								0.16666667,	0.044658199, -0.39433757, 0.39433757, 0.10566243, -0.10566243, -0.10566243,
								-0.39433757, 0.39433757, 0.10566243, 0.16666667, 0.044658199, 0.16666667, 0.62200847,
								-0.10566243, 0.10566243, 0.39433757, -0.39433757, -0.39433757, -0.10566243, 0.10566243,
								0.39433757,	0.62200847,	0.16666667,	0.044658199, 0.16666667, -0.39433757, 0.39433757,
								0.10566243,	-0.10566243, -0.39433757, -0.10566243, 0.10566243,	0.39433757 };
		for (size_t i = 0; i < gaussianPoints; i++) {
			//std::cout << i << "\t" << data[i] << "\t" << refData[i] << "\n";
			ASSERT_NEAR(data[i], refData[i], 1e-6);
		}
		
	}
};

TEST_F(cTestShapeHex8, checkForCorrectShapeFunctionValues) {
	TEST_GETSHAPEFUNCTIONS();
}

TEST_F(cTestShapeHex8, checkForCorrectShapeOfTheElement) {
	TEST_GETELEMENTSHAPE();
}

TEST_F(cTestShapeHex8, checkForCorrectNumberOfFacesOfTheHex8Element) {
	TEST_GETNUMBEROFFACES();
}

TEST_F(cTestShapeHex8, checkForCorrectNumberingOfFacesOfTheHex8Element) {
	TEST_GETINDICESOFFACENODES();
}

TEST_F(cTestShapeHex8, checkForCorrectNumberOfNodesPerFaceOfTheHex8Element) {
	TEST_GETNUMBEROFNODESPERFACE();
}

TEST_F(cTestShapeHex8, checkForCorrectShapeFunctionValuesOfElementsFace) {
	TEST_GETSHAPEFUNCTIONSFACE();
}