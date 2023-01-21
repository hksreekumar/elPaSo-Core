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

/*--------------------------------------------------------------------------
* elPaSo - Test project for elpasoCore
*
* 16.10.2020
* Harikrishnan Sreekumar
* Institut für Akustik, Technische Universität Braunschweig
*---------------------------------------------------------------------------*/

#ifdef HAVE_GTEST
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#endif // HAVE_GTEST

#ifdef HAVE_HDF5
#include "./misc/hdf5/testHandlerHDF5.h"
#include "./misc/hdf5/testReaderHDF5.h"
#include "./misc/hdf5/testInputHDF5.h"
#include "./misc/hdf5/testOutputHDF5.h"
#endif // HAVE_HDF5
#ifdef HAVE_HDF5
#include "./misc/parser/testFemParserHDF5.h"
#endif // HAVE_HDF5

// Analysis

// FEDATA
#include "./fedata/testdof.h"
#include "./fedata/testproblem.h"

// Mathlibrary

// Material
#include "./material/structure/isotrop/linear/viscoelastic/testmaterialstructurevisco.h"
#include "./material/structure/isotrop/linear/elastic/testmaterialstructureisotrop.h"
#include "./material/fluid/linear/elastic/testmaterialfluidideal.h"
#include "./material/structure/isotrop/linear/elastic/testmaterialspring.h"
#include "./material/structure/testmaterialstructure.h"
#include "./material/fluid/testmaterialfluid.h"
#include "./material/testmaterial.h"

// Misc
#include "./misc/testpoint.h"

// Boundary Condition
#include "./bc/testboundaryconditionstructure.h"
#include "./bc/testboundaryconditionfluid.h"

// Basics
#include "./basics/watch/teststopwatch.h"

// Shape
#include "./shape/testshapefunction.h"
#include "./shape/testshapequad4.h"
#include "./shape/testshapehex8.h"

//! @file testmain.cpp
//! @brief Starting point for elPaSo tests
//! @author Harikrishnan Sreekumar
//! @date 16.10.2020
int main(int argc, char* argv[])
{
	PetscInitialize(NULL, NULL, (char*)0, NULL);
	SlepcInitialize(NULL, NULL, (char*)0, NULL);

	::testing::InitGoogleMock(&argc, argv);
	int err_Code = RUN_ALL_TESTS();

	SlepcFinalize();
	PetscFinalize();
	return err_Code;
}