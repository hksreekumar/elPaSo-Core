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

#ifndef INFAM_ENTITIES_H
#define INFAM_ENTITIES_H

#include <iostream>


#include "group.h"

// ---------------------------------------------------------------------------
//    analysistypes
// ---------------------------------------------------------------------------
#include "../analysis/frequency/analysisfrequencybasic.h"

// ---------------------------------------------------------------------------
//    interfaces
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
//    elementtypes
// ---------------------------------------------------------------------------

//****************************************************************************************
//                             StructureLinear
//****************************************************************************************

#include "../element/structure/linear/mass/elementstructuremass.h"

#include "../element/structure/linear/beam/elementstructurebeam.h"
#include "../element/structure/linear/beam/elementstructurebeam3d.h"
#include "../element/structure/linear/beam/elementstructurebeam10.h"
#include "../element/structure/linear/beam/elementstructurebeam12.h"

#include "../element/structure/linear/spring/elementstructurespring.h"
#include "../element/structure/linear/spring/elementstructurespringz.h"
#include "../element/structure/linear/spring/elementstructurespringbc.h"
#include "../element/structure/linear/spring/elementstructurespringbcx.h"
#include "../element/structure/linear/spring/elementstructurespringbcy.h"
#include "../element/structure/linear/spring/elementstructurespringbcz.h"
#include "../element/structure/linear/spring/elementstructurespringbcrx.h"
#include "../element/structure/linear/spring/elementstructurespringbcry.h"
#include "../element/structure/linear/spring/elementstructurespringbcrz.h"

#include "../element/structure/linear/brick/elementstructurebrick8.h"

#include "../element/structure/linear/plate/elementstructuremindlindsg4.h"
#include "../element/structure/linear/plate/elementstructurekirchhoff4.h"

//****************************************************************************************
//                             StructureNonLinear
//****************************************************************************************

//****************************************************************************************
//                             StructurePoro
//****************************************************************************************

//****************************************************************************************
//                             StructurePoro3d
//****************************************************************************************


//****************************************************************************************
//                             Fluid Flow auxiliary elements
//****************************************************************************************


//****************************************************************************************
//                             FluidLinear
//****************************************************************************************
#include "../element/fluid/elementfluid8.h"

//****************************************************************************************
#include "../element/interface/elementinterface.h"
#include "../element/interface/elementinterfacemindlin.h"

#include "../element/elementcoupling.h"


//***********************************************************************
//               Non Conforming
//***********************************************************************

// ---------------------------------------------------------------------------
//    materialdefinitions
// ---------------------------------------------------------------------------
#include "../material/structure/materialstructure.h"
#include "../material/structure/isotrop/materialstructureiso.h"
#include "../material/structure/isotrop/linear/materialstructureisolin.h"
#include "../material/structure/isotrop/linear/elastic/materialstructureisotrop.h"
#include "../material/structure/isotrop/linear/elastic/materialspring.h"
#include "../material/structure/isotrop/linear/elastic/materialmass.h"
#include "../material/structure/isotrop/linear/viscoelastic/materialstructurevisco.h"

#include "../material/fluid/materialfluid.h"
#include "../material/fluid/linear/materialfluidlin.h"
#include "../material/fluid/linear/elastic/materialfluidideal.h"

// ---------------------------------------------------------------------------
//    elementloads
// ---------------------------------------------------------------------------
#include "../element/load/elementloadstructureconst.h"
#include "../element/load/elementloadfluid.h"

// ---------------------------------------------------------------------------
//    nodal boundaryconditions
// ---------------------------------------------------------------------------
#include "../bc/boundaryconditionstructure.h"
#include "../bc/boundaryconditionfluid.h"

// ---------------------------------------------------------------------------
//    nodal forces
// ---------------------------------------------------------------------------
#include "../misc/nodalforce/nodalforcestructure.h"
#include "../misc/nodalforce/nodalforcefluid.h"

// ---------------------------------------------------------------------------
//    nodal moments
// ---------------------------------------------------------------------------
#include "../misc/nodalmoment/nodalmomentstructure.h"

// ----------------------------------------------------------------------------
//     nodal pressures
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
//     nodal values (fluid flow)
// ----------------------------------------------------------------------------

// ---------------------------------------------------------------------------
//    constraints
// ---------------------------------------------------------------------------

void dumpObjectCounters(FILE *os);

#endif
