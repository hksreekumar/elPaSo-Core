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

#include "shapequad4.h"

/*BEGIN_NO_COVERAGE*/
cShapeQuad4::cShapeQuad4()
{
  if (howMany() == 1)
    initializeShapeFunctions();
}


cShapeQuad4::cShapeQuad4(const cShapeQuad4 &other)
{
}


cShapeQuad4::~cShapeQuad4()
{
}
/*END_NO_COVERAGE*/

// ---------------------------------------------------------------------------
//  initialize static class members
// ---------------------------------------------------------------------------
cArray3d cShapeQuad4::N;
cArray3d cShapeQuad4::N_map;



void cShapeQuad4::initializeShapeFunctions(void)
{
  const int nnod = 4;
  const int ngp  = 2;

  // -------------------------------------------------------------------------
  //  allocate memory
  // -------------------------------------------------------------------------
  N.initialize(3, nnod, ngp*ngp);


  // -------------------------------------------------------------------------
  //  evaluate shape functions at Gauss points
  // -------------------------------------------------------------------------
  for (int n=0; n<ngp*ngp; n++)
  {
    cPoint gp( m_GaussPoints.getGaussPoint2D(ngp, n) );

    for (int k=0; k<nnod; k++)
    {
      N(0,k,n) = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nnod, k, N_fun , gp);
      N(1,k,n) = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nnod, k, N_xi  , gp);
      N(2,k,n) = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nnod, k, N_eta , gp);
    }
  }
}

void cShapeQuad4::initializeShapeFunctions_N_map(void) //KR
{
  const int nnod = 4;
  const int ngp  = 2;

  // -------------------------------------------------------------------------
  //  allocate memory
  // -------------------------------------------------------------------------
  N_map.initialize(3, nnod, ngp*ngp);


  // -------------------------------------------------------------------------
  //  evaluate shape functions at Gauss points
  // -------------------------------------------------------------------------
  for (int n=0; n<ngp*ngp; n++)
  {
    cPoint gp( m_GaussPoints.getGaussPoint2D(ngp, n) );

    for (int k=0; k<nnod; k++)
    {
      N_map(0,k,n) = 0.;
      N_map(1,k,n) = 0.;
      N_map(2,k,n) = 0.;
    }
  }
}

// ----------------------------------------------------------------------------
//
//        [3]         [2]
//           o-------o
//           |       |
//           |       |--> n
//           |       |
//           o-------o
//        [0]         [1]
//
// ----------------------------------------------------------------------------
std::vector<short> cShapeQuad4::getIndicesOfFaceNodes(int Face) const
{
  std::vector<short> surface(2);

  switch(Face)
  {
    case 0:  // n = [  0; -1;  0 ]
      surface[0] = 1; surface[1] = 0;
      break;
    case 1:  // n = [  1;  0;  0 ]
      surface[0] = 2; surface[1] = 1;
      break;
    case 2:  // n = [  0;  1;  0 ]
      surface[0] = 3; surface[1] = 2;
      break;
    case 3:  // n = [ -1;  0;  0 ]
      surface[0] = 0; surface[1] = 3;
      break;
    default:
      message("*** invalid value for face : %d\n", Face);
      throw cException("Abort", __FILE__, __LINE__);
  }

  return surface;
}

