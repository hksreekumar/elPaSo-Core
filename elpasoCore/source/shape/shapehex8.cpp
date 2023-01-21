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

#include "shapehex8.h"

/*BEGIN_NO_COVERAGE*/
cShapeHex8::cShapeHex8()
{
  if (howMany() == 1)
    initializeShapeFunctions();
}


cShapeHex8::cShapeHex8(const cShapeHex8 &other)
{
  // empty
}


cShapeHex8::~cShapeHex8()
{
  // empty
}
/*END_NO_COVERAGE*/

// ---------------------------------------------------------------------------
//  initialize static class members
// ---------------------------------------------------------------------------
cArray3d cShapeHex8::N;
cArray3d cShapeHex8::Nface;


void cShapeHex8::initializeShapeFunctions(void)
{
  const int nnod     = 8;
  const int nnodface = 4;
  const int ngp      = 2;

  // -------------------------------------------------------------------------
  //  allocate memory
  // -------------------------------------------------------------------------
  N.initialize(4, nnod, ngp*ngp*ngp);
  Nface.initialize(3, nnodface, ngp*ngp);


  // -------------------------------------------------------------------------
  //  compute shapefunctions
  // -------------------------------------------------------------------------
  for (int n=0; n<ngp*ngp*ngp; n++)
  {
    cPoint gp( m_GaussPoints.getGaussPoint3D(ngp, n) );

    for (int k=0; k<nnod; k++)
    {
      N(0,k,n) = m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nnod, k, N_fun , gp);
      N(1,k,n) = m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nnod, k, N_xi  , gp);
      N(2,k,n) = m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nnod, k, N_eta , gp);
      N(3,k,n) = m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nnod, k, N_zeta, gp);
    }
  }


  for (int n=0; n<ngp*ngp; n++)
  {
    cPoint gp( m_GaussPoints.getGaussPoint2D(ngp, n) );

    for (int k=0; k<nnodface; k++)
    {
      Nface(0,k,n) = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nnodface, k, N_fun , gp);
      Nface(1,k,n) = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nnodface, k, N_xi  , gp);
      Nface(2,k,n) = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nnodface, k, N_eta , gp);
    }
  }
}


// ----------------------------------------------------------------------------
//   Definition of the 6 surface elements. The numbering of the nodes
//   follows this scheme:
//
//        [3]         [2]
//           o-------o
//           |  <--  |
//           |    |  |
//           |  ---  |
//           o-------o
//        [0]         [1]
//
//   Die Reihenfolge der Oberflaechenelemente entspricht der Konvention
//   von MSC Patran. Dort heissen diese Oberflaechen Faces. Wenn nun also
//   im Eingabedatensatz von Face Nr. 3 die Rede ist, so entspricht dies
//   hier surfaces[3-1][..]. Die Normalenvektoren zeigen aus dem Volumen-
//   element heraus.
//           o-------o
//          /       /|
//         /       / |
//        o-------o ---> n
//        |       |  o
//        |       | /
//        |       |/
//        o-------o
// ----------------------------------------------------------------------------
std::vector<short> cShapeHex8::getIndicesOfFaceNodes(int Face) const
{
  std::vector<short> surface(4);

  switch(Face)
  {
    case 0:  // n = [  0; -1;  0 ]
      surface[0] = 0; surface[1] = 1; surface[2] = 5; surface[3] = 4;
      break;
    case 1:  // n = [  0;  1;  0 ]
      surface[0] = 2; surface[1] = 3; surface[2] = 7; surface[3] = 6;
      break;
    case 2:  // n = [  0;  0; -1 ]
      surface[0] = 0; surface[1] = 3; surface[2] = 2; surface[3] = 1;
      break;
    case 3:  // n = [  1;  0;  0 ]
      surface[0] = 1; surface[1] = 2; surface[2] = 6; surface[3] = 5;
      break;
    case 4:  // n = [  0;  0;  1 ]
      surface[0] = 4; surface[1] = 5; surface[2] = 6; surface[3] = 7;
      break;
    case 5:  // n = [ -1;  0;  0 ]
      surface[0] = 0; surface[1] = 4; surface[2] = 7; surface[3] = 3;
      break;
  }

  return surface;
}

