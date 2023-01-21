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

#include "shapefunctions.h"

/*BEGIN_NO_COVERAGE*/
cShapeFunctions::cShapeFunctions()
{
  // empty
}


cShapeFunctions::~cShapeFunctions()
{
  // empty
}
/*END_NO_COVERAGE*/

double cShapeFunctions::evaluateShapeFunction(
  const eElementShape &shape,
  const int &nnod,
  const int &k,
  const eDerivative &deriv,
  const cPoint &gp)
{
  switch (shape)
  {
    case Beam:
      return BeamN(nnod, k, deriv, gp);
      break;
    case Tria:
      return TriaN(nnod, k, deriv, gp);
      break;
    case Quadrilateral:
      return LagrangeQuadN(nnod, k, deriv, gp);
      break;
    case QuadSerendipity:
      return SerendipityQuadN(nnod, k, deriv, gp);
      break;
    case Hexahedron:
      return LagrangeHexN(nnod, k, deriv, gp);
      break;
    case Tetrahedron:
      return TetN(nnod, k, deriv, gp);
      break;
    default:
      break;
  }

  return 0.0;
}



// ---------------------------------------------------------------------------
//   test functions of quadrilateral elements (4 and 9 nodes)
// ---------------------------------------------------------------------------
double cShapeFunctions::LagrangeQuadN(
  const int &nnod, 
  const int &k, 
  const eDerivative &deriv,
  const cPoint &gp)
{
  switch (nnod)
  {
    case 4:
      return LagrangeQuad4(k, deriv, gp);
      break;
    case 9:
      return LagrangeQuad9(k, deriv, gp);
      break;
    default:
      break;
  }

  return 0.0;
}


// ----------------------------------------------------------------------------
//
//             ^ eta
//       3     |     2
//        o----|----o
//        |    |    |
//        |    |    |
//       ------+-------->*xi
//        |    |    |
//        |    |    |
//        o----|----o
//       0     |     1
//
// ----------------------------------------------------------------------------
double cShapeFunctions::LagrangeQuad4(
  const int &k,
  const eDerivative &deriv,
  const cPoint &gp)
{
  double fkt = 0.0;

  switch(deriv)
  {
    case N_fun:
      switch(k)
      {
        case 0:
          fkt = 0.25 * (1.0 - gp[0]) * (1.0 - gp[1]);
          break;
        case 1:
          fkt = 0.25 * (1.0 + gp[0]) * (1.0 - gp[1]);
          break;
        case 2:
          fkt = 0.25 * (1.0 + gp[0]) * (1.0 + gp[1]);
          break;
        case 3:
          fkt = 0.25 * (1.0 - gp[0]) * (1.0 + gp[1]);
          break;
      }
      break;
  
    case N_xi:
      switch(k)
      {
        case 0:
          fkt = -0.25 * (1.0 - gp[1]);
          break;
        case 1:
          fkt =  0.25 * (1.0 - gp[1]);
          break;
        case 2:
          fkt =  0.25 * (1.0 + gp[1]);
          break;
        case 3:
          fkt = -0.25 * (1.0 + gp[1]);
          break;
      }
      break;
  
    case N_eta:
      switch(k)
      {
        case 0:
          fkt = -0.25 * (1.0 - gp[0]);
          break;
        case 1:
          fkt = -0.25 * (1.0 + gp[0]);
          break;
        case 2:
          fkt =  0.25 * (1.0 + gp[0]);
          break;
        case 3:
          fkt =  0.25 * (1.0 - gp[0]);
          break;
      }
      break;

    default:
      throw cException(" **** invalid request for derivative **** ", __FILE__, __LINE__);
      break;
  }

  return fkt;
}


// ----------------------------------------------------------------------------
//
//             ^ eta
//       3     |6    2
//        o----o----o
//        |    |    |
//       7|    |8   |5
//       -o----o----o--->*xi
//        |    |    |
//        |    |    |
//        o----o----o
//       0     |4    1
//
// ----------------------------------------------------------------------------
double cShapeFunctions::LagrangeQuad9(
  const int &k,
  const eDerivative &deriv,
  const cPoint &gp)
{
  double fkt = 0.0;

  switch(deriv)
  {
    case N_fun:
      switch(k)
      {
        case 0:
          fkt =  0.25 * gp[0]*(1.0 - gp[0]) * gp[1]*(1.0 - gp[1]);
          break;
        case 1:
          fkt = -0.25 * gp[0]*(1.0 + gp[0]) * gp[1]*(1.0 - gp[1]);
          break;
        case 2:
          fkt =  0.25 * gp[0]*(1.0 + gp[0]) * gp[1]*(1.0 + gp[1]);
          break;
        case 3:
          fkt = -0.25 * gp[0]*(1.0 - gp[0]) * gp[1]*(1.0 + gp[1]);
          break;
        case 4:
          fkt = -0.50 * (1.0-gp[0]*gp[0]) * gp[1]*(1.0 - gp[1]);;
          break;
        case 5:
          fkt =  0.50 * gp[0]*(1.0+gp[0]) * (1.0 - gp[1]*gp[1]);
          break;
        case 6:
          fkt =  0.50 * (1.0-gp[0]*gp[0]) * gp[1]*(1.0 + gp[1]);
          break;
        case 7:
          fkt = -0.50 * gp[0]*(1.0 - gp[0]) * (1.0 - gp[1]*gp[1]);
          break;
        case 8:
          fkt = (1.0-gp[0]*gp[0])*(1.0-gp[1]*gp[1]);
          break;
      }
      break;


    case N_xi:
      switch(k)
      {
        case 0:
          fkt =  0.25 * (1.0 - 2.0*gp[0]) * gp[1]*(1.0 - gp[1]);
          break;
        case 1:
          fkt = -0.25 * (1.0 + 2.0*gp[0]) * gp[1]*(1.0 - gp[1]);
          break;
        case 2:
          fkt =  0.25 * (1.0 + 2.0*gp[0]) * gp[1]*(1.0 + gp[1]);
          break;
        case 3:
          fkt = -0.25 * (1.0 - 2.0*gp[0]) * gp[1]*(1.0 + gp[1]);
          break;
        case 4:
          fkt = -0.50 * (-2.0*gp[0]) * gp[1]*(1.0 - gp[1]);
          break;
        case 5:
          fkt =  0.50 * (1.0+2.0*gp[0]) * (1.0-gp[1]*gp[1]);
          break;
        case 6:
          fkt =  0.50 * (-2.0*gp[0]) * gp[1]*(1.0 + gp[1]);
          break;
        case 7:
          fkt = -0.50 * (1.0 - 2.0*gp[0]) * (1.0 - gp[1]*gp[1]);
          break;
        case 8:
          fkt = -2.0*gp[0] * (1.0 - gp[1]*gp[1]);
          break;
      }
      break;


    case N_eta:
      switch(k)
      {
        case 0:
          fkt =  0.25 * gp[0]*(1.0 - gp[0]) * (1.0 - 2.0*gp[1]);
          break;
        case 1:
          fkt = -0.25 * gp[0]*(1.0 + gp[0]) * (1.0 - 2.0*gp[1]);
          break;
        case 2:
          fkt =  0.25 * gp[0]*(1.0 + gp[0]) * (1.0 + 2.0*gp[1]);
          break;
        case 3:
          fkt = -0.25 * gp[0]*(1.0 - gp[0]) * (1.0 + 2.0*gp[1]);
          break;
        case 4:
          fkt = -0.50 * (1.0 - gp[0]*gp[0]) * (1.0 - 2.0*gp[1]);
          break;
        case 5:
          fkt =  0.50 * gp[0]*(1.0 + gp[0])*(-2.0*gp[1]);
          break;
        case 6:
          fkt =  0.50 * (1.0 - gp[0]*gp[0]) * (1.0 + 2.0*gp[1]);
          break;
        case 7:
          fkt = -0.50 * gp[0]*(1.0 - gp[0]) * (-2.0 * gp[1]);
          break;
        case 8:
          fkt = (1.0 - gp[0]*gp[0]) * (-2.0*gp[1]);
          break;
      }
      break;


    default:
      throw cException(" **** invalid request for derivative **** ", __FILE__, __LINE__);
      break;
  }

  return fkt;
}

// ---------------------------------------------------------------------------
//   test functions of quadrilateral elements (4 and 9 nodes)
// ---------------------------------------------------------------------------
double cShapeFunctions::HermiteQuadN(
  const int &nnod, 
  const int &k, 
  const eDerivative &deriv,
  const cPoint &gp)
{
  switch (nnod)
  {
    case 4:
      return HermiteQuad4(k, deriv, gp);
      break;
//    case 9:
//      return LagrangeQuad9(k, deriv, gp);
//      break;
    default:
      break;
  }

  return 0.0;
}

// ----------------------------------------------------------------------------
//
//             ^ eta
//       3     |     2
//        o----|----o
//        |    |    |
//        |    |    |
//       ------+-------->*xi
//        |    |    |
//        |    |    |
//        o----|----o
//       0     |     1
//
// ----------------------------------------------------------------------------
double cShapeFunctions::HermiteQuad4(
  const int &k,
  const eDerivative &deriv,
  const cPoint &gp)
{
  double fkt = 0.0;

  switch(deriv)
  {
    case N_fun:
      switch(k)
      {
        case 0:
          fkt = 1./16. * (2.0 - 3.0*gp[0] + gp[0]*gp[0]*gp[0]) * (2.0 - 3.0*gp[1] + gp[1]*gp[1]*gp[1]);
          break;
        case 1:
          fkt = 1./16. * (2.0 + 3.0*gp[0] - gp[0]*gp[0]*gp[0]) * (2.0 - 3.0*gp[1] + gp[1]*gp[1]*gp[1]);
          break;
        case 2:
          fkt = 1./16. * (2.0 + 3.0*gp[0] - gp[0]*gp[0]*gp[0]) * (2.0 + 3.0*gp[1] - gp[1]*gp[1]*gp[1]);
          break;
        case 3:
          fkt = 1./16. * (2.0 - 3.0*gp[0] + gp[0]*gp[0]*gp[0]) * (2.0 + 3.0*gp[1] - gp[1]*gp[1]*gp[1]);
          break;
      }
      break;
  
    case N_xi:
      switch(k)
      {
        case 0:
          fkt = 1./16. * (-3.0 + 3.0*gp[0]*gp[0]) * (2.0 - 3.0*gp[1] + gp[1]*gp[1]*gp[1]);
          break;
        case 1:
          fkt = 1./16. * ( 3.0 + 3.0*gp[0]*gp[0]) * (2.0 - 3.0*gp[1] + gp[1]*gp[1]*gp[1]);
          break;
        case 2:
          fkt = 1./16. * ( 3.0 - 3.0*gp[0]*gp[0]) * (2.0 + 3.0*gp[1] - gp[1]*gp[1]*gp[1]);
          break;
        case 3:
          fkt = 1./16. * (-3.0 + 3.0*gp[0]*gp[0]) * (2.0 + 3.0*gp[1] - gp[1]*gp[1]*gp[1]);
          break;
      }
      break;
  
    case N_eta:
      switch(k)
      {
        case 0:
          fkt = 1./16. * (2.0 - 3.0*gp[0] + gp[0]*gp[0]*gp[0]) * (-3.0 + 3.0*gp[1]*gp[1]);
          break;
        case 1:
          fkt = 1./16. * (2.0 + 3.0*gp[0] - gp[0]*gp[0]*gp[0]) * (-3.0 + 3.0*gp[1]*gp[1]);
          break;
        case 2:
          fkt = 1./16. * (2.0 + 3.0*gp[0] - gp[0]*gp[0]*gp[0]) * ( 3.0 - 3.0*gp[1]*gp[1]);
          break;
        case 3:
          fkt = 1./16. * (2.0 - 3.0*gp[0] + gp[0]*gp[0]*gp[0]) * ( 3.0 - 3.0*gp[1]*gp[1]);
          break;
      }
      break;

    case N_xi2:
      switch(k)
      {
        case 0:
          fkt = 1./16. * ( 6.0*gp[0]) * (2.0 - 3.0*gp[1] + gp[1]*gp[1]*gp[1]);
          break;
        case 1:
          fkt = 1./16. * ( 6.0*gp[0]) * (2.0 - 3.0*gp[1] + gp[1]*gp[1]*gp[1]);
          break;
        case 2:
          fkt = 1./16. * (-6.0*gp[0]) * (2.0 + 3.0*gp[1] - gp[1]*gp[1]*gp[1]);
          break;
        case 3:
          fkt = 1./16. * (-6.0*gp[0]) * (2.0 + 3.0*gp[1] - gp[1]*gp[1]*gp[1]);
          break;
      }
      break;
  
    case N_eta2:
      switch(k)
      {
        case 0:
          fkt = 1./16. * (2.0 - 3.0*gp[0] + gp[0]*gp[0]*gp[0]) * ( 6.0*gp[1]);
          break;
        case 1:
          fkt = 1./16. * (2.0 + 3.0*gp[0] - gp[0]*gp[0]*gp[0]) * ( 6.0*gp[1]);
          break;
        case 2:
          fkt = 1./16. * (2.0 + 3.0*gp[0] - gp[0]*gp[0]*gp[0]) * (-6.0*gp[1]);
          break;
        case 3:
          fkt = 1./16. * (2.0 - 3.0*gp[0] + gp[0]*gp[0]*gp[0]) * (-6.0*gp[1]);
          break;
      }
      break;

    default:
      throw cException(" **** invalid request for derivative **** ", __FILE__, __LINE__);
      break;
  }

  return fkt;
}



// ---------------------------------------------------------------------------
//   serendipity shape functions
// ---------------------------------------------------------------------------
PetscReal cShapeFunctions::SerendipityQuadN(
  const int &nnod, 
  const int &k, 
  const eDerivative &deriv, 
  const cPoint &gp)
{
  switch (nnod)
  {
    case 8:
      return SerendipityQuad8(k, deriv, gp);
      break;
    default:
      break;
  }

  return 0.0;
}


// ----------------------------------------------------------------------------
//
//             ^ eta
//       3     |6    2
//        o----o----o
//        |    |    |
//       7|    |    |5
//       -o----+----o--->*xi
//        |    |    |
//        |    |    |
//        o----o----o
//       0     |4    1
//
// ----------------------------------------------------------------------------
PetscReal cShapeFunctions::SerendipityQuad8(
  const int &k, 
  const eDerivative &deriv, 
  const cPoint &gp)
{
  double fkt = 0.0;

  double xi = gp[0], eta = gp[1];

  switch(deriv)
  {
    case N_fun:
      switch(k)
      {
      case 0:
        fkt = 0.25 * (1. -xi) * (1. - eta) * ( -xi - eta - 1.);
        break;
      case 1:
        fkt = 0.25 * (1. +xi) * (1. - eta) * ( xi - eta - 1.);
        break;
      case 2:
        fkt = 0.25 * (1. +xi) * (1. + eta) * ( xi + eta - 1.);
        break;
      case 3:
        fkt = 0.25 * (1. -xi) * (1. + eta) * ( -xi + eta - 1.);
        break;
      case 4:
        fkt = 0.50 * (1. -xi *xi) * (1. - eta);
        break;
      case 5:
        fkt = 0.50 * (1. +xi) * (1. - eta * eta);
        break;
      case 6:
        fkt = 0.50 * (1. -xi *xi) * (1. + eta);
        break;
      case 7:
        fkt = 0.50 * (1. -xi) * (1. - eta * eta);
        break;
      }
      break;


    case N_xi:
      switch(k)
      {
      case 0:
        fkt = -0.25 * (1. - eta) * (-xi - eta - 1.) - 0.25 * (1. -xi) * (1. - eta);
        break;
      case 1:
        fkt =  0.25 * (1. - eta) * (xi - eta - 1.) + 0.25 * (1. +xi) * (1. - eta);
        break;
      case 2:
        fkt =  0.25 * (1. + eta) * (xi + eta - 1.) + 0.25 * (1. +xi) * (1. + eta);
        break;
      case 3:
        fkt = -0.25 * (1. + eta) * (-xi + eta - 1.) - 0.25 * (1. -xi) * (1. + eta);
        break;
      case 4:
        fkt = -1.00 *xi * (1. - eta);
        break;
      case 5:
        fkt =  0.50 - 0.50 * eta * eta;
        break;
      case 6:
        fkt = -1.00 *xi * (1. + eta);
        break;
      case 7:
        fkt = -0.50 + 0.50 * eta * eta;
        break;
      }
      break;


    case N_eta:
      switch(k)
      {
      case 0:
        fkt = -0.25 * (1. -xi) * (-xi - eta - 1.) - 0.25 * (1. -xi) * (1. - eta);
        break;
      case 1:
        fkt = -0.25 * (1. +xi) * (xi - eta - 1.) - 0.25 * (1. +xi) * (1. - eta);
        break;
      case 2:
        fkt =  0.25 * (1. +xi) * (xi + eta - 1.) + 0.25 * (1. +xi) * (1. + eta);
        break;
      case 3:
        fkt =  0.25 * (1. -xi) * (-xi + eta - 1.) + 0.25 * (1. -xi) * (1. + eta);
        break;
      case 4:
        fkt = -0.50 + 0.50 *xi *xi;
        break;
      case 5:
        fkt = -1.00 * (1. +xi) * eta;
        break;
      case 6:
        fkt =  0.50 - 0.50 *xi *xi;
        break;
      case 7:
        fkt = -1.00 * (1. -xi) * eta;
        break;
      }
      break;

    default:
      throw cException(" **** invalid request for derivative **** ", __FILE__, __LINE__);
      break;
  }

  return fkt;
}


// ---------------------------------------------------------------------------
//   test functions of hexahedrons (8 and 27 nodes)
//   The sequence of the node is corresponding to that one used by
//   MSC Patran
// ---------------------------------------------------------------------------
double cShapeFunctions::LagrangeHexN(
  const int &nnod,
  const int &k,
  const eDerivative &deriv,
  const cPoint &gp)
{
  switch (nnod)
  {
    case 8:
      return LagrangeHex8(k, deriv, gp);
      break;
    case 20:
      return LagrangeHex20(k, deriv, gp);
      break;
    case 27:
      return LagrangeHex27(k, deriv, gp);
      break;
    default:
      break;
  }

  return 0.0;
}


PetscReal cShapeFunctions::LagrangeHex20(
  const int &k,
  const eDerivative &deriv,
  const cPoint &gp)
{
  PetscReal fkt = 0.0;

  switch(deriv)
  {
  // ======================================================================
  //    N(xsi,eta,zeta)
  // ======================================================================
  case N_fun:
    switch(k)
    {
      //corners I
    case 0:
      fkt = (1-gp[0])*(1-gp[1])*(1+gp[2])*(-gp[0]-gp[1]+gp[2]-2)/8;
      break;
    case 1:
      fkt = 0.125*(1-gp[0])*(1+gp[1])*(1+gp[2])*(-gp[0]+gp[1]+gp[2]-2);
      break;
    case 2:
      fkt = 0.125*(1+gp[0])*(1+gp[1])*(1+gp[2])*( gp[0]+gp[1]+gp[2]-2);
      break;
    case 3:
      fkt = 0.125*(1+gp[0])*(1-gp[1])*(1+gp[2])*( gp[0]-gp[1]+gp[2]-2);
      break;
      //corners II
    case 4:
      fkt = 0.125*(1-gp[0])*(1-gp[1])*(1-gp[2])*(-gp[0]-gp[1]-gp[2]-2);
      break;
    case 5:
      fkt = 0.125*(1-gp[0])*(1+gp[1])*(1-gp[2])*(-gp[0]+gp[1]-gp[2]-2);
      break;
    case 6:
      fkt = 0.125*(1+gp[0])*(1+gp[1])*(1-gp[2])*(+gp[0]+gp[1]-gp[2]-2);
      break;
    case 7:
      fkt = 0.125*(1+gp[0])*(1-gp[1])*(1-gp[2])*(+gp[0]-gp[1]-gp[2]-2);
      break;
      //midsides I
    case 9:
      fkt = 0.25*(1-gp[0]*gp[0])*(1+gp[1])*(1+gp[2]);
      break;
    case 11:
      fkt = 0.25*(1-gp[0]*gp[0])*(1-gp[1])*(1+gp[2]);
      break;
    case 17:
      fkt = 0.25*(1-gp[0]*gp[0])*(1+gp[1])*(1-gp[2]);
      break;
    case 19:
      fkt = 0.25*(1-gp[0]*gp[0])*(1-gp[1])*(1-gp[2]);
      break;
      //midsides II
    case 8:
      fkt = 0.25*(1-gp[0])*(1-gp[1]*gp[1])*(1+gp[2]);
      break;
    case 10:
      fkt = 0.25*(1+gp[0])*(1-gp[1]*gp[1])*(1+gp[2]);
      break;
    case 16:
      fkt = 0.25*(1-gp[0])*(1-gp[1]*gp[1])*(1-gp[2]);
      break;
    case 18:
      fkt = 0.25*(1+gp[0])*(1-gp[1]*gp[1])*(1-gp[2]);
      break;
      //midsides III
    case 12:
      fkt = 0.25*(1-gp[0])*(1-gp[1])*(1-gp[2]*gp[2]);
      break;
    case 13:
      fkt = 0.25*(1-gp[0])*(1+gp[1])*(1-gp[2]*gp[2]);
      break;
    case 14:
      fkt = 0.25*(1+gp[0])*(1+gp[1])*(1-gp[2]*gp[2]);
      break;
    case 15:
      fkt = 0.25*(1+gp[0])*(1-gp[1])*(1-gp[2]*gp[2]);
      break;
    }
    break;


  // ======================================================================
  //    N,xsi(xsi,eta,zeta)
  // ======================================================================
  case N_xi:
    switch(k)
    {
      //corners I
    case 0:
      fkt = 0.125*(1-gp[1])*(1+gp[2])*(-(-gp[0]-gp[1]+gp[2]-2) - (1-gp[0]) );
      break;
    case 1:
      fkt = 0.125*(1+gp[1])*(1+gp[2])*(-(-gp[0]+gp[1]+gp[2]-2) - (1-gp[0]) );
      break;
    case 2:
      fkt = 0.125*(1+gp[1])*(1+gp[2])*( (+gp[0]+gp[1]+gp[2]-2) + (1+gp[0]) );
      break;
    case 3:
      fkt = 0.125*(1-gp[1])*(1+gp[2])*( (+gp[0]-gp[1]+gp[2]-2) + (1+gp[0]) );
      break;
      //corners II
    case 4:
      fkt = 0.125*(1-gp[1])*(1-gp[2])*(-(-gp[0]-gp[1]-gp[2]-2) - (1-gp[0]) );
      break;
    case 5:
      fkt = 0.125*(1+gp[1])*(1-gp[2])*(-(-gp[0]+gp[1]-gp[2]-2) - (1-gp[0]) );
      break;
    case 6:
      fkt = 0.125*(1+gp[1])*(1-gp[2])*( ( gp[0]+gp[1]-gp[2]-2) + (1+gp[0]) );
      break;
    case 7:
      fkt = 0.125*(1-gp[1])*(1-gp[2])*( ( gp[0]-gp[1]-gp[2]-2) + (1+gp[0]) );
      break;
      //midsides I
    case 9:
      fkt = -0.5*gp[0]*(1+gp[1])*(1+gp[2]);
      break;
    case 11:
      fkt = -0.5*gp[0]*(1-gp[1])*(1+gp[2]);
      break;
    case 17:
      fkt = -0.5*gp[0]*(1+gp[1])*(1-gp[2]);
      break;
    case 19:
      fkt = -0.5*gp[0]*(1-gp[1])*(1-gp[2]);
      break;
      //midsides II
    case 8:
      fkt = -0.25*(1-gp[1]*gp[1])*(1+gp[2]);
      break;
    case 10:
      fkt =  0.25*(1-gp[1]*gp[1])*(1+gp[2]);
      break;
    case 16:
      fkt = -0.25*(1-gp[1]*gp[1])*(1-gp[2]);
      break;
    case 18:
      fkt =  0.25*(1-gp[1]*gp[1])*(1-gp[2]);
      break;
      //midsides III
    case 12:
      fkt = -0.25*(1-gp[1])*(1-gp[2]*gp[2]);
      break;
    case 13:
      fkt = -0.25*(1+gp[1])*(1-gp[2]*gp[2]);
      break;
    case 14:
      fkt =  0.25*(1+gp[1])*(1-gp[2]*gp[2]);
      break;
    case 15:
      fkt =  0.25*(1-gp[1])*(1-gp[2]*gp[2]);
      break; 
    }
    break;


  // ======================================================================
  //    N,eta(xsi,eta,zeta)
  // ======================================================================
  case N_eta:
    switch(k)
    {
      //corners I
    case 0:
      fkt = 0.125*(1-gp[0])*(1+gp[2])*(-(-gp[0]-gp[1]+gp[2]-2) - (1-gp[1]) );
      break;
    case 1:
      fkt = 0.125*(1-gp[0])*(1+gp[2])*( (-gp[0]+gp[1]+gp[2]-2) + (1+gp[1]) );
      break;
    case 2:
      fkt = 0.125*(1+gp[0])*(1+gp[2])*( ( gp[0]+gp[1]+gp[2]-2) + (1+gp[1]) );
      break;
    case 3:
      fkt = 0.125*(1+gp[0])*(1+gp[2])*(-( gp[0]-gp[1]+gp[2]-2) - (1-gp[1]) );
      break;
      //corners II
    case 4:
      fkt = 0.125*(1-gp[0])*(1-gp[2])*(-(-gp[0]-gp[1]-gp[2]-2) - (1-gp[1]) );
      break;
    case 5:
      fkt = 0.125*(1-gp[0])*(1-gp[2])*( (-gp[0]+gp[1]-gp[2]-2) + (1+gp[1]) );
      break;
    case 6:
      fkt = 0.125*(1+gp[0])*(1-gp[2])*( ( gp[0]+gp[1]-gp[2]-2) + (1+gp[1]) );
      break;
    case 7:
      fkt = 0.125*(1+gp[0])*(1-gp[2])*(-( gp[0]-gp[1]-gp[2]-2) - (1-gp[1]) );
      break;
      //midsides I
    case 9:
      fkt =  0.25*(1-gp[0]*gp[0])*(1+gp[2]);
      break;
    case 11:
      fkt = -0.25*(1-gp[0]*gp[0])*(1+gp[2]);
      break;
    case 17:
      fkt =  0.25*(1-gp[0]*gp[0])*(1-gp[2]);
      break;
    case 19:
      fkt = -0.25*(1-gp[0]*gp[0])*(1-gp[2]);
      break;
      //midsides II
    case 8:
      fkt = -0.5*gp[1]*(1-gp[0])*(1+gp[2]);
      break;
    case 10:
      fkt = -0.5*gp[1]*(1+gp[0])*(1+gp[2]);
      break;
    case 16:
      fkt = -0.5*gp[1]*(1-gp[0])*(1-gp[2]);
      break;
    case 18:
      fkt = -0.5*gp[1]*(1+gp[0])*(1-gp[2]);
      break;
      //midsides III
    case 12:
      fkt = -0.25*(1-gp[0])*(1-gp[2]*gp[2]);
      break;
    case 13:
      fkt =  0.25*(1-gp[0])*(1-gp[2]*gp[2]);
      break;
    case 14:
      fkt =  0.25*(1+gp[0])*(1-gp[2]*gp[2]);
      break;
    case 15:
      fkt = -0.25*(1+gp[0])*(1-gp[2]*gp[2]);
      break;
    }
    break;


  // ======================================================================
  //    N,zeta(xsi,eta,zeta)
  // ======================================================================
  case N_zeta:
    switch(k)
    {
      //corners I
    case 0:
      fkt = 0.125*(1-gp[0])*(1-gp[1])*( (-gp[0]-gp[1]+gp[2]-2) + (1+gp[2]) );
      break;
    case 1:
      fkt = 0.125*(1-gp[0])*(1+gp[1])*( (-gp[0]+gp[1]+gp[2]-2) + (1+gp[2]) );
      break;
    case 2:
      fkt = 0.125*(1+gp[0])*(1+gp[1])*( (+gp[0]+gp[1]+gp[2]-2) + (1+gp[2]) );
      break;
    case 3:
      fkt = 0.125*(1+gp[0])*(1-gp[1])*( (+gp[0]-gp[1]+gp[2]-2) + (1+gp[2]) );
      break;
      //corners II
    case 4:
      fkt = 0.125*(1-gp[0])*(1-gp[1])*(-(-gp[0]-gp[1]-gp[2]-2) - (1-gp[2]) );
      break;
    case 5:
      fkt = 0.125*(1-gp[0])*(1+gp[1])*(-(-gp[0]+gp[1]-gp[2]-2) - (1-gp[2]) );
      break;
    case 6:
      fkt = 0.125*(1+gp[0])*(1+gp[1])*(-(+gp[0]+gp[1]-gp[2]-2) - (1-gp[2]) );
      break;
    case 7:
      fkt = 0.125*(1+gp[0])*(1-gp[1])*(-(+gp[0]-gp[1]-gp[2]-2) - (1-gp[2]) );
      break;
      //midsides I
    case 9:
      fkt =  0.25*(1-gp[0]*gp[0])*(1+gp[1]);
      break;
    case 11:
      fkt =  0.25*(1-gp[0]*gp[0])*(1-gp[1]);
      break;
    case 17:
      fkt = -0.25*(1-gp[0]*gp[0])*(1+gp[1]);
      break;
    case 19:
      fkt = -0.25*(1-gp[0]*gp[0])*(1-gp[1]);
      break;
      //midsides II
    case 8:
      fkt = +0.25*(1-gp[0])*(1-gp[1]*gp[1]);
      break;
    case 10:
      fkt = +0.25*(1+gp[0])*(1-gp[1]*gp[1]);
      break;
    case 16:
      fkt = -0.25*(1-gp[0])*(1-gp[1]*gp[1]);
      break;
    case 18:
      fkt = -0.25*(1+gp[0])*(1-gp[1]*gp[1]);
      break;
      //midsides III
    case 12:
      fkt = -0.5*gp[2]*(1-gp[0])*(1-gp[1]);
      break;
    case 13:
      fkt = -0.5*gp[2]*(1-gp[0])*(1+gp[1]);
      break;
    case 14:
      fkt = -0.5*gp[2]*(1+gp[0])*(1+gp[1]);
      break;
    case 15:
      fkt = -0.5*gp[2]*(1+gp[0])*(1-gp[1]);
      break;
    }
    break;

    default:
      throw cException(" **** invalid request for derivative **** ", __FILE__, __LINE__);
      break;
  }

  return fkt;
}


double cShapeFunctions::LagrangeHex8(
  const int &k,
  const eDerivative &deriv,
  const cPoint &gp)
{
  double fkt = 0.0;

  switch(deriv)
  {
    case N_fun:
      switch(k)
      {
      case 0:
        fkt = (1.-gp[0]) * (1.-gp[1]) * (1.-gp[2]) * 0.125;
        break;
      case 1:
        fkt = (1.+gp[0]) * (1.-gp[1]) * (1.-gp[2]) * 0.125;
        break;
      case 2:
        fkt = (1.+gp[0]) * (1.+gp[1]) * (1.-gp[2]) * 0.125;
        break;
      case 3:
        fkt = (1.-gp[0]) * (1.+gp[1]) * (1.-gp[2]) * 0.125;
        break;
      case 4:
        fkt = (1.-gp[0]) * (1.-gp[1]) * (1.+gp[2]) * 0.125;
        break;
      case 5:
        fkt = (1.+gp[0]) * (1.-gp[1]) * (1.+gp[2]) * 0.125;
        break;
      case 6:
        fkt = (1.+gp[0]) * (1.+gp[1]) * (1.+gp[2]) * 0.125;
        break;
      case 7:
        fkt = (1.-gp[0]) * (1.+gp[1]) * (1.+gp[2]) * 0.125;
        break;
      }
      break;


    case N_xi:
      switch(k)
      {
      case 0:
        fkt = - (1.-gp[1]) * (1.-gp[2]) * 0.125;
        break;
      case 1:
        fkt =   (1.-gp[1]) * (1.-gp[2]) * 0.125;
        break;
      case 2:
        fkt =   (1.+gp[1]) * (1.-gp[2]) * 0.125;
        break;
      case 3:
        fkt = - (1.+gp[1]) * (1.-gp[2]) * 0.125;
        break;
      case 4:
        fkt = - (1.-gp[1]) * (1.+gp[2]) * 0.125;
        break;
      case 5:
        fkt =   (1.-gp[1]) * (1.+gp[2]) * 0.125;
        break;
      case 6:
        fkt =   (1.+gp[1]) * (1.+gp[2]) * 0.125;
        break;
      case 7:
        fkt = - (1.+gp[1]) * (1.+gp[2]) * 0.125;
        break;
      }
      break;


    case N_eta:
      switch(k)
      {
      case 0:
        fkt = - (1.-gp[0]) * (1.-gp[2]) * 0.125;
        break;
      case 1:
        fkt = - (1.+gp[0]) * (1.-gp[2]) * 0.125;
        break;
      case 2:
        fkt =   (1.+gp[0]) * (1.-gp[2]) * 0.125;
        break;
      case 3:
        fkt =   (1.-gp[0]) * (1.-gp[2]) * 0.125;
        break;
      case 4:
        fkt = - (1.-gp[0]) * (1.+gp[2]) * 0.125;
        break;
      case 5:
        fkt = - (1.+gp[0]) * (1.+gp[2]) * 0.125;
        break;
      case 6:
        fkt =   (1.+gp[0]) * (1.+gp[2]) * 0.125;
        break;
      case 7:
        fkt =   (1.-gp[0]) * (1.+gp[2]) * 0.125;
        break;
      }
      break;


    case N_zeta:
      switch(k)
      {
      case 0:
        fkt = - (1.-gp[0]) * (1.-gp[1]) * 0.125;
        break;
      case 1:
        fkt = - (1.+gp[0]) * (1.-gp[1]) * 0.125;
        break;
      case 2:
        fkt = - (1.+gp[0]) * (1.+gp[1]) * 0.125;
        break;
      case 3:
        fkt = - (1.-gp[0]) * (1.+gp[1]) * 0.125;
        break;
      case 4:
        fkt =   (1.-gp[0]) * (1.-gp[1]) * 0.125;
        break;
      case 5:
        fkt =   (1.+gp[0]) * (1.-gp[1]) * 0.125;
        break;
      case 6:
        fkt =   (1.+gp[0]) * (1.+gp[1]) * 0.125;
        break;
      case 7:
        fkt =   (1.-gp[0]) * (1.+gp[1]) * 0.125;
        break;
      }
      break;

    default:
      throw cException(" **** invalid request for derivative **** ", __FILE__, __LINE__);
      break;
  }

  return fkt;
}


double cShapeFunctions::LagrangeHex27(
  const int &k,
  const eDerivative &deriv,
  const cPoint &gp)
{
  double fkt = 0.0;
  switch(deriv)
  {
    // ======================================================================
    //    N(xi,eta,zeta)
    // ======================================================================
    case N_fun:
      switch(k)
      {
      case 0:
        fkt = -gp[0]*(1.-gp[0]) * gp[1]*(1.-gp[1]) * gp[2]*(1.-gp[2])*0.125;
        break;
      case 8:
        fkt =  (1.-gp[0]*gp[0]) * gp[1]*(1.-gp[1]) * gp[2]*(1.-gp[2])*0.250;
        break;
      case 1:
        fkt =  gp[0]*(1.+gp[0]) * gp[1]*(1.-gp[1]) * gp[2]*(1.-gp[2])*0.125;
        break;
      case 11:
        fkt =  gp[0]*(1.-gp[0]) * (1.-gp[1]*gp[1]) * gp[2]*(1.-gp[2])*0.250;
        break;
      case 21:
        fkt = -(1.-gp[0]*gp[0]) *( 1.-gp[1]*gp[1]) * gp[2]*(1.-gp[2])*0.500;
        break;
      case 9:
        fkt = -gp[0]*(1.+gp[0]) * (1.-gp[1]*gp[1]) * gp[2]*(1.-gp[2])*0.250;
        break;
      case 3:
        fkt =  gp[0]*(1.-gp[0]) * gp[1]*(1.+gp[1]) * gp[2]*(1.-gp[2])*0.125;
        break;
      case 10:
        fkt = -(1.-gp[0]*gp[0]) * gp[1]*(1.+gp[1]) * gp[2]*(1.-gp[2])*0.250;
        break;
      case 2:
        fkt = -gp[0]*(1.+gp[0]) * gp[1]*(1.+gp[1]) * gp[2]*(1.-gp[2])*0.125;
        break;
      // ------------------------------------------------------------------
      case 12:
        fkt =  gp[0]*(1.-gp[0]) * gp[1]*(1.-gp[1]) * (1.-gp[2]*gp[2])*0.250;
        break;
      case 25:
        fkt = -(1.-gp[0]*gp[0]) * gp[1]*(1.-gp[1]) * (1.-gp[2]*gp[2])*0.500;
        break;
      case 13:
        fkt = -gp[0]*(1.+gp[0]) * gp[1]*(1.-gp[1]) * (1.-gp[2]*gp[2])*0.250;
        break;
      case 23:
        fkt = -gp[0]*(1.-gp[0]) * (1.-gp[1]*gp[1]) * (1.-gp[2]*gp[2])*0.500;
        break;
      case 20:
        fkt =  (1.-gp[0]*gp[0]) * (1.-gp[1]*gp[1]) * (1.-gp[2]*gp[2]);
        break;
      case 24:
        fkt =  gp[0]*(1.+gp[0]) * (1.-gp[1]*gp[1]) * (1.-gp[2]*gp[2])*0.500;
        break;
      case 15:
        fkt = -gp[0]*(1.-gp[0]) * gp[1]*(1.+gp[1]) * (1.-gp[2]*gp[2])*0.250;
        break;
      case 26:
        fkt =  (1.-gp[0]*gp[0]) * gp[1]*(1.+gp[1]) * (1.-gp[2]*gp[2])*0.500;
        break;
      case 14:
        fkt =  gp[0]*(1.+gp[0]) * gp[1]*(1.+gp[1]) * (1.-gp[2]*gp[2])*0.250;
        break;
      // ------------------------------------------------------------------
      case 4:
        fkt =  gp[0]*(1.-gp[0]) * gp[1]*(1.-gp[1]) * gp[2]*(1.+gp[2])*0.125;
        break;
      case 16:
        fkt = -(1.-gp[0]*gp[0]) * gp[1]*(1.-gp[1]) * gp[2]*(1.+gp[2])*0.250;
        break;
      case 5:
        fkt = -gp[0]*(1.+gp[0]) * gp[1]*(1.-gp[1]) * gp[2]*(1.+gp[2])*0.125;
        break;
      case 19:
        fkt = -gp[0]*(1.-gp[0]) * (1.-gp[1]*gp[1]) * gp[2]*(1.+gp[2])*0.250;
        break;
      case 22:
        fkt =  (1.-gp[0]*gp[0]) * (1.-gp[1]*gp[1]) * gp[2]*(1.+gp[2])*0.500;
        break;
      case 17:
        fkt =  gp[0]*(1.+gp[0]) * (1.-gp[1]*gp[1]) * gp[2]*(1.+gp[2])*0.250;
        break;
      case 7:
        fkt = -gp[0]*(1.-gp[0]) * gp[1]*(1.+gp[1]) * gp[2]*(1.+gp[2])*0.125;
        break;
      case 18:
        fkt =  (1.-gp[0]*gp[0]) * gp[1]*(1.+gp[1]) * gp[2]*(1.+gp[2])*0.250;
        break;
      case 6:
        fkt =  gp[0]*(1.+gp[0]) * gp[1]*(1.+gp[1]) * gp[2]*(1.+gp[2])*0.125;
        break;
      }
      break; // N_fun


    // ======================================================================
    //    N,xi(xi,eta,zeta)
    // ======================================================================
    case N_xi:
      switch(k)
      {
      case 0:
        fkt = - (1.-2.*gp[0]) * gp[1]*(1.-gp[1]) * gp[2]*(1.-gp[2])*0.125;
        break;
      case 8:
        fkt =    (-2.*gp[0])  * gp[1]*(1.-gp[1]) * gp[2]*(1.-gp[2])*0.250;
        break;
      case 1:
        fkt =   (1.+2.*gp[0]) * gp[1]*(1.-gp[1]) * gp[2]*(1.-gp[2])*0.125;
        break;
      case 11:
        fkt =   (1.-2.*gp[0]) * (1.-gp[1]*gp[1]) * gp[2]*(1.-gp[2])*0.250;
        break;
      case 21:
        fkt = -  (-2.*gp[0])  * (1.-gp[1]*gp[1]) * gp[2]*(1.-gp[2])*0.500;
        break;
      case 9:
        fkt = - (1.+2.*gp[0]) * (1.-gp[1]*gp[1]) * gp[2]*(1.-gp[2])*0.250;
        break;
      case 3:
        fkt =   (1.-2.*gp[0]) * gp[1]*(1.+gp[1]) * gp[2]*(1.-gp[2])*0.125;
        break;
      case 10:
        fkt = -  (-2.*gp[0])  * gp[1]*(1.+gp[1]) * gp[2]*(1.-gp[2])*0.250;
        break;
      case 2:
        fkt = - (1.+2.*gp[0]) * gp[1]*(1.+gp[1]) * gp[2]*(1.-gp[2])*0.125;
        break;
      // ----------------------------------------------------------------------
      case 12:
        fkt =   (1.-2.*gp[0]) * gp[1]*(1.-gp[1]) * (1.-gp[2]*gp[2])*0.250;
        break;
      case 25:
        fkt = -  (-2.*gp[0])  * gp[1]*(1.-gp[1]) * (1.-gp[2]*gp[2])*0.500;
        break;
      case 13:
        fkt = - (1.+2.*gp[0]) * gp[1]*(1.-gp[1]) * (1.-gp[2]*gp[2])*0.250;
        break;
      case 23:
        fkt = - (1.-2.*gp[0]) * (1.-gp[1]*gp[1]) * (1.-gp[2]*gp[2])*0.500;
        break;
      case 20:
        fkt =    (-2.*gp[0])  * (1.-gp[1]*gp[1]) * (1.-gp[2]*gp[2]);
        break;
      case 24:
        fkt =   (1.+2.*gp[0]) * (1.-gp[1]*gp[1]) * (1.-gp[2]*gp[2])*0.500;
        break;
      case 15:
        fkt = - (1.-2.*gp[0]) * gp[1]*(1.+gp[1]) * (1.-gp[2]*gp[2])*0.250;
        break;
      case 26:
        fkt =    (-2.*gp[0])  * gp[1]*(1.+gp[1]) * (1.-gp[2]*gp[2])*0.500;
        break;
      case 14:
        fkt =   (1.+2.*gp[0]) * gp[1]*(1.+gp[1]) * (1.-gp[2]*gp[2])*0.250;
        break;
      // ------------------------------------------------------------------
      case 4:
        fkt =   (1.-2.*gp[0]) * gp[1]*(1.-gp[1]) * gp[2]*(1.+gp[2])*0.125;
        break;
      case 16:
        fkt = -  (-2.*gp[0])  * gp[1]*(1.-gp[1]) * gp[2]*(1.+gp[2])*0.250;
        break;
      case 5:
        fkt = - (1.+2.*gp[0]) * gp[1]*(1.-gp[1]) * gp[2]*(1.+gp[2])*0.125;
        break;
      case 19:
        fkt = - (1.-2.*gp[0]) * (1.-gp[1]*gp[1]) * gp[2]*(1.+gp[2])*0.250;
        break;
      case 22:
        fkt =    (-2.*gp[0])  * (1.-gp[1]*gp[1]) * gp[2]*(1.+gp[2])*0.500;
        break;
      case 17:
        fkt =   (1.+2.*gp[0]) * (1.-gp[1]*gp[1]) * gp[2]*(1.+gp[2])*0.250;
        break;
      case 7:
        fkt = - (1.-2.*gp[0]) * gp[1]*(1.+gp[1]) * gp[2]*(1.+gp[2])*0.125;
        break;
      case 18:
        fkt =    (-2.*gp[0])  * gp[1]*(1.+gp[1]) * gp[2]*(1.+gp[2])*0.250;
        break;
      case 6:
        fkt =   (1.+2.*gp[0]) * gp[1]*(1.+gp[1]) * gp[2]*(1.+gp[2])*0.125;
        break;
      }
      break; // N_xi


    // ======================================================================
    //    N,eta(xi,eta,zeta)
    // ======================================================================
    case N_eta:
      switch(k)
      {
      case 0:
        fkt = -gp[0]*(1.-gp[0]) * (1.-2.*gp[1]) * gp[2]*(1.-gp[2])*0.125;
        break;
      case 8:
        fkt =  (1.-gp[0]*gp[0]) * (1.-2.*gp[1]) * gp[2]*(1.-gp[2])*0.250;
        break;
      case 1:
        fkt =  gp[0]*(1.+gp[0]) * (1.-2.*gp[1]) * gp[2]*(1.-gp[2])*0.125;
        break;
      case 11:
        fkt =  gp[0]*(1.-gp[0]) *  (-2.*gp[1])  * gp[2]*(1.-gp[2])*0.250;
        break;
      case 21:
        fkt = -(1.-gp[0]*gp[0]) *  (-2.*gp[1])  * gp[2]*(1.-gp[2])*0.500;
        break;
      case 9:
        fkt = -gp[0]*(1.+gp[0]) *  (-2.*gp[1])  * gp[2]*(1.-gp[2])*0.250;
        break;
      case 3:
        fkt =  gp[0]*(1.-gp[0]) * (1.+2.*gp[1]) * gp[2]*(1.-gp[2])*0.125;
        break;
      case 10:
        fkt = -(1.-gp[0]*gp[0]) * (1.+2.*gp[1]) * gp[2]*(1.-gp[2])*0.250;
        break;
      case 2:
        fkt = -gp[0]*(1.+gp[0]) * (1.+2.*gp[1]) * gp[2]*(1.-gp[2])*0.125;
        break;
      // ------------------------------------------------------------------
      case 12:
        fkt =  gp[0]*(1.-gp[0]) * (1.-2.*gp[1]) * (1.-gp[2]*gp[2])*0.250;
        break;
      case 25:
        fkt = -(1.-gp[0]*gp[0]) * (1.-2.*gp[1]) * (1.-gp[2]*gp[2])*0.500;
        break;
      case 13:
        fkt = -gp[0]*(1.+gp[0]) * (1.-2.*gp[1]) * (1.-gp[2]*gp[2])*0.250;
        break;
      case 23:
        fkt = -gp[0]*(1.-gp[0]) *  (-2.*gp[1])  * (1.-gp[2]*gp[2])*0.500;
        break;
      case 20:
        fkt =  (1.-gp[0]*gp[0]) *  (-2.*gp[1])  * (1.-gp[2]*gp[2]);
        break;
      case 24:
        fkt =  gp[0]*(1.+gp[0]) *  (-2.*gp[1])  * (1.-gp[2]*gp[2])*0.500;
        break;
      case 15:
        fkt = -gp[0]*(1.-gp[0]) * (1.+2.*gp[1]) * (1.-gp[2]*gp[2])*0.250;
        break;
      case 26:
        fkt =  (1.-gp[0]*gp[0]) * (1.+2.*gp[1]) * (1.-gp[2]*gp[2])*0.500;
        break;
      case 14:
        fkt =  gp[0]*(1.+gp[0]) * (1.+2.*gp[1]) * (1.-gp[2]*gp[2])*0.250;
        break;
      // ------------------------------------------------------------------
      case 4:
        fkt =  gp[0]*(1.-gp[0]) * (1.-2.*gp[1]) * gp[2]*(1.+gp[2])*0.125;
        break;
      case 16:
        fkt = -(1.-gp[0]*gp[0]) * (1.-2.*gp[1]) * gp[2]*(1.+gp[2])*0.250;
        break;
      case 5:
        fkt = -gp[0]*(1.+gp[0]) * (1.-2.*gp[1]) * gp[2]*(1.+gp[2])*0.125;
        break;
      case 19:
        fkt = -gp[0]*(1.-gp[0]) *  (-2.*gp[1])  * gp[2]*(1.+gp[2])*0.250;
        break;
      case 22:
        fkt =  (1.-gp[0]*gp[0]) *  (-2.*gp[1])  * gp[2]*(1.+gp[2])*0.500;
        break;
      case 17:
        fkt =  gp[0]*(1.+gp[0]) *  (-2.*gp[1])  * gp[2]*(1.+gp[2])*0.250;
        break;
      case 7:
        fkt = -gp[0]*(1.-gp[0]) * (1.+2.*gp[1]) * gp[2]*(1.+gp[2])*0.125;
        break;
      case 18:
        fkt =  (1.-gp[0]*gp[0]) * (1.+2.*gp[1]) * gp[2]*(1.+gp[2])*0.250;
        break;
      case 6:
        fkt =  gp[0]*(1.+gp[0]) * (1.+2.*gp[1]) * gp[2]*(1.+gp[2])*0.125;
        break;
      }
      break; // N_eta


    // ======================================================================
    //    N,zeta(xi,eta,zeta)
    // ======================================================================
    case N_zeta:
      switch(k)
      {
      case 0:
        fkt = -gp[0]*(1.-gp[0]) * gp[1]*(1.-gp[1]) * (1.-2.*gp[2])*0.125;
        break;
      case 8:
        fkt =  (1.-gp[0]*gp[0]) * gp[1]*(1.-gp[1]) * (1.-2.*gp[2])*0.250;
        break;
      case 1:
        fkt =  gp[0]*(1.+gp[0]) * gp[1]*(1.-gp[1]) * (1.-2.*gp[2])*0.125;
        break;
      case 11:
        fkt =  gp[0]*(1.-gp[0]) * (1.-gp[1]*gp[1]) * (1.-2.*gp[2])*0.250;
        break;
      case 21:
        fkt = -(1.-gp[0]*gp[0]) * (1.-gp[1]*gp[1]) * (1.-2.*gp[2])*0.500;
        break;
      case 9:
        fkt = -gp[0]*(1.+gp[0]) * (1.-gp[1]*gp[1]) * (1.-2.*gp[2])*0.250;
        break;
      case 3:
        fkt =  gp[0]*(1.-gp[0]) * gp[1]*(1.+gp[1]) * (1.-2.*gp[2])*0.125;
        break;
      case 10:
        fkt = -(1.-gp[0]*gp[0]) * gp[1]*(1.+gp[1]) * (1.-2.*gp[2])*0.250;
        break;
      case 2:
        fkt = -gp[0]*(1.+gp[0]) * gp[1]*(1.+gp[1]) * (1.-2.*gp[2])*0.125;
        break;
      // ------------------------------------------------------------------
      case 12:
        fkt =  gp[0]*(1.-gp[0]) * gp[1]*(1.-gp[1]) *  (-2.*gp[2]) *0.250;
        break;
      case 25:
        fkt = -(1.-gp[0]*gp[0]) * gp[1]*(1.-gp[1]) *  (-2.*gp[2]) *0.500;
        break;
      case 13:
        fkt = -gp[0]*(1.+gp[0]) * gp[1]*(1.-gp[1]) *  (-2.*gp[2]) *0.250;
        break;
      case 23:
        fkt = -gp[0]*(1.-gp[0]) * (1.-gp[1]*gp[1]) *  (-2.*gp[2]) *0.500;
        break;
      case 20:
        fkt =  (1.-gp[0]*gp[0]) * (1.-gp[1]*gp[1]) *  (-2.*gp[2]);
        break;
      case 24:
        fkt =  gp[0]*(1.+gp[0]) * (1.-gp[1]*gp[1]) *  (-2.*gp[2]) *0.500;
        break;
      case 15:
        fkt = -gp[0]*(1.-gp[0]) * gp[1]*(1.+gp[1]) *  (-2.*gp[2]) *0.250;
        break;
      case 26:
        fkt =  (1.-gp[0]*gp[0]) * gp[1]*(1.+gp[1]) *  (-2.*gp[2]) *0.500;
        break;
      case 14:
        fkt =  gp[0]*(1.+gp[0]) * gp[1]*(1.+gp[1]) *  (-2.*gp[2]) *0.250;
        break;
      // ------------------------------------------------------------------
      case 4:
        fkt =  gp[0]*(1.-gp[0]) * gp[1]*(1.-gp[1]) * (1.+2.*gp[2])*0.125;
        break;
      case 16:
        fkt = -(1.-gp[0]*gp[0]) * gp[1]*(1.-gp[1]) * (1.+2.*gp[2])*0.250;
        break;
      case 5:
        fkt = -gp[0]*(1.+gp[0]) * gp[1]*(1.-gp[1]) * (1.+2.*gp[2])*0.125;
        break;
      case 19:
        fkt = -gp[0]*(1.-gp[0]) * (1.-gp[1]*gp[1]) * (1.+2.*gp[2])*0.250;
        break;
      case 22:
        fkt =  (1.-gp[0]*gp[0]) * (1.-gp[1]*gp[1]) * (1.+2.*gp[2])*0.500;
        break;
      case 17:
        fkt =  gp[0]*(1.+gp[0]) * (1.-gp[1]*gp[1]) * (1.+2.*gp[2])*0.250;
        break;
      case 7:
        fkt = -gp[0]*(1.-gp[0]) * gp[1]*(1.+gp[1]) * (1.+2.*gp[2])*0.125;
        break;
      case 18:
        fkt =  (1.-gp[0]*gp[0]) * gp[1]*(1.+gp[1]) * (1.+2.*gp[2])*0.250;
        break;
      case 6:
        fkt =  gp[0]*(1.+gp[0]) * gp[1]*(1.+gp[1]) * (1.+2.*gp[2])*0.125;
        break;
      }
      break; // N_zeta

    default:
      throw cException(" **** invalid request for derivative **** ", __FILE__, __LINE__);
      break;
  }

  return fkt;
}


// ---------------------------------------------------------------------------
//   test functions for tria elements (3,6 and 9 nodes)
// ---------------------------------------------------------------------------
PetscReal cShapeFunctions::TriaN(
  const int &nnod,
  const int &k,
  const eDerivative &deriv,
  const cPoint &gp)
{
  switch (nnod)
  {
    case 3:
      return Tria3(k, deriv, gp);
      break;
    case 6:
      return Tria6(k, deriv, gp);
      break;
    case 9:
      return Tria9(k, deriv, gp);
      break;

  }

  return 0.0;
}


// ----------------------------------------------------------------------------
//
//           ^ eta
//           |
//           |
//         2 o
//           |.
//           | .
//           |  .
//           |   .
//           |    .
//           o-----o  --->*xsi
//          0       1
//
// ----------------------------------------------------------------------------
PetscReal cShapeFunctions::Tria3(
  const int &k,
  const eDerivative &deriv,
  const cPoint &gp)
{
  PetscReal fkt = 0.0;

  switch(deriv)
  {
    case N_fun:
    {
      switch(k)
      {
        case 0:
          fkt = 1.0 - gp[0] - gp[1];
          break;
        case 1:
          fkt = gp[0];
          break;
        case 2:
          fkt = gp[1];
          break;
      }
      break;
    }

    case N_xi:
    {
      switch(k)
      {
        case 0:
          fkt = -1.0;
          break;
        case 1:
          fkt =  1.0;
          break;
        case 2:
          fkt =  0.0;
          break;
      }
      break;
    }

    case N_eta:
    {
      switch(k)
      {
        case 0:
          fkt = -1.0;
          break;
        case 1:
          fkt =  0.0;
          break;
        case 2:
          fkt =  1.0;
          break;
      }
      break;
    }

    default:
      throw cException(" **** invalid request for derivative **** ", __FILE__, __LINE__);
  }

  return fkt;
}



// ----------------------------------------------------------------------------
//
//           ^ eta
//           |
//           |
//         2 o
//           |.
//           | .
//         5 o  o 4
//           |   .
//           |    .
//           o--o--o  --->*xsi
//          0   3   1
//
// ----------------------------------------------------------------------------
PetscReal cShapeFunctions::Tria6(
  const int &k,
  const eDerivative &deriv,
  const cPoint &gp)
{
  PetscReal fkt = 0.0;

  switch(deriv)
  {
    case N_fun:
    {
      switch(k)
      {
        case 0:
          fkt = (2.*gp[0] - 1.)*gp[0];
          break;
        case 1:
          fkt = (2.*gp[1] - 1.)*gp[1];
          break;
        case 2:
          fkt = (2.*gp[2] - 1.)*gp[2];
          break;
        case 3:
          fkt = 4.*gp[0]*gp[1];
          break;
        case 4:
          fkt = 4.0*gp[1]*gp[2];
          break;
        case 5:
          fkt = 4.0*gp[0]*gp[2];
          break;

      }
      break;
    }

    case N_xi:
    {
      switch(k)
      {
        case 0:
          fkt = 1. - 4.*gp[0];
          break;
        case 1:
          fkt = 4.*gp[1] - 1.;
          break;
        case 2:
          fkt = 0.;
          break;
        case 3:
          fkt = 4.*(gp[0] - gp[1]);
          break;
        case 4:
          fkt = 4.*gp[2];
          break;
        case 5:
          fkt = -4.*gp[2];
          break;
      }
      break;
    }
       
   case N_eta:
    {
      switch(k)
      {
        case 0:
          fkt = 1. - 4.*gp[0];
          break;
        case 1:
          fkt = 0.;
          break;
        case 2:
          fkt = 4.*gp[2] - 1.;
          break;
        case 3:
          fkt = -4.0*gp[1];
          break;
        case 4:
          fkt = 4.*gp[1];
          break;
        case 5:
          fkt = 4.*(gp[0] - gp[2]);
          break;
      }
      break;
    }

    default:
      throw cException(" **** invalid request for derivative **** ", __FILE__, __LINE__);
  }

  return fkt;
}


// ----------------------------------------------------------------------------
//
//           ^ eta
//         2 o
//           | .
//           |  .
//         7 o   o 6
//           |    .
//           |     o 5
//         8 o      .
//           |       .
//           o--o--o--o  --->*xsi
//          0   3  4  1     
//
// ----------------------------------------------------------------------------
PetscReal cShapeFunctions::Tria9(
  const int &k,
  const eDerivative &deriv,
  const cPoint &gp)
{
  PetscReal fkt = 0.0;

  switch(deriv)
  {
    case N_fun:
    {
      switch(k)
      {
        case 0:
          fkt = 0.5*(3.*gp[0] - 1.)*(3.*gp[0] - 2.)*gp[0];
          break;
        case 1:
          fkt = 0.5*(3.*gp[1] - 1.)*(3.*gp[1] - 2.)*gp[1];
          break;
        case 2:
          fkt = 0.5*(3.*gp[2] - 1.)*(3.*gp[2] - 2.)*gp[2];
          break;
        case 3:
          fkt = 4.5*gp[0]*gp[1]*(3.*gp[0] - 1.);
          break;
        case 4:
          fkt = 4.5*gp[0]*gp[1]*(3.*gp[1] - 1.);
          break;
        case 5:
          fkt = 4.5*gp[1]*gp[2]*(3.*gp[1] - 1.);
          break;
        case 6:
          fkt = 4.5*gp[1]*gp[2]*(3.*gp[2] - 1.);
          break;
        case 7:
          fkt = 4.5*gp[0]*gp[2]*(3.*gp[2] - 1.);
          break;
        case 8:
          fkt = 4.5*gp[0]*gp[2]*(3.*gp[0] - 1.);
          break;

      }
      break;
    }

    case N_xi:
    {
      switch(k)
      {
        case 0:
          fkt = -0.5*(27.*(gp[0]*gp[0]) - 18.*gp[0] + 2.);
          break;
        case 1:
          fkt = 0.5*(27.*(gp[1]*gp[1]) - 18.*gp[1] + 2.);; 
          break;
        case 2:
          fkt = 0.;
          break;
        case 3:
          fkt = 4.5*(3.*(gp[0]*gp[0]) - gp[0] - 6.*(gp[0]*gp[1]) + gp[1]);
          break;
        case 4:
          fkt = 4.5*(-3.*(gp[1]*gp[1]) + gp[1] + 6.*(gp[0]*gp[1]) - gp[0]);
          break;
        case 5:
          fkt = 4.5*(6.*(gp[1]*gp[2]) - gp[2]);
          break;
        case 6:
          fkt = 4.5*(3.*(gp[2]*gp[2]) - gp[2]);
          break;
        case 7:
          fkt = -4.5*(3.*(gp[2]*gp[2]) - gp[2]);
          break;
        case 8:
          fkt = -4.5*(6.*(gp[0]*gp[2]) - gp[2]);
          break;

      }
      break;
    }

    case N_eta:
    {
      switch(k)
      {
        case 0:
          fkt = -0.5*(27.*(gp[0]*gp[0]) - 18.*gp[0] + 2.);
          break;
        case 1:
          fkt =  0.;
          break;
        case 2:
          fkt = 0.5*(27.*(gp[2]*gp[2]) - 18.*gp[2] + 2.); 
          break;
        case 3:
          fkt = -4.5*(6.*(gp[0]*gp[1]) - gp[1]);
          break;
        case 4:
          fkt = -4.5*(3.*(gp[1]*gp[1]) - gp[1]);
          break;
        case 5:
          fkt = 4.5*(3.*(gp[1]*gp[1]) - gp[1]);
          break;
        case 6:
          fkt = 4.5*(6.*(gp[1]*gp[2]) - gp[1]);
          break;
        case 7:
          fkt = 4.5*(-3.*(gp[2]*gp[2]) + gp[2] + 6.*(gp[0]*gp[2]) - gp[0]);
          break;
        case 8:
          fkt = 4.5*(3.*(gp[0]*gp[0]) - gp[0] - 6.*(gp[0]*gp[2]) + gp[2]);
          break;
      }
      break;
    }

    default:
      throw cException(" **** invalid request for derivative **** ", __FILE__, __LINE__);
  }

  return fkt;
}



// ---------------------------------------------------------------------------
//   test functions for tetraehedron elements (4 or 10 or 16 nodes)
//   Thomas Huber
//   the tetrahedrons's local coordinates are L1, L2, L3 and L4  
// ---------------------------------------------------------------------------
PetscReal cShapeFunctions::TetN(
  const int &nnod,
  const int &k,
  const eDerivative &deriv,
  const cPoint &gp)
{
  switch (nnod)
  {
    case 4:
      return Tet4(k, deriv, gp);
      break;
    case 10:
      return Tet10(k, deriv, gp);
      break;
    case 16:
      return Tet16(k, deriv, gp);
      break;
    case 41:
      return Tet4L(k, deriv, gp);
      break;
    case 101:
      return Tet10L(k, deriv, gp);
      break;
    case 161:
      return Tet16L(k, deriv, gp);
      break;
  }

  return 0.0;
 }


PetscReal cShapeFunctions::Tet4(
   const int &k,
   const eDerivative &deriv,
   const cPoint &gp)
{
  PetscReal fkt = 0.0;

  switch(deriv) 
  {
    case N_fun:
      switch (k)
      {
        case 0:
          fkt = 1.0 - gp[0] - gp[1]- gp[2];
          break;
        case 1:
          fkt = gp[0];
          break;
        case 2:
          fkt = gp[1]; 
          break;
        case 3:
          fkt = gp[2];
          break; 
      }
      break;

    case N_xi:
      switch (k)
      {
      case 0:
        fkt = -1.0;
        break;
      case 1:
        fkt = 1.0;
        break;
      case 2:
        fkt = 0.0;
        break;
      case 3:
        fkt = 0.0;
        break;
      }
      break;
    case N_eta:
      switch (k)
      {
      case 0:
        fkt = -1.0;
        break;
      case 1:
        fkt = 0.0;
        break;
      case 2:
        fkt = 1.0;
        break;
      case 3:
        fkt = 0.0;
        break;
      }
      break;

    case N_zeta:
      switch (k)
      {
      case 0:
        fkt = -1.0;
        break;
      case 1:
        fkt = 0.0;
        break;
      case 2:
        fkt = 0.0;
        break;
      case 3:
        fkt = 1.0;
        break;
      }
      break;

    default:
     throw cException(" **** invalid request for derivative **** ", __FILE__, __LINE__);
      break;
  }

  return fkt;
}

PetscReal cShapeFunctions::Tet4L(
   const int &k,
   const eDerivative &deriv,
   const cPoint &gp)
{
  PetscReal fkt = 0.0;

  switch(deriv) 
  {
    case N_fun:
      switch (k)
      {
        case 0:
          fkt = gp[0];
          break;
        case 1:
          fkt = gp[1];
          break;
        case 2:
          fkt = gp[2]; 
          break;
        case 3:
          fkt = gp[3];
          break; 
      }
      break;

    case N_xi:
      switch (k)
      {
      case 0:
        fkt = -1.0;
        break;
      case 1:
        fkt = 1.0;
        break;
      case 2:
        fkt = 0.0;
        break;
      case 3:
        fkt = 0.0;
        break;
      }
      break;

    case N_eta:
      switch (k)
      {
      case 0:
        fkt = -1.0;
        break;
      case 1:
        fkt = 0.0;
        break;
      case 2:
        fkt = 1.0;
        break;
      case 3:
        fkt = 0.0;
        break;
      }
      break;

    case N_zeta:
      switch (k)
      {
      case 0:
        fkt = -1.0;
        break;
      case 1:
        fkt = 0.0;
        break;
      case 2:
        fkt = 0.0;
        break;
      case 3:
        fkt = 1.0;
        break;
      }
      break;

    default:
      throw cException(" **** invalid request for derivative **** ", __FILE__, __LINE__);
      break;
  }

  return fkt;
}

PetscReal cShapeFunctions::Tet10(
   const int &k,
   const eDerivative &deriv,
   const cPoint &gp)
{
  PetscReal fkt = 0.0;

  switch(deriv) 
  {
    case N_fun:
      switch (k)
      {
        case 0:
          fkt = (1.- gp[0] - gp[1] - gp[2])*(2.*(1. - gp[0] - gp[1] - gp[2]) - 1.);
          break;
        case 1:
          fkt = gp[0]*(2.*gp[0] - 1.);
          break;
        case 2:
          fkt = gp[1]*(2.*gp[1] - 1.); 
          break;
        case 3:
          fkt = gp[2]*(2.*gp[2] - 1.);
          break;
        case 4:
          fkt = 4.*(1. - gp[0] - gp[1] - gp[2])*gp[0];
          break;
        case 5:
          fkt = 4.*gp[0]*gp[1];
          break;
        case 6:
          fkt = 4.*(1.- gp[0] - gp[1] - gp[2])*gp[1]; 
          break;
        case 7:
          fkt = 4.*(1. - gp[0] - gp[1] - gp[2])*gp[2];
          break;
        case 8:
          fkt = 4.*gp[0]*gp[2];
          break;
        case 9:
          fkt = 4.*gp[1]*gp[2];
          break;
      }
      break;

    case N_xi:
      switch (k)
      {
        case 0:
          fkt = 4.*(gp[0] + gp[1] + gp[2]) - 3.;
          break;
        case 1:
          fkt = 4.*gp[0] - 1.;
          break;
        case 2:
          fkt = 0.;
          break;
        case 3:
          fkt = 0.;
          break;
        case 4:
          fkt = 4. - 4.*gp[1] - 4.*gp[2] - 8.*gp[0];
          break;
        case 5:
          fkt = 4.*gp[1];
          break;
        case 6:
          fkt = -4.*gp[1];
          break;
        case 7:
          fkt = -4.*gp[2];
          break;
        case 8:
          fkt = 4.*gp[2];
          break;
        case 9:
          fkt = 0.;
          break;
      }
      break;

   case N_eta:
      switch (k)
      {
        case 0:
          fkt = 4.*(gp[0] + gp[1] + gp[2]) - 3.;
          break;
        case 1:
          fkt = 0.;
          break;
        case 2:
          fkt = 4.*gp[1] - 1.;
          break;
        case 3:
          fkt = 0.;
          break;
        case 4:
          fkt = -4.*gp[0];
          break;
        case 5:
          fkt = 4.*gp[0];
          break;
        case 6:
          fkt = 4. - 8.*gp[1] - 4.*gp[2] - 4.*gp[0];
          break;
        case 7:
          fkt = -4.*gp[2];
          break;
        case 8:
          fkt = 0.;
          break;
        case 9:
          fkt = 4.*gp[2];
          break;
      }
      break;

   case N_zeta:
      switch (k)
      {
        case 0:
          fkt = 4.*(gp[0] + gp[1] + gp[2]) - 3.;
          break;
        case 1:
          fkt = 0.;
          break;
        case 2:
          fkt = 0.;
          break;
        case 3:
          fkt = 4.*gp[2] - 1.;
          break;
        case 4:
          fkt = -4.*gp[0];
          break;
        case 5:
          fkt = 0.;
          break;
        case 6:
          fkt = -4.*gp[1];
          break;
        case 7:
          fkt = 4. - 4.*gp[1] - 8.*gp[2] - 4.*gp[0];
          break;
        case 8:
          fkt = 4.*gp[0];
          break;
        case 9:
          fkt = 4.*gp[1];
          break;
      }
      break;

    default:
      throw cException(" **** invalid request for derivative **** ", __FILE__, __LINE__);
      break;
  }
  return fkt;
}  

PetscReal cShapeFunctions::Tet10L(
   const int &k,
   const eDerivative &deriv,
   const cPoint &gp)
{
  PetscReal fkt = 0.0;

  switch(deriv) 
  {
    case N_fun:
      switch (k)
      {
        case 0:
          fkt = (2.*gp[0] - 1.)*gp[0];
          break;
        case 1:
          fkt = (2.*gp[1] - 1.)*gp[1];
          break;
        case 2:
          fkt = (2.*gp[2] - 1.)*gp[2];
          break;
        case 3:
          fkt = (2.*gp[3] - 1.)*gp[3];
          break;
        case 4:
          fkt = 4.*gp[0]*gp[1];
          break;
        case 5:
          fkt = 4.*gp[1]*gp[2];
          break;
        case 6:
          fkt = 4.*gp[0]*gp[2]; 
          break;
        case 7:
          fkt = 4.*gp[0]*gp[3];
          break;
        case 8:
          fkt = 4.*gp[1]*gp[3];
          break;
        case 9:
          fkt = 4.*gp[2]*gp[3];
          break;
      }
      break;

    case N_xi:
      switch (k)
      {
        case 0:
          fkt = 1. - 4.*gp[0];
          break;
        case 1:
          fkt = 4.*gp[1] - 1.;
          break;
        case 2:
          fkt = 0.;
          break;
        case 3:
          fkt = 0.;
          break;
        case 4:
          fkt = 4.*(gp[0] - gp[1]);
          break;
        case 5:
          fkt = 4.*gp[2];
          break;
        case 6:
          fkt = -4.*gp[2];
          break;
        case 7:
          fkt = -4.*gp[3];
          break;
        case 8:
          fkt = 4.*gp[3];
          break;
        case 9:
          fkt = 0.;
          break;
      }
      break;

   case N_eta:
      switch (k)
      {
        case 0:
          fkt = 1. - 4.*gp[0];
          break;
        case 1:
          fkt = 0.;
          break;
        case 2:
          fkt = 4.*gp[2] - 1.;
          break;
        case 3:
          fkt = 0.;
          break;
        case 4:
          fkt = -4.*gp[1];
          break;
        case 5:
          fkt = 4.*gp[1];
          break;
        case 6:
          fkt = 4.*(gp[0] - gp[2]) ;
          break;
        case 7:
          fkt = -4.*gp[3];
          break;
        case 8:
          fkt = 0.;
          break;
        case 9:
          fkt = 4.*gp[3];
          break;
      }
      break;

   case N_zeta:
      switch (k)
      {
        case 0:
          fkt = 1. - 4.*gp[0];
          break;
        case 1:
          fkt = 0.;
          break;
        case 2:
          fkt = 0.;
          break;
        case 3:
          fkt = 4.*gp[3] - 1.;
          break;
        case 4:
          fkt = -4.*gp[1];
          break;
        case 5:
          fkt = 0.;
          break;
        case 6:
          fkt = -4.*gp[2];
          break;
        case 7:
          fkt = 4.*(gp[0] - gp[3]);
          break;
        case 8:
          fkt = 4.*gp[1];
          break;
        case 9:
          fkt = 4.*gp[2];
          break;
      }
      break;

    default:
      throw cException(" **** invalid request for derivative **** ", __FILE__, __LINE__);
      break;
  }
  return fkt;
}  

PetscReal cShapeFunctions::Tet16(
   const int &k,
   const eDerivative &deriv,
   const cPoint &gp)
{
  PetscReal fkt = 0.0;

  switch(deriv) 
  {
    case N_fun:
      switch (k)
      {
        case 0:
          fkt = 0.5*(1. - gp[0] - gp[1] - gp[2])*((3.*(1. - gp[0] - gp[1] - gp[2])) - 1.)*((3.*(1.- gp[0] - gp[1] - gp[2])) - 2.);
          break;
        case 1:
          fkt = 0.5*gp[0]*(3.*gp[0] - 1.)*(3.*gp[0] - 2.);
          break;
        case 2:
          fkt = 0.5*gp[1]*(3.*gp[1] - 1.)*(3.*gp[1] - 2.);
          break;
        case 3:
          fkt = 0.5*gp[2]*(3.*gp[2] - 1.)*(3.*gp[2] - 2.);
          break;
        case 4:
          fkt = 4.5*(1. - gp[0] - gp[1] - gp[2])*gp[0]*((3.*(1. - gp[0] - gp[1] - gp[2])) - 1.);
          break;
        case 5:
          fkt = 4.5*(1. - gp[0] - gp[1] - gp[2])*gp[0]*(3.*gp[0] - 1.);
          break;
        case 6:
          fkt = 4.5*gp[1]*gp[0]*(3.*gp[0] - 1.);
          break;
        case 7:
          fkt = 4.5*gp[1]*gp[0]*(3.*gp[1] - 1.);
          break;
        case 8:
          fkt = 4.5*(1. - gp[0] - gp[1] - gp[2])*gp[1]*(3.*gp[1] - 1.);
          break;
        case 9:
          fkt = 4.5*(1. - gp[0] - gp[1] - gp[2])*gp[1]*((3.*(1. - gp[0] - gp[1] - gp[2])) - 1.);
          break;
        case 10:
          fkt = 4.5*(1. - gp[0] - gp[1] - gp[2])*gp[2]*((3.*(1. - gp[0] - gp[1] - gp[2])) - 1.);
          break;
        case 11:
          fkt = 4.5*gp[2]*gp[0]*(3.*gp[0] - 1.);
          break;
        case 12:
          fkt = 4.5*gp[1]*gp[2]*(3.*gp[1] - 1.);
          break;
        case 13:
          fkt = 4.5*(1. - gp[0] - gp[1] - gp[2])*gp[2]*(3.*gp[2] - 1.);
          break;
        case 14:
          fkt = 4.5*gp[2]*gp[0]*(3.*gp[2] - 1.);
          break;
        case 15:
          fkt = 4.5*gp[1]*gp[2]*(3.*gp[2] - 1.);
          break; 
      }
      break;

    case N_xi:
      switch (k)
      {
        case 0:
          fkt = -0.5*(27*((1 - gp[0] - gp[1] - gp[2])*(1 - gp[0] - gp[1] - gp[2]))-18*(1 - gp[0] - gp[1] - gp[2]) + 2.);
          break;
        case 1:
          fkt = 0.5*(27*gp[0]*gp[0] - 18*gp[0] + 2);
          break;
        case 2:
          fkt = 0.0;
          break;
        case 3:
          fkt = 0.0;
          break;
        case 4:
          fkt = 4.5*(3*(1 - gp[0] - gp[1] - gp[2])*(1 - gp[0] - gp[1] - gp[2]) - (1 - gp[0] - gp[1] - gp[2]) - 6*(1 - gp[0] - gp[1] - gp[2])*gp[0] + gp[0]);
          break;
        case 5:
          fkt = 4.5*(-3*gp[0]*gp[0] + gp[0] + 6*(1 - gp[0] - gp[1] - gp[2])*gp[0] - (1 - gp[0] - gp[1] - gp[2]));
          break;
        case 6:
          fkt = 4.5*(6*gp[0]*gp[1] - gp[1]);
          break;
        case 7:
          fkt = 4.5*(3*gp[1]*gp[1] - gp[1]);
          break;
        case 8:
          fkt = -4.5*(3*gp[1]*gp[1] - gp[1]);
          break;
        case 9:
          fkt = -4.5*(6*(1 - gp[0] - gp[1] - gp[2])*gp[1] - gp[1]);
          break;
        case 10:
          fkt = -4.5*(6*(1 - gp[0] - gp[1] - gp[2])*gp[2] - gp[2]);
          break;
        case 11:
          fkt = 4.5*(6*gp[0]*gp[2] - gp[2]);
          break;
        case 12:
          fkt = 0.0;
          break;
        case 13:
          fkt = -4.5*(3*gp[2]*gp[2] - gp[2]);
          break;
        case 14:
          fkt = 4.5*(3*gp[2]*gp[2] - gp[2]);
          break;
        case 15:
          fkt = 0.0;
          break;
      }
      break;

   case N_eta:
      switch (k)
      {
        case 0:
          fkt = -0.5*(27*((1 - gp[0] - gp[1] - gp[2])*(1 - gp[0] - gp[1] - gp[2]))-18*(1 - gp[0] - gp[1] - gp[2]) + 2.);
          break;
        case 1:
          fkt = 0.0;
          break;
        case 2:
          fkt = 0.5*(27*gp[1]*gp[1] - 18*gp[1] + 2);
          break;
        case 3:
          fkt = 0.0;
          break;
        case 4:
          fkt = -4.5*(6*(1 - gp[0] - gp[1] - gp[2])*gp[0] - gp[0]);
          break;
        case 5:
          fkt =-4.5*(3*gp[0]*gp[0] - gp[0]);
          break;
        case 6:
          fkt = 4.5*(3*gp[0]*gp[0] - gp[0]);
          break;
        case 7:
          fkt = 4.5*(6*gp[0]*gp[1] - gp[0]);
          break;
        case 8:
          fkt = 4.5*(-3*gp[1]*gp[1] + gp[1] + 6*(1 - gp[0] - gp[1] - gp[2])*gp[1] - (1 - gp[0] - gp[1] - gp[2]));
          break;
        case 9:
          fkt = 4.5*(3*(1 - gp[0] - gp[1] - gp[2])*(1 - gp[0] - gp[1] - gp[2]) - (1 - gp[0] - gp[1] - gp[2]) - 6*(1 - gp[0] - gp[1] - gp[2])*gp[1] + gp[1]);
          break;
        case 10:
          fkt = -4.5*(6*(1 - gp[0] - gp[1] - gp[2])*gp[2] - gp[2]);
          break;
        case 11:
          fkt = 0.0;
          break;
        case 12:
          fkt = 4.5*(6*gp[1]*gp[2] - gp[2]);
          break;
        case 13:
          fkt = -4.5*(3*gp[2]*gp[2] - gp[2]);
          break;
        case 14:
          fkt = 0.0;
          break;
        case 15:
          fkt = 4.5*(3*gp[2]*gp[2] - gp[2]);
          break;
      }
      break;

   case N_zeta:
      switch (k)
      {
        case 0:
          fkt = -0.5*(27*((1 - gp[0] - gp[1] - gp[2])*(1 - gp[0] - gp[1] - gp[2]))-18*(1 - gp[0] - gp[1] - gp[2]) + 2.);
          break;
        case 1:
          fkt = 0.0;
          break;
        case 2:
          fkt = 0.0;
          break;
        case 3:
          fkt = 0.5*(27*gp[2]*gp[2] - 18*gp[2] + 2);
          break;
        case 4:
          fkt = -4.5*(6*(1 - gp[0] - gp[1] - gp[2])*gp[0] - gp[0]);
          break;
        case 5:
          fkt = -4.5*(3*gp[0]*gp[0] - gp[0]);
          break;
        case 6:
          fkt = 0.0;
          break;
        case 7:
          fkt = 0.0;
          break;
        case 8:
          fkt = -4.5*(3*gp[1]*gp[1] - gp[1]);
          break;
        case 9:
          fkt = -4.5*(6*(1 - gp[0] - gp[1] - gp[2])*gp[1] - gp[1]);
          break;
        case 10:
          fkt = 4.5*(3*(1 - gp[0] - gp[1] - gp[2])*(1 - gp[0] - gp[1] - gp[2]) - (1 - gp[0] - gp[1] - gp[2]) - 6*(1 - gp[0] - gp[1] - gp[2])*gp[2] + gp[2]);
          break;
        case 11:
          fkt = 4.5*(3*gp[0]*gp[0] - gp[0]);
          break;
        case 12:
          fkt = 4.5*(3*gp[1]*gp[1] - gp[1]);
          break;
        case 13:
          fkt = 4.5*(-3*gp[2]*gp[2] + gp[2] + 6*(1 - gp[0] - gp[1] - gp[2])*gp[2] - (1 - gp[0] - gp[1] - gp[2]));
          break;
        case 14:
          fkt = 4.5*(6*gp[0]*gp[2] - gp[0]);
          break;
        case 15:
          fkt = 4.5*(6*gp[1]*gp[2] - gp[1]);
          break;
      }
      break;

    default:
      throw cException(" **** invalid request for derivative **** ", __FILE__, __LINE__);
      break;
  }
  return fkt;
}  

PetscReal cShapeFunctions::Tet16L(
   const int &k,
   const eDerivative &deriv,
   const cPoint &gp)
{
  PetscReal fkt = 0.0;
  switch(deriv) 
   {
    case N_fun:
      switch (k)
      {
        case 0:
          fkt = 0.5*(3.*gp[0] - 1.)*(3.*gp[0] - 2.)*gp[0];
          break;
        case 1:
          fkt = 0.5*(3.*gp[1] - 1.)*(3.*gp[1] - 2.)*gp[1];
          break;
        case 2:
          fkt = 0.5*(3.*gp[2] - 1.)*(3.*gp[2] - 2.)*gp[2];
          break;
        case 3:
          fkt = 0.5*(3.*gp[3] - 1.)*(3.*gp[3] - 2.)*gp[3];
          break;
        case 4:
          fkt = 4.5*(gp[0]*gp[1])*(3.*gp[0] - 1.);
          break;
        case 5:
          fkt = 4.5*(gp[0]*gp[1])*(3.*gp[1] - 1.);
          break;
        case 6:
          fkt = 4.5*(gp[1]*gp[2])*(3.*gp[1] - 1.);
          break;
        case 7:
          fkt = 4.5*(gp[1]*gp[2])*(3.*gp[2] - 1.); 
          break;
        case 8:
          fkt = 4.5*(gp[0]*gp[2])*(3.*gp[2] - 1.);
          break;
        case 9:
          fkt = 4.5*(gp[0]*gp[2])*(3.*gp[0] - 1.);
          break;
        case 10:
          fkt = 4.5*(gp[0]*gp[3])*(3.*gp[0] - 1.);
          break;
        case 11:
          fkt = 4.5*(gp[1]*gp[3])*(3.*gp[1] - 1.);
          break;
        case 12:
          fkt = 4.5*(gp[2]*gp[3])*(3.*gp[2] - 1.);
          break;
        case 13:
          fkt = 4.5*(gp[0]*gp[3])*(3.*gp[3] - 1.);
          break;
        case 14:
          fkt = 4.5*(gp[1]*gp[3])*(3.*gp[3] - 1.);
          break;
        case 15:
          fkt = 4.5*(gp[2]*gp[3])*(3.*gp[3] - 1.);
          break;
      }
      break;
   
   /*case N_xi:
      switch (k)
      {
        case 0:
          fkt = -((3*gp[1] + 3*gp[0] + 3*gp[2] - 1)*(3*gp[1] + 3*gp[0] + 3*gp[2] - 2))/2 - 3*(3*gp[1] + 3*gp[0] + 3*gp[2] - 1)*(gp[1]/2 + gp[0]/2 + gp[2]/2 - 1/2) - 3*(3*gp[1] + 3*gp[0] + 3*gp[2] - 2)*(gp[1]/2 + gp[0]/2 + gp[2]/2 - 1/2);
          break;
        case 1:
          fkt = (3*gp[0]*(3*gp[0] - 1))/2 + (3*gp[0]*(3*gp[0] - 2))/2 + ((3*gp[0] - 1)*(3*gp[0] - 2))/2;
          break;
        case 2:
          fkt = 0.;
          break;
        case 3:
          fkt = 0.;
          break;
        case 4:
          fkt = (3*gp[1] + 3*gp[0] + 3*gp[2] - 2)*((9*gp[1])/2 + (9*gp[0])/2 + (9*gp[2])/2 - 9/2) + (9*gp[0]*(3*gp[1] + 3*gp[0] + 3*gp[2] - 2))/2 + 3*gp[0]*((9*gp[1])/2 + (9*gp[0])/2 + (9*gp[2])/2 - 9/2);
          break;
        case 5:
          fkt = -(9*gp[0]*(3*gp[0] - 1))/2 - (3*gp[0] - 1)*((9*gp[1])/2 + (9*gp[0])/2 + (9*gp[2])/2 - 9/2) - 3*gp[0]*((9*gp[1])/2 + (9*gp[0])/2 + (9*gp[2])/2 - 9/2);
          break;
        case 6:
          fkt = (27*gp[1]*gp[0])/2 + (9*gp[1]*(3*gp[0] - 1))/2;
          break;
        case 7:
          fkt = (9*gp[1]*(3*gp[1] - 1))/2;
          break;
        case 8:
          fkt = -(9*gp[1]*(3*gp[1] - 1))/2;
          break;
        case 9:
          fkt = (9*gp[1]*(3*gp[1] + 3*gp[0] + 3*gp[2] - 2))/2 + 3*gp[1]*((9*gp[1])/2 + (9*gp[0])/2 + (9*gp[2])/2 - 9/2);
          break;
        case 10:
          fkt = (9*gp[2]*(3*gp[1] + 3*gp[0] + 3*gp[2] - 2))/2 + 3*gp[2]*((9*gp[1])/2 + (9*gp[0])/2 + (9*gp[2])/2 - 9/2);
          break;
        case 11:
          fkt = (9*gp[2]*(3*gp[0] - 1))/2 + (27*gp[0]*gp[2])/2;
          break;
        case 12:
          fkt = 0.;
          break;
        case 13:
          fkt = -(9*gp[2]*(3*gp[2] - 1))/2;
          break;
        case 14:
          fkt = (9*gp[2]*(3*gp[2] - 1))/2;
          break;
        case 15:
          fkt = 0.;
          break;
      }
      break;

     case N_eta:
      switch (k)
      {
        case 0:
          fkt = -((3*gp[1] + 3*gp[0] + 3*gp[2] - 1)*(3*gp[1] + 3*gp[0] + 3*gp[2] - 2))/2 - 3*(3*gp[1] + 3*gp[0] + 3*gp[2] - 1)*(gp[1]/2 + gp[0]/2 + gp[2]/2 - 1/2) - 3*(3*gp[1] + 3*gp[0] + 3*gp[2] - 2)*(gp[1]/2 + gp[0]/2 + gp[2]/2 - 1/2);
          break;
        case 1:
          fkt = 0.;
          break;
        case 2:
          fkt = ((3*gp[1] - 1)*(3*gp[1] - 2))/2 + (3*gp[1]*(3*gp[1] - 1))/2 + (3*gp[1]*(3*gp[1] - 2))/2;
          break;
        case 3:
          fkt = 0.;
          break;
        case 4:
          fkt = (9*gp[0]*(3*gp[1] + 3*gp[0] + 3*gp[2] - 2))/2 + 3*gp[0]*((9*gp[1])/2 + (9*gp[0])/2 + (9*gp[2])/2 - 9/2);
          break;
        case 5:
          fkt = -(9*gp[0]*(3*gp[0] - 1))/2;
          break;
        case 6:
          fkt = (9*gp[0]*(3*gp[0] - 1))/2;
          break;
        case 7:
          fkt = (27*gp[1]*gp[0])/2 + (9*gp[0]*(3*gp[1] - 1))/2;
          break;
        case 8:
          fkt = -(3*gp[1] - 1)*((9*gp[1])/2 + (9*gp[0])/2 + (9*gp[2])/2 - 9/2) - 3*gp[1]*((9*gp[1])/2 + (9*gp[0])/2 + (9*gp[2])/2 - 9/2) - (9*gp[1]*(3*gp[1] - 1))/2;
          break;
        case 9:
          fkt = (3*gp[1] + 3*gp[0] + 3*gp[2] - 2)*((9*gp[1])/2 + (9*gp[0])/2 + (9*gp[2])/2 - 9/2) + (9*gp[1]*(3*gp[1] + 3*gp[0] + 3*gp[2] - 2))/2 + 3*gp[1]*((9*gp[1])/2 + (9*gp[0])/2 + (9*gp[2])/2 - 9/2);
          break;
        case 10:
          fkt = (9*gp[2]*(3*gp[1] + 3*gp[0] + 3*gp[2] - 2))/2 + 3*gp[2]*((9*gp[1])/2 + (9*gp[0])/2 + (9*gp[2])/2 - 9/2);
          break;
        case 11:
          fkt = 0.;
          break;
        case 12:
          fkt = (27*gp[1]*gp[2])/2 + (9*gp[2]*(3*gp[1] - 1))/2;
          break;
        case 13:
          fkt = -(9*gp[2]*(3*gp[2] - 1))/2;
          break;
        case 14:
          fkt = 0.;
          break;
        case 15:
          fkt = (9*gp[2]*(3*gp[2] - 1))/2;
          break;
      }
      break;

  case N_zeta:
      switch (k)
      {
        case 0:
          fkt = -((3*gp[1] + 3*gp[0] + 3*gp[2] - 1)*(3*gp[1] + 3*gp[0] + 3*gp[2] - 2))/2 - 3*(3*gp[1] + 3*gp[0] + 3*gp[2] - 1)*(gp[1]/2 + gp[0]/2 + gp[2]/2 - 1/2) - 3*(3*gp[1] + 3*gp[0] + 3*gp[2] - 2)*(gp[1]/2 + gp[0]/2 + gp[2]/2 - 1/2);
          break;
        case 1:
          fkt = 0.;
          break;
        case 2:
          fkt = 0.;
          break;
        case 3:
          fkt = (3*gp[2]*(3*gp[2] - 1))/2 + (3*gp[2]*(3*gp[2] - 2))/2 + ((3*gp[2] - 1)*(3*gp[2] - 2))/2;
          break;
        case 4:
          fkt = (9*gp[0]*(3*gp[1] + 3*gp[0] + 3*gp[2] - 2))/2 + 3*gp[0]*((9*gp[1])/2 + (9*gp[0])/2 + (9*gp[2])/2 - 9/2);
          break;
        case 5:
          fkt = -(9*gp[0]*(3*gp[0] - 1))/2;
          break;
        case 6:
          fkt = 0.;
          break;
        case 7:
          fkt = 0.;
          break;
        case 8:
          fkt = -(9*gp[1]*(3*gp[1] - 1))/2;
          break;
        case 9:
          fkt = (9*gp[1]*(3*gp[1] + 3*gp[0] + 3*gp[2] - 2))/2 + 3*gp[1]*((9*gp[1])/2 + (9*gp[0])/2 + (9*gp[2])/2 - 9/2);
          break;
        case 10:
          fkt = (3*gp[1] + 3*gp[0] + 3*gp[2] - 2)*((9*gp[1])/2 + (9*gp[0])/2 + (9*gp[2])/2 - 9/2) + (9*gp[2]*(3*gp[1] + 3*gp[0] + 3*gp[2] - 2))/2 + 3*gp[2]*((9*gp[1])/2 + (9*gp[0])/2 + (9*gp[2])/2 - 9/2);
          break;
        case 11:
          fkt = (9*gp[0]*(3*gp[0] - 1))/2;
          break;
        case 12:
          fkt = (9*gp[1]*(3*gp[1] - 1))/2;
          break;
        case 13:
          fkt = -(9*gp[2]*(3*gp[2] - 1))/2 - (3*gp[2] - 1)*((9*gp[1])/2 + (9*gp[0])/2 + (9*gp[2])/2 - 9/2) - 3*gp[2]*((9*gp[1])/2 + (9*gp[0])/2 + (9*gp[2])/2 - 9/2);
          break;
        case 14:
          fkt = (9*gp[0]*(3*gp[2] - 1))/2 + (27*gp[0]*gp[2])/2;
          break;
        case 15:
          fkt = (27*gp[1]*gp[2])/2 + (9*gp[1]*(3*gp[2] - 1))/2;
          break;
      }
      break;*/
case N_xi:
      switch (k)
      {
        case 0:
          fkt = -0.5*(27.*(gp[0]*gp[0]) - 18.*gp[0] + 2.);
          break;
        case 1:
          fkt = 0.5*(27.*(gp[1]*gp[1]) - 18.*gp[1] + 2.);
          break;
        case 2:
          fkt = 0.;
          break;
        case 3:
          fkt = 0.;
          break;
        case 4:
          fkt = 4.5*(3.*(gp[0]*gp[0]) - gp[0] - 6.*(gp[0]*gp[1]) + gp[1]);
          break;
        case 5:
          fkt = 4.5*(-3.*(gp[1]*gp[1]) + gp[1] + 6.*(gp[0]*gp[1]) - gp[0]);
          break;
        case 6:
          fkt = 4.5*(6.*(gp[1]*gp[2]) - gp[2]);
          break;
        case 7:
          fkt = 4.5*(3.*(gp[2]*gp[2]) - gp[2]);
          break;
        case 8:
          fkt = -4.5*(3.*(gp[2]*gp[2]) - gp[2]);
          break;
        case 9:
          fkt = -4.5*(6.*(gp[0]*gp[2]) - gp[2]);
          break;
        case 10:
          fkt = -4.5*(6.*(gp[0]*gp[3]) - gp[3]);
          break;
        case 11:
          fkt = 4.5*(6.*(gp[1]*gp[3]) - gp[3]);;
          break;
        case 12:
          fkt = 0.;
          break;
        case 13:
          fkt = -4.5*(3.*(gp[3]*gp[3]) - gp[3]);
          break;
        case 14:
          fkt = 4.5*(3.*(gp[3]*gp[3]) - gp[3]);
          break;
        case 15:
          fkt = 0.;
          break;
      }
      break;

     case N_eta:
      switch (k)
      {
        case 0:
          fkt = -0.5*(27.*(gp[0]*gp[0]) - 18.*gp[0] + 2.);
          break;
        case 1:
          fkt = 0.;
          break;
        case 2:
          fkt = 0.5*(27.*(gp[2]*gp[2]) - 18.*gp[2] + 2.);
          break;
        case 3:
          fkt = 0.;
          break;
        case 4:
          fkt = -4.5*(6.*(gp[0]*gp[1]) - gp[1]);
          break;
        case 5:
          fkt = -4.5*(3.*(gp[1]*gp[1]) - gp[1]);
          break;
        case 6:
          fkt = 4.5*(3.*(gp[1]*gp[1]) - gp[1]);
          break;
        case 7:
          fkt = 4.5*(6.*(gp[1]*gp[2]) - gp[1]);
          break;
        case 8:
          fkt = 4.5*(-3.*(gp[2]*gp[2]) + gp[2] + 6.*(gp[0]*gp[2]) - gp[0]);
          break;
        case 9:
          fkt = 4.5*(3.*(gp[0]*gp[0]) - gp[0] - 6.*(gp[0]*gp[2]) + gp[2]);
          break;
        case 10:
          fkt = -4.5*(6.*(gp[0]*gp[3]) - gp[3]);
          break;
        case 11:
          fkt = 0.;
          break;
        case 12:
          fkt = 4.5*(6.*(gp[2]*gp[3]) - gp[3]);
          break;
        case 13:
          fkt = -4.5*(3.*(gp[3]*gp[3]) - gp[3]);
          break;
        case 14:
          fkt = 0.;
          break;
        case 15:
          fkt = 4.5*(3.*(gp[3]*gp[3]) - gp[3]);
          break;
      }
      break;

  case N_zeta:
      switch (k)
      {
        case 0:
          fkt = -0.5*(27.*(gp[0]*gp[0]) - 18.*gp[0] + 2.);
          break;
        case 1:
          fkt = 0.;
          break;
        case 2:
          fkt = 0.;
          break;
        case 3:
          fkt = 0.5*(27.*(gp[3]*gp[3]) - 18.*gp[3] + 2.);
          break;
        case 4:
          fkt = -4.5*(6.*(gp[0]*gp[1]) - gp[1]);
          break;
        case 5:
          fkt = -4.5*(3.*(gp[1]*gp[1]) - gp[1]);
          break;
        case 6:
          fkt = 0.;
          break;
        case 7:
          fkt = 0.;
          break;
        case 8:
          fkt = -4.5*(3.*(gp[2]*gp[2]) - gp[2]);
          break;
        case 9:
          fkt = -4.5*(6.*(gp[0]*gp[2]) - gp[2]);
          break;
        case 10:
          fkt = 4.5*(3.*(gp[0]*gp[0]) - gp[0] - 6.*(gp[0]*gp[3]) + gp[3]);
          break;
        case 11:
          fkt = 4.5*(3.*(gp[1]*gp[1]) - gp[1]);
          break;
        case 12:
          fkt = 4.5*(3.*(gp[2]*gp[2]) - gp[2]);
          break;
        case 13:
          fkt = 4.5*(-3.*(gp[3]*gp[3]) + gp[3] + 6.*(gp[0]*gp[3]) - gp[0]);
          break;
        case 14:
          fkt = 4.5*(6.*(gp[1]*gp[3]) - gp[1]);
          break;
        case 15:
          fkt = 4.5*(6.*(gp[2]*gp[3]) - gp[2]);
          break;
      }
      break;


 /*case N_L1:
      switch (k)
      {
        case 0:
          fkt = 0.5*(27.*(gp[0]*gp[0]) - 18.*gp[0] + 2.);
          break;
        case 1:
          fkt = 0.;
          break;
        case 2:
          fkt = 0.;
          break;
        case 3:
          fkt = 0.;
          break;
        case 4:
          fkt = 4.5*(6.*gp[0]*gp[1] - gp[1]);
          break;
        case 5:
          fkt = 4.5*(3.*(gp[1]*gp[1]) - gp[1]);
          break;
        case 6:
          fkt = 0.;
          break;
        case 7:
          fkt = 0.;
          break;
        case 8:
          fkt = 4.5*(3.*(gp[2]*gp[2]) - gp[2]);
          break;
        case 9:
          fkt = 4.5*(6.*(gp[0]*gp[2]) - gp[2]);
          break;
        case 10:
          fkt = 4.5*(6.*(gp[0]*gp[3]) - gp[3]);
          break;
        case 11:
          fkt = 0.;
          break;
        case 12:
          fkt = 0.;
          break;
        case 13:
          fkt = 4.5*(3.*(gp[3]*gp[3]) - gp[3]);
          break;
        case 14:
          fkt = 0.;
          break;
        case 15:
          fkt = 0.;
          break;
      }
      break;
   case N_L2:
      switch (k)
      {
        case 0:
          fkt = 0.;
          break;
        case 1:
          fkt = 0.5*(27.*(gp[1]*gp[1]) - 18.*gp[1] + 2.);;
          break;
        case 2:
          fkt = 0.;
          break;
        case 3:
          fkt = 0.;
          break;
        case 4:
          fkt = 4.5*(3.*(gp[0]*gp[0]) - gp[0]);
          break;
        case 5:
          fkt = 4.5*(6.*(gp[0]*gp[1]) - gp[0]);
          break;
        case 6:
          fkt = 4.5*(6.*(gp[1]*gp[2]) - gp[2]);
          break;
        case 7:
          fkt = 4.5*(3.*(gp[2]*gp[2]) - gp[2]);
          break;
        case 8:
          fkt = 0.;
          break;
        case 9:
          fkt = 0.;
          break;
        case 10:
          fkt = 0.;
          break;
        case 11:
          fkt = 4.5*(6.*(gp[1]*gp[3]) - gp[3]);
          break;
        case 12:
          fkt = 0.;
          break;
        case 13:
          fkt = 0.;
          break;
        case 14:
          fkt = 4.5*(3.*(gp[3]*gp[3]) - gp[3]);
          break;
        case 15:
          fkt = 0.;
          break;
      }
      break;
   case N_L3:
      switch (k)
      {
        case 0:
          fkt = 0.;
          break;
        case 1:
          fkt = 0.;
          break;
        case 2:
          fkt = 0.5*(27.*(gp[2]*gp[2]) - 18.*gp[2] + 2.);
          break;
        case 3:
          fkt = 0.;
          break;
        case 4:
          fkt = 0.;
          break;
        case 5:
          fkt = 0.;
          break;
        case 6:
          fkt = 4.5*(3.*(gp[1]*gp[1]) - gp[1]);
          break;
        case 7:
          fkt = 4.5*(6.*(gp[1]*gp[2]) - gp[1]);
          break;
        case 8:
          fkt = 4.5*(6.*(gp[0]*gp[2]) - gp[0]);
          break;
        case 9:
          fkt = 4.5*(3.*(gp[0]*gp[0]) - gp[0]);
          break;
        case 10:
          fkt = 0.;
          break;
        case 11:
          fkt = 0.;
          break;
        case 12:
          fkt = 4.5*(6.*(gp[2]*gp[3]) - gp[3]);
          break;
        case 13:
          fkt = 0.;
          break;
        case 14:
          fkt = 0.;
          break;
        case 15:
          fkt = 4.5*(3.*(gp[3]*gp[3]) - gp[3]);
          break;
      }
      break;
    case N_L4:
      switch (k)
      {
        case 0:
          fkt = 0.;
          break;
        case 1:
          fkt = 0.;
          break;
        case 2:
          fkt = 0.;
          break;
        case 3:
          fkt = 0.5*(27.*(gp[3]*gp[3]) - 18.*gp[3] + 2.);
          break;
        case 4:
          fkt = 0.;
          break;
        case 5:
          fkt = 0.;
          break;
        case 6:
          fkt = 0.;
          break;
        case 7:
          fkt = 0.;
          break;
        case 8:
          fkt = 0.;
          break;
        case 9:
          fkt = 0.;
          break;
        case 10:
          fkt = 4.5*(3.*(gp[0]*gp[0]) - gp[0]);
          break;
        case 11:
          fkt = 4.5*(3.*(gp[1]*gp[1]) - gp[1]);
          break;
        case 12:
          fkt = 4.5*(3.*(gp[2]*gp[2]) - gp[2]);
          break;
        case 13:
          fkt = 4.5*(6.*(gp[0]*gp[3]) - gp[0]);
          break;
        case 14:
          fkt = 4.5*(6.*(gp[1]*gp[3]) - gp[1]);
          break;
        case 15:
          fkt = 4.5*(6.*(gp[2]*gp[3]) - gp[2]);
          break;
      }
      break;*/
      default:
      throw cException(" **** invalid request for derivative **** ", __FILE__, __LINE__);
      break;
  }
  return fkt;
}


PetscReal cShapeFunctions::BeamN(
  const int &nnod,
  const int &k,
  const eDerivative &deriv, 
  const cPoint &gp)
{
  switch (nnod)
  {
    case 2:
      return Beam2(k, deriv, gp);
      break;
    case 3:
      return Beam3(k, deriv, gp);
      break;
    default:
      break;
  }

  return 0.0;
}

// ----------------------------------------------------------------------------
//
//       -o----+----o--->*xi
//       0           1
//
// ----------------------------------------------------------------------------
PetscReal cShapeFunctions::Beam2(
  const int &k,
  const eDerivative &deriv,
  const cPoint &gp)
{
  double fkt = 0.0;

  switch(deriv)
  {
    case N_fun:
      switch(k)
      {
        case 0:
          fkt = 0.5 * (1.0 - gp[0]);
          break;
        case 1:
          fkt = 0.5 * (1.0 + gp[0]);
          break;
      }
      break;

    case N_xi:
      switch(k)
      {
        case 0:
          fkt = -0.5;
          break;
        case 1:
          fkt =  0.5;
          break;
      }
      break;

    default:
      throw cException(" **** invalid request for derivative **** ", __FILE__, __LINE__);
      break;
  }

  return fkt;
}


/*BEGIN_NO_COVERAGE*/
// ----------------------------------------------------------------------------
//
//       -o----+----o--->*xi
//       0     1     2
//
// ----------------------------------------------------------------------------
PetscReal cShapeFunctions::Beam3(
  const int &k,
  const eDerivative &deriv,
  const cPoint &gp)
{
  PetscReal fkt = 0.;

  throw cException("cShapeFunctions::Beam3 noch nicht implementiert", __FILE__, __LINE__);

  return fkt;
}
/*END_NO_COVERAGE*/

