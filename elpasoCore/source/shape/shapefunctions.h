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

#ifndef INFAM_SHAPEFUNCTION_H
#define INFAM_SHAPEFUNCTION_H

#include "../misc/point.h"
#include "../misc/log/logging.h"


/// known derivatives of the functions
enum eDerivative
{
  N_fun    =  0,    /*!< N(xsi,eta) */
  N_xi     =  1,    /*!< dN/dxsi */
  N_eta    =  2,    /*!< dN/deta */
  N_zeta   =  3,    /*!< dN/dzeta */
  N_L1     =  4,    /*!< dN/dL1 */
  N_L2     =  5,    /*!< dN/dL2 */
  N_L3     =  6,    /*!< dN/dL3 */
  N_L4     =  7,    /*!< dN/dL4 */
  N_xi2    =  8,    /*!< dN/dxsi */
  N_eta2   =  9,    /*!< dN/deta */
  N_xi_eta = 10,    /*!< dN/dxi,deta */
};


/**
 * @brief Shape/testfunctions
 * @author Dirk Clasen
 * @date 28.09.2004
 *
 * This class provides shape functions used to integrate the element matrices.
 */
class cShapeFunctions :
  public virtual cLogging
{
protected:

  PetscReal LagrangeQuadN(const int &nnod, const int &k, const eDerivative &deriv,	const cPoint &gp);
  PetscReal LagrangeQuad4(const int &k,	const eDerivative &deriv,	const cPoint &gp);
  PetscReal LagrangeQuad9(const int &k,	const eDerivative &deriv,	const cPoint &gp);

  PetscReal HermiteQuadN(const int &nnod, const int &k, const eDerivative &deriv,	const cPoint &gp);
  PetscReal HermiteQuad4(const int &k,	const eDerivative &deriv,	const cPoint &gp);

  PetscReal SerendipityQuadN(const int &nnod, const int &k, const eDerivative &deriv, const cPoint &gp);
  PetscReal SerendipityQuad8(const int &k, const eDerivative &deriv, const cPoint &gp);

  PetscReal LagrangeHexN(const int &nnod, const int &k, const eDerivative &deriv, const cPoint &gp);
  PetscReal LagrangeHex8(const int &k, const eDerivative &deriv, const cPoint &gp);
  PetscReal LagrangeHex20(const int &k, const eDerivative &deriv, const cPoint &gp);
  PetscReal LagrangeHex27(const int &k,	const eDerivative &deriv,	const cPoint &gp);

  PetscReal TriaN(const int &nnod, const int &k, const eDerivative &deriv, const cPoint &gp);
  PetscReal Tria3(const int &k, const eDerivative &deriv, const cPoint &gp);
  PetscReal Tria6(const int &k, const eDerivative &deriv, const cPoint &gp);
  PetscReal Tria9(const int &k, const eDerivative &deriv, const cPoint &gp);

  PetscReal TetN(const int &nnod, const int &k, const eDerivative &deriv, const cPoint &gp);
  PetscReal Tet4(const int &k, const eDerivative &deriv, const cPoint &gp);
  PetscReal Tet4L(const int &k, const eDerivative &deriv, const cPoint &gp);
  PetscReal Tet10(const int &k, const eDerivative &deriv, const cPoint &gp);
  PetscReal Tet10L(const int &k, const eDerivative &deriv, const cPoint &gp);
  PetscReal Tet16(const int &k, const eDerivative &deriv, const cPoint &gp);
  PetscReal Tet16L(const int &k, const eDerivative &deriv, const cPoint &gp);

  PetscReal BeamN(const int &nnod, const int &k, const eDerivative &deriv, const cPoint &gp);
  PetscReal Beam2(const int &k, const eDerivative &deriv, const cPoint &gp);
  PetscReal Beam3(const int &k, const eDerivative &deriv, const cPoint &gp);


public:
  cShapeFunctions();
  ~cShapeFunctions();

  /**
   * Evaluate a single shape-/testfunction
   * @param shape shape of the element (Tria, Quadrilateral, Hexahedron, Tetrahedron)
   * @param nnod  number of nodes within element
   * @param k     evaluate shapefunction of node k
   * @param deriv denotes the derivative (N_fun, N_xi, N_eta, N_zeta)
   * @param gp    Gauss point for which the shapefunction has to be evaluated
   */
  PetscReal evaluateShapeFunction(
    const eElementShape &shape,
    const int &nnod,
    const int &k,
    const eDerivative &deriv,
    const cPoint &gp);
};

#endif
