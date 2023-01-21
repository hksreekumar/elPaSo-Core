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

#ifndef INFAM_GAUSSPOINTS_H
#define INFAM_GAUSSPOINTS_H

#include "point.h"

namespace infam {

//! integration points for 1D
const double gp4[10] = {0.0,  // 1
                        0.577350269189626,
                        -0.577350269189626,  // 2
                        0.774596669241483,
                        -0.000000000000000,
                        -0.774596669241483,  // 3
                        0.861136311594053,
                        0.339981043584856,
                        -0.339981043584856,  // 4
                        -0.861136311594053};

//! integration weights for 1D
const double gw4[10] = {
    2.000000000000000,                                        // 1
    1.000000000000000, 1.000000000000000,                     // 2
    0.555555555555556, 0.888888888888889, 0.555555555555556,  // 3
    0.347854845137454, 0.652145154862546, 0.652145154862546,  // 4
    0.347854845137454};

// ---------------------------------------------------------------------------
//   Stuetzstellen und Gewichte fuer Tria-Elemente
// ---------------------------------------------------------------------------
const int anztria = 3;            // Anzahl der Integrationslevel
const int arten[3] = {1, 3, 4};   // Anzahl der Gauss-Pkt.
const int sttria[3] = {0, 1, 4};  // Start-Position in Array

//! integration points used for triangles as well as tetrahedrons
const double gp3[24] = {
    0.3333333333333330,   0.3333333333333330, 0.3333333333333330,  // 1

    0.500000000000000000, 0.5000000000000000, 0.0000000000000000,  // 3
    0.500000000000000000, 0.0000000000000000, 0.5000000000000000,
    0.000000000000000000, 0.5000000000000000, 0.5000000000000000,

    0.3333333333333330,   0.3333333333333330, 0.3333333333333340,  // 4
    0.6000000000000000,   0.2000000000000000, 0.2000000000000000,
    0.2000000000000000,   0.6000000000000000, 0.2000000000000000,
    0.2000000000000000,   0.2000000000000000, 0.6000000000000000,
};

//! integration weights used for triangles as well as tetrahedrons
const double gw3[8] = {
    1.0000000000000000,  // 1

    0.3333333333333329,  0.3333333333333329, 0.3333333333333329,  // 3

    -0.5625000000000000,  // 4
    0.5208333333333330,  0.5208333333333330, 0.5208333333333330,
};

// Zienkiewicz, FEM Its basis & fundamentals 6th edition, page 166
// mech I 057
//! integration points used for tetrahedrons in xi,eta,zeta coordinates system
const int anztetra = 6;                        // Anzahl der Integrationslevel
const int anzgp[6] = {1, 4, 5, 6, 7, 8};       // Anzahl der Gauss-Pkt.
const int sttetra[6] = {0, 1, 5, 30, 34, 50};  // Start-Position in Array

const double gp3t[70] = {
    0.2500000000000000,
    0.2500000000000000,
    0.2500000000000000,  // 1

    0.5854102000000000,
    0.1381960000000000,
    0.1381960000000000,  // 2
    0.1381960000000000,
    0.5854102000000000,
    0.1381960000000000,
    0.1381960000000000,
    0.1381960000000000,
    0.5854102000000000,
    0.1381960000000000,
    0.1381960000000000,
    0.1381960000000000,

    0.2500000000000000,
    0.2500000000000000,
    0.2500000000000000,  // 3
    0.1666666666666667,
    0.1666666666666667,
    0.1666666666666667,
    0.5000000000000000,
    0.1666666666666667,
    0.1666666666666667,
    0.1666666666666667,
    0.5000000000000000,
    0.1666666666666667,
    0.1666666666666667,
    0.1666666666666667,
    0.5000000000000000,

    0.2500000000000000,
    0.2500000000000000,
    0.2500000000000000,
    0.2500000000000000,  // 1

    0.5854102000000000,
    0.1381960000000000,
    0.1381960000000000,
    0.1381960000000000,  // 2
    0.1381960000000000,
    0.5854102000000000,
    0.1381960000000000,
    0.1381960000000000,
    0.1381960000000000,
    0.1381960000000000,
    0.5854102000000000,
    0.1381960000000000,
    0.1381960000000000,
    0.1381960000000000,
    0.1381960000000000,
    0.5854102000000000,

    // 3
    0.2500000000000000,
    0.2500000000000000,
    0.2500000000000000,
    0.2500000000000000,
    0.5000000000000000,
    0.1666666666666667,
    0.1666666666666667,
    0.1666666666666667,
    0.1666666666666667,
    0.5000000000000000,
    0.1666666666666667,
    0.1666666666666667,
    0.1666666666666667,
    0.1666666666666667,
    0.5000000000000000,
    0.1666666666666667,
    0.1666666666666667,
    0.1666666666666667,
    0.1666666666666667,
    0.5000000000000000,

};

//! integration weights used for tetrahedrons

const double gw3t[20] = {
    // integration weights for gauss points in xi,eta,zeta coordinate system

    1., 0.250000000000000, 0.250000000000000, 0.250000000000000,
    0.250000000000000,  // 2
    -0.800000000000000, 0.450000000000000, 0.450000000000000, 0.450000000000000,
    0.450000000000000,  // 3
    // integration weights for gauss points in L1,L2,L3,L4 coordinate system
    1.0, 0.250000000000000, 0.250000000000000, 0.250000000000000,
    0.250000000000000,  // 2
    -0.80000000000000, 0.450000000000000, 0.450000000000000, 0.450000000000000,
    0.450000000000000,  // 3

};

}  // namespace infam

//! @brief stores the Gauss points used for numerical integration
//! @author Dirk Clasen
//! @date 25.04.2005
//! This class provides the Gauss points as well as their weights
//! used for numerical integration of the elementmatrices.
class cGaussPoints {
 private:
  //! returns the first index of Gauss point coordinates of the specified
  //! integrationlevel within the arrays gpt4/gwt4
  //! @param level integrationlevel
  //! @return first index within array
  inline int Start(int level) const;

  //! Gauss points for hexahedrons
  //! @param level number of Gauss points (per direction)
  //! @param x [0..level-1]
  //! @param y [0..level-1]
  //! @param z [0..level-1]
  //! @return soughtfor Gauss point
  cPoint getGaussPoint3D(int level, int x, int y, int z) const;

  //! Gauss points for tetrahedrons
  //! @param level number of Gauss points (per direction)
  //! @param x [0..level-1]
  //! @param y [0..level-1]
  //! @param z [0..level-1]
  //! @return soughtfor Gauss point
  // cPoint getGaussPointTet3D(int level, int x, int y, int z) const;

  //! integration weights for hexahedrons
  //! @param level number of Gauss points (per direction)
  //! @param x [0..level-1]
  //! @param y [0..level-1]
  //! @param z [0..level-1]
  //! @return weight
  PetscReal getGaussWeight3D(int level, int x, int y, int z) const;

  //! integration weights for tetrahedrons
  //! @param level number of Gauss points (per direction)
  //! @param x [0..level-1]
  //! @param y [0..level-1]
  //! @param z [0..level-1]
  //! @return weight
  // PetscReal getGaussWeightTet3D(int level, int x, int y, int z) const;

  //! Gauss points for quadrilaterals
  //! @param level number of Gauss points (per direction)
  //! @param x [0..level-1]
  //! @param y [0..level-1]
  //! @return soughtfor Gauss point
  cPoint getGaussPoint2D(int level, int x, int y) const;

  //! integration weights for quadrilaterals
  //! @param level number of Gauss points (per direction)
  //! @param x [0..level-1]
  //! @param y [0..level-1]
  //! @return weight
  PetscReal getGaussWeight2D(int level, int x, int y) const;

 public:
  cGaussPoints();
  ~cGaussPoints();

  //! Gauss points used for numerical integration of hexahedrons.
  //! These points are numbered sequentially, their "coordinates" are
  //! determined within this routine by a private member of this class.
  //! @param level integration level
  //! @param nr number of current point
  //! @return soughtfor Gauss point
  cPoint getGaussPoint3D(int level, int nr) const;

  //! Integration weights used for numerical integration of hexahedrons.
  //! These points are numbered sequentially, their "coordinates" are
  //! determined within this routine by a private member of this class.
  //! @param level integration level
  //! @param nr number of current point
  //! @return soughtfor integration weights
  PetscReal getGaussWeight3D(int level, int nr) const;

  //! Gauss points used for numerical integration of tetrahedrons.
  //! These points are numbered sequentially, their "coordinates" are
  //! determined within this routine by a private member of this class.
  //! @param level integration level
  //! @param nr number of current point
  //! @return soughtfor Gauss point
  cPoint getGaussPointTet3D(int level, int nr) const;

  //! Gauss points used for numerical integration of tetrahedrons in L1,L2,L3,L4
  //! coordinates. These points are numbered sequentially, their "coordinates"
  //! are determined within this routine by a private member of this class.
  //! @param level integration level
  //! @param nr number of current point
  //! @return soughtfor Gauss point
  cPoint getGaussPointTet3DL(int level, int nr) const;

  //! Integration weights used for numerical integration of tetrahedrons.
  //! These points are numbered sequentially, their "coordinates" are
  //! determined within this routine by a private member of this class.
  //! @param level integration level
  //! @param nr number of current point
  //! @return soughtfor integration weights
  PetscReal getGaussWeightTet3D(int level, int nr) const;

  //! Integration weights used for numerical integration of tetrahedrons in
  //! L1,L2,L3,L4 coordinates. These points are numbered sequentially, their
  //! "coordinates" are determined within this routine by a private member of
  //! this class.
  //! @param level integration level
  //! @param nr number of current point
  //! @return soughtfor integration weights
  PetscReal getGaussWeightTet3DL(int level, int nr) const;

  //! Gauss points used for numerical integration of quadrilaterals.
  //! These points are numbered sequentially, their "coordinates" are
  //! determined within this routine by a private member of this class.
  //! @param level integration level
  //! @param nr number of current point
  //! @return soughtfor Gauss point
  cPoint getGaussPoint2D(int level, int nr) const;

  //! Integration weights used for numerical integration of quadrilaterals.
  //! These points are numbered sequentially, their "coordinates" are
  //! determined within this routine by a private member of this class.
  //! @param level integration level
  //! @param nr number of current point
  //! @return soughtfor integration weights
  PetscReal getGaussWeight2D(int level, int nr) const;

  cPoint getGaussPoint1D(int level, int nr) const;

  PetscReal getGaussWeight1D(int level, int nr) const;

  //! Gauss points used for numerical integration of triangles.
  //! These points are numbered sequentially, their "coordinates" are
  //! determined within this routine by a private member of this class.
  //! @param level integration level
  //! @param nr number of current point
  //! @return soughtfor Gauss point
  cPoint getGaussPointTria(int level, int nr) const;

  //! Integration weights used for numerical integration of triangles.
  //! These points are numbered sequentially, their "coordinates" are
  //! determined within this routine by a private member of this class.
  //! @param level integration level
  //! @param nr number of current point
  //! @return soughtfor integration weights
  PetscReal getGaussWeightTria(int level, int nr) const;

  //! Gauss points used for numerical integration of tetrahedron.
  //! These points are numbered sequentially, their "coordinates" are
  //! determined within this routine by a private member of this class.
  //! @param level integration level
  //! @param nr number of current point
  //! @return soughtfor Gauss point
  // cPoint getGaussPointTetra(int level, int nr) const;

  //! Integration weights used for numerical integration of tetrahedron.
  //! These points are numbered sequentially, their "coordinates" are
  //! determined within this routine by a private member of this class.
  //! @param level integration level
  //! @param nr number of current point
  //! @return soughtfor integration weights
  // PetscReal getGaussWeightTetra(int level, int nr) const;
};

#endif
