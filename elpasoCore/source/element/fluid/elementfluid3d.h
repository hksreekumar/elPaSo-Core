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

#ifndef INFAM_ELEMENT_FLUID_3D_H
#define INFAM_ELEMENT_FLUID_3D_H

#include "elementfluid.h"

//!  @brief base class for 3d fluid elements
//!  @author Dirk Clasen
//!  @date 25.04.2007
class cElementFluid3d : public cElementFluid,
                        private tCounter<cElementFluid3d> {
 private:
  //! @brief computes the surface integral on a face of a hexaedron. The result
  //! is obtained in global coordinates and has to be transformed to local
  //! coordinate system
  //! @param Face surface on which the integral will be computed. valid range:
  //! [0..5]
  //! @param Nface shape functions of element's face
  //! @param C results of the integration
  void evaluateSurfaceIntegral(int Face, cArray3d &Nface, cMatrix &C);

  //! @brief evaluate impedance boundary condition
  //! This function will simply call subEvaluateImpedance
  //! @param omega angular frequency in s^{-1}
  //! @param KM    dynamic stiffnessmatrix
  void evaluateImpedance(const PetscReal &omega, cElementMatrix &KM);

  //! @brief evaluated shapefunctions for element's faces
  virtual cArray3d getShapeFunctionsFace(void) const = 0;

  //! @brief get the GaussPoints for the element's faces
  void getGaussPointsFace(int Face, std::vector<cPoint> &eps, int &ngps,
                          std::vector<PetscReal> &weights, const int nnod_face);

 protected:
  //! @brief return the number of nodes that describe one face of the element
  virtual short getNumberOfNodesPerFace(void) const = 0;

 public:
  //! @brief Constructor
  cElementFluid3d(short NumberOfNodes, short NumberOfGaussPoints);
  //! @brief Copy constructor
  cElementFluid3d(const cElementFluid3d &other);
  //! @brief Destructor
  virtual ~cElementFluid3d();

  //! @brief return the number of faces of this elements
  virtual short getNumberOfFaces(void) const = 0;

  //! @brief return the normal vector on one face of the 2d fluid element.
  //! The normal vector points outward of the element.
  cVector getNormalOnFace(const short &face) const;

  //! @brief returns the number of class' instances
  static size_t howMany(void) { return tCounter<cElementFluid3d>::howMany(); }

  //! @brief assembly of element's stiffnessmatrix
  //! @param KM assembled stiffnessmatrix
  //! @param x solution vector
  //! @param dx change of x
  void assembleStiffnessMatrix(cElementMatrix &KM, Vec *x, Vec *dx);

  //! @brief assembly of element's massmatrix
  //! @param KM assembled massmatrix
  void assembleMassMatrix(cElementMatrix &MM);

  //! @brief assembly of element's loadvector
  //! @param LV assembled loadvector
  //! @param KM the assembled stiffnessmatrix
  //! @param x solution vector
  //! @param dx change of x
  void assembleLoadVector(cElementVector &LV, cElementMatrix &KM, Vec *x,
                          Vec *dx);

  //! @brief assemble frequency dependent loadvector
  //! @param omega  angular frequency omega = 2 \pi f [s^{-1}]
  //! @param LV  dynamic elementloadvector
  void assembleDynamicLoadVector(const PetscReal &omega, cElementVector &LV,
                                 cElementMatrix &EM);

  //! @brief Numbering of one element's face
  //! The normal vector points outward.
  //! @param Face number of the face. Valid range: [0..5]
  virtual std::vector<short> getIndicesOfFaceNodes(int Face) const = 0;

  //! @brief evaluate area of a face
  //! @param face number of face
  PetscReal getFaceArea(const short &face);

  //! @brief compute stresses element wise
  //! @param fullSolution the full solution vector, each process holds a copy
  //! @param Stresses a vector that holds the global stresses
  void computeStressesCells(Vec &fullSolution, Vec &stresses,
                            Vec &stressesSec2) {
    trace("ERROR: cElementFluid::computeStressesCells() not implemented, yet.");
    ExitApp();
  }

  //! @brief compute stresses at the nodes
  //! @param fullSolution the full solution vector, each process holds a copy
  //! @param Stresses a vector that holds the global stresses
  void computeStressesNodes(Vec &fullSolution, Vec &stresses);

  //! @brief compute gradient of pore fluids pressure. This will be needed
  //! in order to compute the fluid velocity.
  void computeGradP(cElementVector &gradP, const cElementVector &valP);

  //! @brief compute gradient of pore fluids pressure at a specific point
  void computeGradPPoint(cPoint &ep, cArray3d &Nk, cElementVector &gradPPoint,
                         const cElementVector &valP);

  //! @brief computes the surface integral of the sound power on a face of a
  //! hexaedron based on the reactive sound intensity.
  //! @param fullSolution the full solution vector, each process holds a copy
  //! @param Face surface on which the integral will be computed. valid range:
  //! [0..5]
  //! @param Ps_face_int results of the integration
  void calculateReactiveSoundPowerSurfaceIntegral(Vec fullSolution, int Face,
                                                  PetscReal &Ps_face_int);

  //! @brief computes the surface integral of the sound power on a face of a
  //! hexaedron based on the active sound intensity.
  //! @param fullSolution the full solution vector, each process holds a copy
  //! @param Face surface on which the integral will be computed. valid range:
  //! [0..5]
  //! @param Ps_face_int results of the integration
  void calculateActiveSoundPowerSurfaceIntegral(Vec fullSolution, int Face,
                                                PetscReal &Ps_face_int);
  void calculateActiveSoundPowerSurfaceIntegralGivenVn(Vec fullSolution,
                                                       int Face,
                                                       PetscReal &Ps_face_int,
                                                       PetscReal GivenVn);

  //! @brief computes the average sound power on a face of a hexaedron based on
  //! the reactive sound intensity.
  //! @param fullSolution the full solution vector, each process holds a copy
  //! @param Face surface on which the integral will be computed. valid range:
  //! [0..5]
  //! @param Ps_face_average results of the average sound power
  void calculateReactiveSoundPowerSurfaceAverage(Vec fullSolution, int Face,
                                                 PetscReal &Ps_face_average);

  //! @brief computes the average sound power on a face of a hexaedron based on
  //! the active sound intensity.
  //! @param fullSolution the full solution vector, each process holds a copy
  //! @param Face surface on which the integral will be computed. valid range:
  //! [0..5]
  //! @param Ps_face_average results of the average sound power
  void calculateActiveSoundPowerSurfaceAverage(Vec fullSolution, int Face,
                                               PetscReal &Ps_face_average);

  //! @brief computes the energy density and the intensity at the nodes
  void computeEDandIntNodes(Vec FullSolution, Vec &EnergyDensity,
                            Vec &intensity);

  //! computes the phase shift of pressure and velocity at the nodes
  void computePhaseShiftNodes(Vec FullSolution, Vec &phaseShiftAverageVec,
                              Vec &phaseShiftVec);

  //! @brief writes this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &write(std::ostream &os) const = 0;

  //! @brief writes this object in XML
  //! @param os outputstream
  //! @return der modified outputstream
  virtual std::ostream &writeXml(std::ostream &os) const = 0;
};

//! @brief overloaded output operator
inline std::ostream &operator<<(std::ostream &os,
                                const cElementFluid3d &other) {
  return other.write(os);
}

#endif
