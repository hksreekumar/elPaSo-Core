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

#ifndef INFAM_ELEMENT_STRUCTURE_MINDLIN_H
#define INFAM_ELEMENT_STRUCTURE_MINDLIN_H

#include "../elementstructurelinear.h"

//! @brief base class for plate elements based upon Mindlins' theory
//! @author Dirk Clasen
//! @date 09.08.2005
class cElementStructureMindlin : public cElementStructureLinear,
                                 private tCounter<cElementStructureMindlin> {
 public:
  //! @brief Constructor
  cElementStructureMindlin(short NumberOfNodes, short NumberOfGaussPoints);
  //! @brief Copy constructor
  cElementStructureMindlin(const cElementStructureMindlin &other);
  //! @brief Destructor
  virtual ~cElementStructureMindlin();

  //! @brief return the type of element (sometimes better to use this
  //! instead of typeid or dynamic_cast<>
  eElementType getElementType(void) const { return MindlinPlate; }

  //! @brief return the number of instances of this object
  static size_t howMany(void) {
    return tCounter<cElementStructureMindlin>::howMany();
  }

  //! @brief return vector that contains elements degrees of freedom
  std::vector<eKnownDofs> getDofs(void) const;

  //! @brief compute elements' massmatrix
  //! @param KM assembled massmatrix of this element
  void assembleMassMatrix(cElementMatrix &MM);

  //! @brief compute element's loadvector
  //! @param LV element's loadvector
  void assembleLoadVector(cElementVector &LV, cElementMatrix &KM, Vec *x,
                          Vec *dx);

  //! @brief compute stresses element wise
  //! @param fullSolution the full solution vector, each process holds a copy
  //! @param Stresses a vector that holds the global stresses
  void computeStressesCells(Vec &fullSolution, Vec &stresses,
                            Vec &stressesSec2);

  //! @brief get stresses element wise without adding to global stresses (needed
  //! for shell element)
  //! @param fullSolution the full solution vector, each process holds a copy
  //! @param stress2d a vector that holds the local stresses
  //! @param stress2dSec2 a vector that holds the local stresses for section 2
  void getStressesCells(cElementVector &relevantLocalSolutionDSG4,
                        cElementVector &stress2d, cElementVector &stress2dSec2);

  //! @brief compute stresses at the nodes
  //! @param FullSolution the full solution vector, each process holds a copy
  //! @param Stresses a vector that holds the global stresses
  void computeStressesNodes(Vec &fullSolution, Vec &stresses) {
    trace(
        "ERROR: cElementStructureMindlin::computeStressesNodes() not "
        "implemented, yet.");
    ExitApp();
  }

  //! @brief returns local coordinates of element's nodes
  void getLocalCoordinatesOfElementsNodes(cMatrix &xyz) {
    trace("ERROR: getLocalCoordinatesOfElementsNodes() not implemented, yet.");
    ExitApp();
  }

  //! @brief write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &write(std::ostream &os) const = 0;

  //! @brief write this object in XML
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &writeXml(std::ostream &os) const = 0;
};

//! @brief overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cElementStructureMindlin &other) {
  return other.write(os);
}

#endif
