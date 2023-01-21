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

#ifndef INFAM_ELEMENT_STRUCTURE_KIRCHHOFF_H
#define INFAM_ELEMENT_STRUCTURE_KIRCHHOFF_H

#include "../elementstructurelinear.h"

//! @brief base class for Kirchhoff plate elements
//! @author Dirk Clasen
//! @date 24.08.2005
class cElementStructureKirchhoff
    : public cElementStructureLinear,
      private tCounter<cElementStructureKirchhoff> {
 public:
  //! @brief Constructor
  cElementStructureKirchhoff(short NumberOfNodes, short NumberOfGaussPoints);
  //! @brief Copy constructor
  cElementStructureKirchhoff(const cElementStructureKirchhoff &other);
  //! @brief Destructor
  virtual ~cElementStructureKirchhoff();

  //! @brief return the type of element (sometimes better to use this
  //! instead of typeid or dynamic_cast<>
  eElementType getElementType(void) const { return KirchhoffPlate; }

  //! @brief return the number of instances of this object
  static size_t howMany(void) {
    return tCounter<cElementStructureKirchhoff>::howMany();
  }

  //! @brief return vector that contains elements degrees of freedom
  std::vector<eKnownDofs> getDofs(void) const;

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
  void getStressesCells(Vec &fullSolution, cElementVector &stresses,
                        cElementVector &stressesSec2) {
    trace(
        "ERROR: cElementStructureKirchhoff::getStressesCells() not "
        "implemented, yet.");
  }
  //! @brief compute stresses at the nodes
  //! @param FullSolution the full solution vector, each process holds a copy
  //! @param Stresses a vector that holds the global stresses
  void computeStressesNodes(Vec &fullSolution, Vec &stresses) {
    trace(
        "ERROR: cElementStructureKirchhoff::computeStressesNodes() not "
        "implemented, yet.");
    ExitApp();
  }

  //! @brief Compute the boundaryvelocity of the shell.
  //! @param n_vis number of values in vis
  //! @param vis vector of length 3 - absolut velocity in x(0),y(1) and z(2) -
  //! direction
  //! @param FullSolution Solution vector
  void computeVi(PetscReal *n_vis, std::vector<PetscReal> *vis,
                 const Vec &Solution);

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
                                const cElementStructureKirchhoff &other) {
  return other.write(os);
}

#endif
