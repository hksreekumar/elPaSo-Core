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

#ifndef INFAM_ELEMENTINTERFACE_KIRCH_H
#define INFAM_ELEMENTINTERFACE_KIRCH_H

#include "elementinterfacestructural.h"

//! @brief coupling acoustic domains and Kirchhoff plate elements
//! @author Dirk Clasen
//! @date 16.08.2006
class cElementInterfaceKirchhoff
    : public cElementInterfaceStructural,
      private tCounter<cElementInterfaceKirchhoff> {
 private:
  cMaterialFluid *m_MaterialFluid;

 protected:
  //! @brief compute the coupling matrix
  void computeCouplingMatrix(cElementMatrix &C);

  //! @brief return the density of the fluid. This is needed by
  //! addContributionFrequencyDomain() and casting to fluid material
  //! does not work for coupling of Mindlin plate and 3d poroelastic
  //! elements
  PetscScalar getRhoF(void) const { return m_MaterialFluid->getRhoOmega(); }

  //! @brief returns the element type
  eElementType getElementType(void) const { return ShellAcoustic; }

 public:
  //! @brief Constructor
  cElementInterfaceKirchhoff(short NumberOfNodes);
  //! @brief Copy constructor
  cElementInterfaceKirchhoff(const cElementInterfaceKirchhoff &other);
  //! @brief Destructor
  ~cElementInterfaceKirchhoff();

  //! @brief assing a material of one of the coupled domains
  void setMaterial(cMaterial *ptrMaterial, const int &domain);

  //! @brief return one of the materials of the coupled domains
  cMaterial *getMaterial(const int &domain) const;

  //! @brief denotes the shape of this element
  eElementShape getElementShape(void) const;

  //! @brief returns the number of instances of this objecttype
  static size_t howMany(void) {
    return tCounter<cElementInterfaceKirchhoff>::howMany();
  }

  //! @brief return a vector that tells us which structural domains of
  //! freedom are coupled to the acoustic domain.
  std::vector<eKnownDofs> getDofsStructure(void) const;

  //! @brief for correct computation of the coupling term we need to
  //! find the lower left corner of the element and rotate
  //! the element to that node (only import for elements with 4 nodes)
  void findLowerCorner(void);

  //! @brief empty function
  inline void computeLoadVector(cElementVector &LV) {}

  //! @brief write this object to a stream
  //! @param os outputstream
  //! @return der modified outputstream
  std::ostream &write(std::ostream &os) const;

  //! @brief write this object in XML format to a stream
  //! @param os outputstream
  //! @return der modified outputstream
  std::ostream &writeXml(std::ostream &os) const;
};

//! @brief overloaded outputoperator
inline std::ostream &write(std::ostream &os,
                           const cElementInterfaceKirchhoff &other) {
  return other.write(os);
}

#endif
