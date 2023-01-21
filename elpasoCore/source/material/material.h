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

#ifndef INFAM_MATERIAL_H
#define INFAM_MATERIAL_H

#include "../element/history.h"
#include "../misc/counter.h"
#include "../misc/id.h"
#include "../misc/log/logging.h"
#include "../misc/mytypes.h"
#include "../misc/point.h"

//! @brief virtual base class for all materials
//! @author Dirk Clasen
//! @date 22.06.2005
class cMaterial : public cId,
                  private tCounter<cMaterial>,
                  public virtual cLogging {
 private:
  static PetscReal m_Omega;  ///< current angular frequency (used for
                             ///< frequencydepent materials)

 protected:
  PetscReal m_Rho;           ///< material's density
  std::string m_Identifier;  ///< identifier for material
  cHistory *m_history;  ///< history informations can be stored here, if needed.
                        ///< m_history is NULL per default!
  bool m_nl;  ///< default: false becomes true if element is linear element

 public:
  cMaterial();
  cMaterial(const cMaterial &other);
  virtual ~cMaterial();

  //! return the number of instances of this object
  static size_t howMany(void) { return tCounter<cMaterial>::howMany(); }

  //! returns material's density
  inline double getRho(void) const { return m_Rho; }

  //! assigns a new density to this material
  //! @param Rho the new density
  inline void setRho(double Rho) { m_Rho = Rho; }

  //! assigns a new angular frequency \f$ \omega \f$
  static void setOmega(const PetscReal &Omega);

  //! read the current angular frequency \f$ \omega \f$
  inline PetscReal getOmega(void) const { return m_Omega; }

  //! update the material paramters that are frequency dependent
  virtual void updateMaterial(void) = 0;

  //! set stresses
  void setHistory(cHistory *history) { m_history = history; }

  //! return material's identifier
  std::string getIdentifier(void) const { return m_Identifier; }

  //! assign a new identifier to the material
  void setIdentifier(const std::string &identifier) {
    m_Identifier = identifier;
  }

  //! return rho_x, rho_y and rho_z values for a given point with
  //! x,y,z-coordinates
  virtual void computeRhoValues(cPoint &point, cMatrix &rho_comp) = 0;

  //! return rho_c for a given point with x,y,z-coordinates
  //! and returns the bulkmodulus
  virtual PetscReal getRhoLambda(cPoint &point, PetscReal &rho_c) = 0;

  //! reads a single object of a stream
  //! @param is inputstream
  //! @return modified outputstream
  virtual std::istream &read(std::istream &is) = 0;

  //! write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &write(std::ostream &os) const = 0;

  //! write this object in XML
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &writeXml(std::ostream &os) const = 0;

  //! @see testmaterial.h
  bool isNonlinearElement() const { return m_nl; }
};

//! overloaded inputoperator
inline std::istream &operator>>(std::istream &is, cMaterial &other) {
  return other.read(is);
}

//! overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os, const cMaterial &other) {
  return other.write(os);
}

#endif
