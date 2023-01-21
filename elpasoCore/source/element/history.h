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

#ifndef INFAM_HISTORY_H
#define INFAM_HISTORY_H

#include <iostream>

#include "../misc/log/logging.h"
#include "../misc/mytypes.h"

//! @brief class, to store history of parameters
//! @author Marco Schauer
//! @date 05.06.2012
class cHistory : public virtual cLogging {
 private:
  PetscInt m_id;
  PetscReal m_X;
  PetscReal m_epsV;
  cElementVector *m_stress;
  cElementVector *m_strain;
  cElementVector *m_stressOld;
  cElementVector *m_strainOld;
  cElementVector *m_stressElastic;
  cElementVector *m_strainElastic;
  cElementVector *m_strainPlastic;
  cElementVector *m_backStress;

  bool m_yielding;
  bool m_idealPlastic;
  bool m_elasticValues;

 public:
  cHistory();
  cHistory(const cHistory &other);
  ~cHistory();

  //! set id
  void setId(PetscInt id) { m_id = id; }

  //! return history id
  PetscInt getId() { return m_id; }

  //! set Cap-Position
  void setX(PetscReal X) { m_X = X; }

  //! return Cap-Position
  PetscReal getX() { return m_X; }

  //! set volumetric strain
  void setEpsV(PetscReal depsV) { m_epsV += depsV; }

  //! return volumetric strain
  PetscReal getEpsV() { return m_epsV; }

  //! set stresses
  void setStress(cElementVector &stress);

  //! set strains
  void setStrain(cElementVector &strain);

  //! save stress at beginning of yielding
  void storeElasticStress();

  //! get absolute stress is element has yielded
  void getAbsoluteStress(cElementVector &stress);

  //! set back stress for kinematic hardening
  void setBackStress(cElementVector &stress);

  //! return stresses
  cElementVector *getStresses() { return m_stress; }
  //! return strain
  cElementVector *getStrain() { return m_strain; }
  //! return old stresses
  cElementVector *getStressOld() { return m_stressOld; }
  //! return old strain
  cElementVector *getStrainOld() { return m_strainOld; }
  //! return backStress
  cElementVector *getBackStress() { return m_backStress; }

  //! set yielding
  void setYielding(bool yielding) { m_yielding = yielding; }

  //! get yielding
  bool getYielding() { return m_yielding; }

  //! set elastic stress computation
  void setIdealPlastic(bool idealPlastic) { m_idealPlastic = idealPlastic; }

  //! get elastic stress computation
  bool getIdealPlastic() { return m_idealPlastic; }

  //! set true if elastic values exist
  void setElasticValues(bool arg) { m_elasticValues = arg; }

  //! check if element has elastic values
  bool hasElasticValues() { return m_elasticValues; }

  //! check compuation for convergence on element level
  bool checkConvergence();
};

#endif
