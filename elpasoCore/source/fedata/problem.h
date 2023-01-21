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

#ifndef INFAM_PROBLEM_H
#define INFAM_PROBLEM_H

#include "../analysis/analysis.h"
#include "../misc/log/logging.h"
#include "mesh.h"

//! @brief description of the problem
//! @author Dirk Clasen
//! @date 26.08.2004
//!
//! This class stores a pointer to the mesh and a pointer to the
//! analyis object. It is only needed because the xmlio parser only
//! accepts one paramter in function calls - but sometimes the mesh
//! as well as the analysis object is needed.
class cProblem : public virtual cLogging {
 private:
  cMesh *m_Mesh;           ///< pointer to FE-Mesh
  cAnalysis *m_Analysis;   ///< pointer to analysis parameters
  std::string m_Filename;  ///< name of the inputfile

 public:
  cProblem();
  cProblem(const cProblem &other);
  ~cProblem();

  //! return a pointer to the mesh
  inline cMesh *getMesh(void) const { return m_Mesh; }

  //! return the pointer to the analysis object
  inline cAnalysis *getAnalysis(void) const { return m_Analysis; }

  //! assign a new analysis object
  inline void setAnalysis(cAnalysis *ptr) { m_Analysis = ptr; }

  //! return the name of the inputfile
  inline std::string getFilename(void) const { return m_Filename; }

  //! save the name of the inputfile
  inline void setFilename(const std::string &Filename) {
    m_Filename = Filename;
    m_Mesh->setFilename(Filename);
  }
};

#endif
