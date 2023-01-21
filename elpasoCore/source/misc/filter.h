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

#ifndef INFAM_FRINKFILTER_H
#define INFAM_FRINKFILTER_H

#include <sstream>

#include "../fedata/mesh.h"
#include "mytypes.h"

//! @brief filter results of single node and write them to a text file
//! @author Dirk Clasen
//! @date 04.09.2007
class cFilter : public virtual cLogging {
 private:
  bool m_WantToFilter;  //! tells us if filtering is requested
  PetscInt m_NumNodes;  //! number of nodes to filter
  FILE **m_PltFiles;    //! pointer to the outputfiles

  std::vector<cNode *>
      m_PtrNodes;  //! pointer to nodes whose results will be filtered

 protected:
  bool filter_init;

 public:
  cFilter();
  ~cFilter();

  bool checkForFilter(void) const { return m_WantToFilter; }

  void initializeFilter(cMesh &myMesh);

  void filterFromSolution(PetscReal step, Vec &sol);
};

#endif
