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

#include "filter.h"

cFilter::cFilter()
{
  m_NumNodes = 10;
  filter_init = false;

  m_PltFiles = new FILE*[m_NumNodes];
  m_PtrNodes.resize( m_NumNodes, NULL );
}


cFilter::~cFilter()
{
  if (m_NumNodes != 0) {
    for (int k=0; k<m_NumNodes; k++) {
      PetscFClose(PETSC_COMM_WORLD, m_PltFiles[k]);
      m_PtrNodes[k] = NULL;
    }
  }

  delete[] m_PltFiles;
}


void cFilter::initializeFilter(cMesh &myMesh)
{
  PetscBool  found = PETSC_FALSE;
  PetscInt   *m_IdsOfNodesToFilter = new PetscInt[m_NumNodes];

  PetscOptionsGetIntArray(PETSC_NULL,"", "-filter", m_IdsOfNodesToFilter , &m_NumNodes, &found);

  if (found == PETSC_TRUE) {
    std::list<PetscInt> setRows;
    std::stringstream   text;
    std::string         basename = myMesh.getFilename();

    m_WantToFilter = true;

    // --- truncate .ak3 suffix
    if (basename.rfind(cstInputFileSuffix) == (basename.length()-cstInputFileSuffix.length()))
      basename.erase(basename.length()-cstInputFileSuffix.length());

    message("  Filtering of nodal results requested...\n");

    for (int k=0; k<m_NumNodes; k++) {
      message("  creating file for node %d ...\n", m_IdsOfNodesToFilter[k]);
      text.str("");
      text << basename << "." << m_IdsOfNodesToFilter[k] << ".flt";
      PetscFOpen(PETSC_COMM_WORLD, text.str().c_str(), "w", &(m_PltFiles[k]));

      m_PtrNodes[k] = myMesh.getNode( m_IdsOfNodesToFilter[k] );

      // --- write file header
      PetscFPrintf(PETSC_COMM_WORLD, m_PltFiles[k], "#  node %d\n", m_PtrNodes[k]->getId() );
#ifdef PETSC_USE_COMPLEX
      PetscFPrintf(PETSC_COMM_WORLD, m_PltFiles[k], "#  complex numbers\n");
#else
      PetscFPrintf(PETSC_COMM_WORLD, m_PltFiles[k], "#  real numbers\n");
#endif
      PetscFPrintf(PETSC_COMM_WORLD, m_PltFiles[k], "#  degrees of freedom:\n");
      PetscFPrintf(PETSC_COMM_WORLD, m_PltFiles[k], "#  step     ");
      for (int dof=0; dof < cstNumberOfKnownDofs; dof++) {
        if (m_PtrNodes[k]->checkIfActive(dof) == true) {
          PetscFPrintf(PETSC_COMM_WORLD, m_PltFiles[k], "%s ", sKnownDofsIdentifiers[dof].c_str());
        }
      }
      PetscFPrintf(PETSC_COMM_WORLD, m_PltFiles[k], "\n");
    }
	filter_init = true;
  }
  else {
    message("  no filtering of nodal results requested.\n");
  }

  delete[] m_IdsOfNodesToFilter;
}


void cFilter::filterFromSolution(PetscReal step, Vec &sol)
{
  if (PetscGlobalRank == 0) {
    PetscInt    row;
    PetscScalar value;

    for (int k=0; k<m_NumNodes; k++) {
      PetscFPrintf(PETSC_COMM_SELF, m_PltFiles[k], "%14.5e ", step);

      for (int dof=0; dof<cstNumberOfKnownDofs; dof++) {
        row = m_PtrNodes[k]->getGlobalRow( dof );

        if (row != -1) {
          VecGetValues(sol, 1, &row, &value);
#ifdef PETSC_USE_COMPLEX
          PetscFPrintf(PETSC_COMM_SELF, m_PltFiles[k], "%14.5e %14.5e", value.real(), value.imag());
#else
          PetscFPrintf(PETSC_COMM_SELF, m_PltFiles[k], "%14.5e ", value);
#endif
        }
      }

      PetscFPrintf(PETSC_COMM_SELF, m_PltFiles[k], "\n");
    }
  }
}

