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

#include "outputfile.h"

#include "../../fedata/entities.h"


cOutputfile::cOutputfile() {
  this->m_WantOutput = 0;
  modify_FilenameRoot = false;
  FilenameRoot_add.clear();

  this->m_InitPostProcess = false;
  this->m_numberofDofs = 0;
  this->m_dispX = 0;
  this->m_dispW = 0;
  this->m_DispVector.resize(cstNumberOfKnownDofs);
  for (int i = 0; i < cstNumberOfKnownDofs; ++i) this->m_DispVector[i] = false;

  m_OutSolution = NULL;
  m_OutVelocity = NULL;
  m_OutAcceleration = NULL;
  m_OutStressesCells = NULL;
  m_OutStressesCellsSec2 = NULL;
  m_OutStressesNodes = NULL;
}

cOutputfile::cOutputfile(cOutputfile &other) {
  this->m_WantOutput = other.m_WantOutput;
  modify_FilenameRoot = other.modify_FilenameRoot;
  FilenameRoot_add.assign(other.FilenameRoot_add);

  this->m_InitPostProcess = other.m_InitPostProcess;
  this->m_numberofDofs = other.m_numberofDofs;
  this->m_dispX = other.m_dispX;
  this->m_dispW = other.m_dispW;
  this->m_DispVector.resize(cstNumberOfKnownDofs);
  for (int i = 0; i < cstNumberOfKnownDofs; ++i)
    this->m_DispVector[i] = other.m_DispVector[i];
}

cOutputfile::~cOutputfile() {}

void cOutputfile::readResFile(const std::string &FileNameRoot, cMesh &MyMesh,
                              PetscInt &Rows, PetscInt &Steps) {
  char str[2000];

  std::string NameResFile = FileNameRoot + ".res";
  std::ifstream input;
  PetscInt NumNodes, IsComplex;

  // --- try to open the file which stores the assignment
  //     between rows and node's dof
  message("  reading '%s'\n", NameResFile.c_str());
  input.open(NameResFile.c_str());

  if (!input) {
    message("ERROR while opening file: '%s'\n", NameResFile.c_str());
    ExitApp();
  }

  // --- ignore 4 lines of comments
  input.getline(str, 2000);
  input.getline(str, 2000);
  input.getline(str, 2000);
  input.getline(str, 2000);

  input >> NumNodes >> IsComplex >> Rows >> Steps;

  // --- get end of line, afterwards ignore 3 lines of comment
  input.getline(str, 2000);
  input.getline(str, 2000);
  input.getline(str, 2000);
  input.getline(str, 2000);

  // --- assign rows to node's dofs
  message("  assign rows to node's degrees of freedom ... ");
  PetscInt SingleNodeId, SingleDof, SingleRow;

  for (PetscInt k = 0; k < Rows; k++) {
    input >> SingleNodeId >> SingleDof >> SingleRow;
    MyMesh.getNode(SingleNodeId)->setGlobalRow(SingleDof, SingleRow);
    MyMesh.getNode(SingleNodeId)->activateDof(cstAllDofs[SingleDof]);
  }

  input.close();

  PetscInt NodeCounter = 0;
  for (ItMapNodes it = MyMesh.getFirstNode(); it != MyMesh.getLastNode();
       it++) {
    // --- initialize sequential numbering of the nodes.
    //     This information will be used when the stresses are
    //     computed on element level.
    it->second->setGlobalSeqId(NodeCounter);
    NodeCounter++;
  }

  trace("    finished ! ");
}

void cOutputfile::readStpGzFile(const std::string &FileNameRoot, PetscInt Rows,
                                int step, PetscReal &StepVal, Vec &solution) {
  // --- create the filename
  std::stringstream Temp;
  Temp << ".";
  Temp.width(8);
  Temp.fill('0');
  Temp << step << ".stp.gz";

  // --- open zipped file that contains the solution vector
  std::string StpFileName = FileNameRoot + Temp.str();
  igzstream result(StpFileName.c_str());

  message("  resultsfile to read : %s\n", StpFileName.c_str());

  // --- timestep / frequency value
  // PetscReal StepVal;
  result >> StepVal;

  // --- read the nodal values
  PetscScalar CurrentValue;

#ifdef PETSC_USE_COMPLEX
  PetscReal RealPart, ImagPart;
#endif

  for (PetscInt k = 0; k < Rows; k++) {
#ifdef PETSC_USE_COMPLEX
    result >> RealPart >> ImagPart;
    CurrentValue = std::complex<PetscReal>(RealPart, ImagPart);
#else
    result >> CurrentValue;
#endif
    VecSetValues(solution, 1, &k, &CurrentValue, ADD_VALUES);
  }

  result.close();
}

void cOutputfile::initPostProcess() {
  // structure elements ---
  if (tCounter<cElementStructure>::howMany() != 0) {
    /*if ( tCounter<cElementStructureDisc>::howMany() != 0 )
    {
      m_DispVector[0]=true; //disp_x1
      m_DispVector[1]=true; //disp_x2
    }*/
    if (tCounter<cElementStructureBeam>::howMany() !=
        0)  // ||
            // tCounter<cElementStructureCable>::howMany()     !=0 ||
            // tCounter<cElementStructureDiscDrill>::howMany() !=0  )
    {
      m_DispVector[0] = true;  // disp_x1
      m_DispVector[1] = true;  // disp_x2
      m_DispVector[5] = true;  // disp_w3
    }

    if (tCounter<cElementStructureBeam3D>::howMany() !=
        0)  //||
            // tCounter<cElementStructureBeam3DNL>::howMany()        !=0 ||
            // tCounter<cElementStructurePlaneShell>::howMany()      !=0 ||
            // tCounter<cElementStructurePlaneShellDrill>::howMany() !=0  )
    {
      for (int i = 0; i < 6; ++i)
        m_DispVector[i] = true;  // disp_x123, disp_w123
    }

    if (  // tCounter<cElementStructureCable>::howMany() !=0 ||
        tCounter<cElementStructureBrick>::howMany() != 0)  //||
    // tCounter<cElementStructureTetra>::howMany() !=0    )
    {
      for (int i = 0; i < 3; ++i) m_DispVector[i] = true;  // disp_x123
    }

    /*if( tCounter<cElementStructureMindlin>::howMany() !=0 )
    {
      m_DispVector[2]=true; //disp_x3
      m_DispVector[3]=true; //disp_w1
      m_DispVector[4]=true; //disp_w2
    }*/
    if (tCounter<cElementStructureKirchhoff>::howMany() != 0) {
      m_DispVector[2] = true;   // disp_x3
      m_DispVector[13] = true;  // disp_dwdx
      m_DispVector[14] = true;  // disp_dwdy
      m_DispVector[15] = true;  // disp_dwdxy
    }
    if (tCounter<cElementStructureSpring>::howMany() != 0) {
      for (int i = 0; i < 6; ++i)
        m_DispVector[i] = true;  // disp_x123, disp_w123
    }
    if (tCounter<cElementStructureSpringBCx>::howMany() != 0) {
      m_DispVector[0] = true;  // disp_x1
    }
    if (tCounter<cElementStructureSpringBCy>::howMany() != 0) {
      m_DispVector[1] = true;  // disp_x2
    }
    if (tCounter<cElementStructureSpringBCz>::howMany() != 0) {
      m_DispVector[2] = true;  // disp_x3
    }
    if (tCounter<cElementStructureSpringBCrx>::howMany() != 0) {
      m_DispVector[3] = true;  // disp_w1
    }
    if (tCounter<cElementStructureSpringBCry>::howMany() != 0) {
      m_DispVector[4] = true;  // disp_w2
    }
    if (tCounter<cElementStructureSpringBCrz>::howMany() != 0) {
      m_DispVector[5] = true;  // disp_w3
    }
    if (tCounter<cElementStructureMass>::howMany() != 0) {
      for (int i = 0; i < 3; ++i) m_DispVector[i] = true;  // disp_x123
    }
  }  //--- structure elements

  // fluid elements ---
  if (tCounter<cElementFluid>::howMany() != 0) {
    m_DispVector[16] = true;  // fluid
  }                           //--- fluid elements

  m_InitPostProcess = true;
  for (int i = 0; i < cstNumberOfKnownDofs; ++i)
    if (m_DispVector[i] == true) m_numberofDofs++;
  // for vector output
  for (int i = 0; i < 3; ++i)
    if (m_DispVector[i] == true) m_dispX++;
  for (int i = 3; i < 6; ++i)
    if (m_DispVector[i] == true) m_dispW++;
  if (m_dispX > 1) m_numberofDofs = m_numberofDofs - m_dispX + 1;
  if (m_dispW > 1) m_numberofDofs = m_numberofDofs - m_dispW + 1;
}

void cOutputfile::createFilenameRoot(std::string &FilenameRoot, cMesh *myMesh) {
  FilenameRoot = myMesh->getFilename();
  if (FilenameRoot.rfind(cstInputFileSuffix) ==
      (FilenameRoot.length() - cstInputFileSuffix.length()))
    FilenameRoot.erase(FilenameRoot.length() - cstInputFileSuffix.length());

  if (modify_FilenameRoot) FilenameRoot.append(FilenameRoot_add);
}
