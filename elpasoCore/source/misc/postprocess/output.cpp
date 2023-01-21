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

#include "output.h"

#include "./../../basics/mpi/mpitools.h"
#include "./../../element/fluid/elementfluid.h"
#include "./../../element/structure/elementstructure.h"
#include "outputlev.h"
#include "outputstp.h"


cOutput::cOutput() {
  m_elementsAtNodesCounted = false;

  m_CreateAK3 = false;
  m_LogEnergy = false;
  m_LogVelocity = false;
}

cOutput::~cOutput() {
  std::list<cOutputfile *>::iterator it;
  for (it = m_OutToFile.begin(); it != m_OutToFile.end(); it++) {
    cOutputfile *ptrOutputfile = *it;
    delete ptrOutputfile;
  }
}

void cOutput::setOutSolution(Vec *vec) {
  std::list<cOutputfile *>::iterator it;
  for (it = m_OutToFile.begin(); it != m_OutToFile.end(); it++) {
    cOutputfile *ptrOutputfile = *it;
    ptrOutputfile->setOutSolution(vec);
  }
}

void cOutput::setOutVelocity(Vec *vec) {
  std::list<cOutputfile *>::iterator it;
  for (it = m_OutToFile.begin(); it != m_OutToFile.end(); it++) {
    cOutputfile *ptrOutputfile = *it;
    ptrOutputfile->setOutVelocity(vec);
  }
}

void cOutput::setOutAcceleration(Vec *vec) {
  std::list<cOutputfile *>::iterator it;
  for (it = m_OutToFile.begin(); it != m_OutToFile.end(); it++) {
    cOutputfile *ptrOutputfile = *it;
    ptrOutputfile->setOutAcceleration(vec);
  }
}

void cOutput::setOutStressesCells(Vec *vec) {
  std::list<cOutputfile *>::iterator it;
  for (it = m_OutToFile.begin(); it != m_OutToFile.end(); it++) {
    cOutputfile *ptrOutputfile = *it;
    ptrOutputfile->setOutStressesCells(vec);
  }
}

void cOutput::setOutStressesCellsSec2(Vec *vec) {
  std::list<cOutputfile *>::iterator it;
  for (it = m_OutToFile.begin(); it != m_OutToFile.end(); it++) {
    cOutputfile *ptrOutputfile = *it;
    ptrOutputfile->setOutStressesCellsSec2(vec);
  }
}

void cOutput::setOutStressesNodes(Vec *vec) {
  std::list<cOutputfile *>::iterator it;
  for (it = m_OutToFile.begin(); it != m_OutToFile.end(); it++) {
    cOutputfile *ptrOutputfile = *it;
    ptrOutputfile->setOutStressesNodes(vec);
  }
}

void cOutput::createFilenameRoot(std::string &FilenameRoot, cMesh *myMesh) {
  FilenameRoot = myMesh->getFilename();
  if (FilenameRoot.rfind(cstInputFileSuffix) ==
      (FilenameRoot.length() - cstInputFileSuffix.length()))
    FilenameRoot.erase(FilenameRoot.length() - cstInputFileSuffix.length());
}

//! inserst a new output-to-file-object
void cOutput::insertOutput(cOutputfile *m_out) { m_OutToFile.push_back(m_out); }

/*
//! get an pointer to an output-to-file-object
cOutputfile* cOutput::getOutput(PetscInt id) {

}
*/

void cOutput::writeResults(cMesh &MyMesh, const int &NrStep,
                           const double &Step) {
  if (!m_elementsAtNodesCounted) {
    countElementsAtNodes(MyMesh);
    m_elementsAtNodesCounted = true;
  }

  std::list<cOutputfile *>::iterator it;
  for (it = m_OutToFile.begin(); it != m_OutToFile.end(); it++) {
    cOutputfile *ptrOutputfile = *it;
    ptrOutputfile->writeResultsToFile(MyMesh, NrStep, Step);
  }
}

void cOutput::writeResults(std::list<PetscReal> &RealParts,
                           std::list<PetscReal> &ImagParts, cMesh &mesh) {
  // --- check if the user really wants to have this output
  if (wantResFile() != true) return;

  ogzstream raus;
  std::stringstream text;
  std::string CurrentName;

  // --- compose current filename
  CurrentName = mesh.getFilename();
  if (CurrentName.rfind(cstInputFileSuffix) ==
      (CurrentName.length() - cstInputFileSuffix.length()))
    CurrentName.erase(CurrentName.length() - cstInputFileSuffix.length());

  text << ".";
  text.width(8);
  text.fill('0');
  text << 1 << ".stp.gz";
  CurrentName += text.str();
  text.str("");

  message("  writing %s ...\n", CurrentName.c_str());

  // --- open file
  raus.open(CurrentName.c_str());

  // --- writing header of outputfile
  raus << "#  No.      Real           Imag" << std::endl;
  raus << "#-----------------------------------------------------------"
       << std::endl;

  PetscInt Count = 1;
  std::list<PetscReal>::const_iterator ItReal;
  std::list<PetscReal>::const_iterator ItImag = ImagParts.begin();

  for (ItReal = RealParts.begin(); ItReal != RealParts.end();
       ItReal++, ItImag++, Count++) {
    raus << std::setw(5) << Count;
    raus.setf(std::ios::scientific);
    raus << std::setw(15) << *ItReal;
    raus << std::setw(15) << *ItImag;
    raus << std::endl;
    raus.unsetf(std::ios::scientific);
  }

  // --- we're finished, close file
  raus.close();
}

// output function to output 4 columns, made for
void cOutput::writeResults(std::list<PetscReal> &RealParts,
                           std::list<PetscReal> &ImagParts,
                           std::list<PetscReal> &RealParts2,
                           std::list<PetscReal> &ImagParts2, cMesh &mesh) {
  // --- check if the user really wants to have this output
  if (wantResFile() != true) return;

  ogzstream raus;
  std::stringstream text;
  std::string CurrentName;

  // --- compose current filename
  CurrentName = mesh.getFilename();
  if (CurrentName.rfind(cstInputFileSuffix) ==
      (CurrentName.length() - cstInputFileSuffix.length()))
    CurrentName.erase(CurrentName.length() - cstInputFileSuffix.length());

  text << ".";
  text.width(8);
  text.fill('0');
  text << 1 << ".stp.gz";
  CurrentName += text.str();
  text.str("");

  message("  writing %s ...\n", CurrentName.c_str());

  // --- open file
  raus.open(CurrentName.c_str());

  // --- writing header of outputfile
  raus
      << "#  No.      Real           Imag       Real2(Vecnorm)   Imag2(Vecnorm)"
      << std::endl;
  raus << "#-----------------------------------------------------------"
       << std::endl;

  PetscInt Count = 1;
  std::list<PetscReal>::const_iterator ItReal;
  std::list<PetscReal>::const_iterator ItImag = ImagParts.begin();
  std::list<PetscReal>::const_iterator ItReal2 = RealParts2.begin();
  std::list<PetscReal>::const_iterator ItImag2 = ImagParts2.begin();

  for (ItReal = RealParts.begin(); ItReal != RealParts.end();
       ItReal++, ItImag++, ItReal2++, ItImag2++, Count++) {
    raus << std::setw(5) << Count;
    raus.setf(std::ios::scientific);
    raus << std::setw(15) << *ItReal;
    raus << std::setw(15) << *ItImag;
    raus << std::setw(15) << *ItReal2;
    raus << std::setw(15) << *ItImag2;
    raus << std::endl;
    raus.unsetf(std::ios::scientific);
  }

  // --- we're finished, close file
  raus.close();
}

//! check, if we have to create an info-file that contains the
//! mesh information
bool cOutput::wantResFile() {
  bool wantFile = false;
  std::list<cOutputfile *>::iterator it;
  for (it = m_OutToFile.begin(); it != m_OutToFile.end(); it++) {
    cOutputLEV *ptrOutputLEV = dynamic_cast<cOutputLEV *>(*it);
    if (ptrOutputLEV != NULL) {
      if (ptrOutputLEV->wantOutput() > 0) {
        wantFile = true;
      }
    }
    cOutputSTP *ptrOutputSTP = dynamic_cast<cOutputSTP *>(*it);
    if (ptrOutputSTP != NULL) {
      if (ptrOutputSTP->wantOutput() > 0) {
        wantFile = true;
      }
    }
  }

  return wantFile;
}

void cOutput::writeInfoFile(cMesh &mesh, const int &steps,
                            const int &NumberOfDofs) {
  // --- check if the user really wants to have this output
  if (wantResFile() != true) return;

  trace("  creating info-file ... ");

  PetscInt Counter = 0;
  PetscInt MaxNumber = cstNumberOfKnownDofs * mesh.getNumberOfNodes();
  PetscInt *GlobalRowsLocal = new PetscInt[MaxNumber];
  PetscInt *GlobalRowsGlobal = new PetscInt[MaxNumber];  // only used on rank 0

  for (PetscInt k = 0; k < MaxNumber; k++) {
    GlobalRowsLocal[k] = 0;
    GlobalRowsGlobal[k] = 0;
  }

  trace("   processing local nodes ... ");
  for (ItMapNodes it = mesh.getFirstNode(); it != mesh.getLastNode(); it++) {
    for (int dof = 0; dof < cstNumberOfKnownDofs; dof++) {
      if (it->second->checkIfActive(dof) == true) {
        // +1 , because MPI_Reduce will not work for row 0
        GlobalRowsLocal[Counter * cstNumberOfKnownDofs + dof] =
            it->second->getGlobalRow(dof) + 1;
      }
    }
    Counter++;
  }

  message("   adjust numbering of rows globally ... ");
  MPI_Reduce(&(GlobalRowsLocal[0]), &(GlobalRowsGlobal[0]), MaxNumber, MPIU_INT,
             MPI_BOR, 0, PETSC_COMM_WORLD);
  trace(" ok! ");

  // --- rank 0 will do the whole work
  if (PetscGlobalRank != 0) {
    delete[] GlobalRowsLocal;
    delete[] GlobalRowsGlobal;
    return;
  }

  trace("   starting to write on rank 0 ... ");

  // --- composing filename
  std::string CurrentName(mesh.getFilename());
  if (CurrentName.rfind(cstInputFileSuffix) ==
      (CurrentName.length() - cstInputFileSuffix.length()))
    CurrentName.erase(CurrentName.length() - cstInputFileSuffix.length());
  CurrentName += ".res";

  // --- opening file
  std::fstream raus;
  raus.open(CurrentName.c_str(), std::ios::out);

  raus
      << "#--------------------------------------------------------------------"
      << std::endl;
  raus << "# results" << std::endl;
  raus
      << "#--------------------------------------------------------------------"
      << std::endl;
  raus << "#    nnod   complex     rows   steps" << std::endl;
  raus << std::setw(9) << mesh.getNumberOfNodes();
#ifdef PETSC_USE_COMPLEX
  raus << std::setw(4) << "1";
#else
  raus << std::setw(4) << "0";
#endif
  raus << std::setw(15) << NumberOfDofs;
  raus << std::setw(9) << steps << std::endl;

  // --- no we'll write what's located in the rows of the solution vector
  raus
      << "#--------------------------------------------------------------------"
      << std::endl;
  raus << "#  node    dof  row" << std::endl;
  raus
      << "#--------------------------------------------------------------------"
      << std::endl;
  Counter = 0;
  for (ItMapNodes it = mesh.getFirstNode(); it != mesh.getLastNode(); it++) {
    for (int dof = 0; dof < cstNumberOfKnownDofs; dof++) {
      if (GlobalRowsGlobal[Counter * cstNumberOfKnownDofs + dof] > 0) {
        raus << it->first << " " << dof << " ";
        raus << GlobalRowsGlobal[Counter * cstNumberOfKnownDofs + dof] - 1
             << std::endl;
      }
    }
    Counter++;
  }

  delete[] GlobalRowsLocal;
  delete[] GlobalRowsGlobal;

  // --- close file
  raus.close();

  trace("  info-file written successfully");
}

/**
 * Writes a new inputfile
 */
void cOutput::createXmlOutputOfCurrentMesh(cMesh *myMesh,
                                           std::ostream &myAnalysis,
                                           std::ostream &myOutput,
                                           PetscInt step) {
  std::ofstream ak3Output;
  std::string filenameRoot;
  std::string filename;

  std::stringstream ssTemp;
  std::stringstream headOutput;
  std::stringstream nodeOutput;
  std::stringstream matOutput;
  std::stringstream groupOutput;
  std::stringstream bcOutput;
  std::stringstream nodalLoadsOutput;
  std::stringstream loadedNodesOutput;
  std::stringstream elementLoadsOutput;
  std::stringstream nodalPressuresOutput;
  std::stringstream footOutput;

  // --- join filename of outputfile
  cOutput::createFilenameRoot(filenameRoot, myMesh);
  ssTemp << "." << step << ".ak3";
  filename = filenameRoot + ssTemp.str();

  message("writing '%s'\n", filename.c_str());

  // -------------------------------------------------------------------------
  //   First, rank 0 is starting to write the file
  // -------------------------------------------------------------------------
  if (PetscGlobalRank == 0) {
    // open file
    ak3Output.open(filename.c_str());
    // write head information of file
    headOutput << "<FEM>" << std::endl;
    headOutput << "<Revision>6</Revision>" << std::endl;
    headOutput << "<File>" << filename << "</File>" << std::endl;
    headOutput << "<Desc>" << myMesh->getDescription() << "</Desc>"
               << std::endl;

    std::stringstream buf;
    buf << myAnalysis.rdbuf();
    buf << myOutput.rdbuf();
    headOutput << buf.str();

    // -----------------------------------------------------------------------
    //   output of nodes
    // -----------------------------------------------------------------------
    nodeOutput << "<Nodes N=\"" << myMesh->getNumberOfNodes() << "\">"
               << std::endl;
    for (ItMapNodes it = myMesh->getFirstNode(); it != myMesh->getLastNode();
         it++) {
      it->second->writeXml(nodeOutput);
      nodeOutput << std::endl;
    }
    nodeOutput << "</Nodes>" << std::endl;

    // -----------------------------------------------------------------------
    //   output of materials
    // -----------------------------------------------------------------------
    matOutput << "<Materials N=\"" << myMesh->getNumberOfMaterials() << "\">"
              << std::endl;
    for (ItMapMaterials it = myMesh->getFirstMaterial();
         it != myMesh->getLastMaterial(); it++) {
      it->second->writeXml(matOutput);
      matOutput << std::endl;
    }
    matOutput << "</Materials>" << std::endl;
  }
  // -----------------------------------------------------------------------
  //   output of the elements
  //   they are written in groups defined in MSC Patran
  // -----------------------------------------------------------------------
  ItMapGroups itGroups;
  for (itGroups = myMesh->getFirstGroup(); itGroups != myMesh->getLastGroup();
       itGroups++) {
    groupOutput << "<Elements";
    groupOutput << " N=\"" << itGroups->second->getNumberOfElementsInGroup()
                << "\"";
    groupOutput << " GroupId=\"" << itGroups->second->getId() << "\"";
    groupOutput << " Group=\"" << itGroups->second->getIdentifier() << "\"";
    groupOutput << " Material=\"" << itGroups->second->getMaterialId() << "\">"
                << std::endl;

    std::set<PetscInt>::const_iterator itElSet;
    for (itElSet = itGroups->second->getFirstElementId();
         itElSet != itGroups->second->getLastElementId(); itElSet++) {
      myMesh->getElement(*itElSet)->writeXml(groupOutput);
      groupOutput << std::endl;
    }
    groupOutput << "</Elements>" << std::endl;
  }
  // collect groups from all ranks
  for (int proc = 1; proc < PetscGlobalSize; proc++) {
    // PetscInt RemoteNumber = 0;
    if (PetscGlobalRank == proc) {
      cMpiTools::sendStringValue(groupOutput.str(), 0, 12300 + proc,
                                 PETSC_COMM_WORLD);
    } else if (PetscGlobalRank == 0) {
      std::string recvbuf;
      recvbuf =
          cMpiTools::receiveStringValue(proc, 12300 + proc, PETSC_COMM_WORLD);
      groupOutput << recvbuf;
    }
  }

  if (PetscGlobalRank == 0) {
    // -----------------------------------------------------------------------
    //   output of boundaryconditions
    // -----------------------------------------------------------------------
    bcOutput << "<NodeBCs N=\"" << myMesh->getNumberOfBoundaryConditions()
             << "\">" << std::endl;
    for (ItMapBoundaryConditions it = myMesh->getFirstBoundaryCondition();
         it != myMesh->getLastBoundaryCondition(); it++)
      it->second->writeXml(bcOutput);
    bcOutput << "</NodeBCs>" << std::endl;

    // -----------------------------------------------------------------------
    //   output of nodes to whom a boundary condition is applied
    // -----------------------------------------------------------------------
    bcOutput << "<FixedNodes>" << std::endl;
    for (ItMapNodes it = myMesh->getFirstNode(); it != myMesh->getLastNode();
         it++) {
      for (ItNodeBCs itb = it->second->getFirstBC();
           itb != it->second->getLastBC(); itb++) {
        bcOutput << "<Fixed>";
        bcOutput << "<Node>" << it->second->getId() << "</Node>";
        bcOutput << "<BC>" << itb->second->getId() << "</BC>";
        bcOutput << "</Fixed>" << std::endl;
      }
    }
    bcOutput << "</FixedNodes>" << std::endl;

    // -----------------------------------------------------------------------
    //   output of nodal forces
    // -----------------------------------------------------------------------
    nodalLoadsOutput << "<NodalForces N=\"" << myMesh->getNumberOfNodalForces()
                     << "\">" << std::endl;
    for (ItMapNodalForces it = myMesh->getFirstNodalForce();
         it != myMesh->getLastNodalForce(); it++)
      it->second->writeXml(nodalLoadsOutput);
    nodalLoadsOutput << "</NodalForces>" << std::endl;

    // -----------------------------------------------------------------------
    //   output of nodal moments
    // -----------------------------------------------------------------------
    nodalLoadsOutput << "<NodalMoments N=\""
                     << myMesh->getNumberOfNodalMoments() << "\">" << std::endl;
    for (ItMapNodalMoments it = myMesh->getFirstNodalMoment();
         it != myMesh->getLastNodalMoment(); it++)
      it->second->writeXml(nodalLoadsOutput);
    nodalLoadsOutput << "</NodalMoments>" << std::endl;

    // -----------------------------------------------------------------------
    //   output of nodes loaded with force
    // -----------------------------------------------------------------------
    loadedNodesOutput << "<NodesWithLoads>" << std::endl;
    for (ItMapNodes it = myMesh->getFirstNode(); it != myMesh->getLastNode();
         it++) {
      if (it->second->getNodalForce() != 0) {
        loadedNodesOutput << "<Loaded>";
        loadedNodesOutput << "<Node>" << it->second->getId() << "</Node>";
        loadedNodesOutput << "<Load>" << it->second->getNodalForce()->getId()
                          << "</Load>";
        loadedNodesOutput << "</Loaded>" << std::endl;
      }
    }
    loadedNodesOutput << "</NodesWithLoads>" << std::endl;

    // -----------------------------------------------------------------------
    //   output of nodes loaded with moment
    // -----------------------------------------------------------------------
    loadedNodesOutput << "<NodesWithMoments>" << std::endl;
    for (ItMapNodes it = myMesh->getFirstNode(); it != myMesh->getLastNode();
         it++) {
      if (it->second->getNodalMoment() != 0) {
        loadedNodesOutput << "<Loaded>";
        loadedNodesOutput << "<Node>" << it->second->getId() << "</Node>";
        loadedNodesOutput << "<Moment>" << it->second->getNodalMoment()->getId()
                          << "</Moment>";
        loadedNodesOutput << "</Loaded>" << std::endl;
      }
    }
    loadedNodesOutput << "</NodesWithMoments>" << std::endl;

    // -----------------------------------------------------------------------
    //   output of elementloads
    // -----------------------------------------------------------------------
    elementLoadsOutput << "<ElemLoads N=\"" << myMesh->getNumberOfElementLoads()
                       << "\">" << std::endl;
    for (ItMapElementLoads it = myMesh->getFirstElementLoad();
         it != myMesh->getLastElementLoad(); it++) {
      it->second->writeXml(elementLoadsOutput);
      elementLoadsOutput << std::endl;
    }
    elementLoadsOutput << "</ElemLoads>" << std::endl;

    // -----------------------------------------------------------------------
    //   output of elements with load applied
    // -----------------------------------------------------------------------
    elementLoadsOutput << "<LoadedElems>" << std::endl;
    for (ItMapElements itE = myMesh->getFirstElement();
         itE != myMesh->getLastElement(); itE++) {
      for (ItLoadsOnElementMap it = itE->second->getFirstElementLoad();
           it != itE->second->getLastElementLoad(); it++) {
        elementLoadsOutput << "<LoadedElem>";
        elementLoadsOutput << "<Id>" << itE->second->getId()
                           << "</Id>";  // no. of element with load applied
        elementLoadsOutput << "<Load>" << it->second->getId()
                           << "</Load>";  // no. of elementload
        elementLoadsOutput
            << "<Face>" << it->first
            << "</Face>";  // face of element to whom load is applied
        elementLoadsOutput << "</LoadedElem>";
        elementLoadsOutput << std::endl;
      }
    }
    elementLoadsOutput << "</LoadedElems>" << std::endl;
  }

  // -------------------------------------------------------------------------
  //   now we're almost finished
  // -------------------------------------------------------------------------
  if (PetscGlobalRank == 0) {
    // raus << "</InterfaceElements>" << std::endl;
    footOutput << "</FEM>" << std::endl;

    ak3Output << headOutput.str();
    ak3Output << nodeOutput.str();
    ak3Output << matOutput.str();
    ak3Output << groupOutput.str();
    ak3Output << bcOutput.str();
    ak3Output << nodalLoadsOutput.str();
    ak3Output << loadedNodesOutput.str();
    ak3Output << elementLoadsOutput.str();
    ak3Output << nodalPressuresOutput.str();
    ak3Output << footOutput.str();

    ak3Output.close();
  }

  MPI_Barrier(MPI_COMM_WORLD);

  message("file '%s' successfully written.\n", filename.c_str());
}

//! write the Mesh to a file in the different formats.
//! @param MyMesh given mesh
void cOutput::writeMesh(cMesh &MyMesh) {
  std::list<cOutputfile *>::iterator it;
  for (it = m_OutToFile.begin(); it != m_OutToFile.end(); it++) {
    cOutputfile *ptrOutputfile = *it;
    ptrOutputfile->writeMeshToFile(MyMesh);
  }
}

//! write the Mesh to a file in the different formats separeted by groups.
//! @param MyMesh given mesh
void cOutput::writeMeshGroup(cMesh &MyMesh) {
  std::list<cOutputfile *>::iterator it;
  for (it = m_OutToFile.begin(); it != m_OutToFile.end(); it++) {
    cOutputfile *ptrOutputfile = *it;
    ptrOutputfile->writeMeshToFileGroup(MyMesh);
  }
}

void cOutput::writeSolutionfromStepFile(cMesh &MyMesh, PetscInt &step) {
  PetscReal StepVal = 0.;
  PetscInt Rows = 0;
  PetscInt Steps = 0;
  std::string FilenameRoot;

  // --- join filename of outputfile
  cOutput::createFilenameRoot(FilenameRoot, &MyMesh);

  // --- read in the .res file
  cOutputfile::readResFile(FilenameRoot, MyMesh, Rows, Steps);

  Vec solution;
  VecCreate(PETSC_COMM_WORLD, &solution);
  VecSetSizes(solution, PETSC_DECIDE, Rows);
  VecSetFromOptions(solution);

  // --- read in the .stp file
  cOutputfile::readStpGzFile(FilenameRoot, Rows, step, StepVal, solution);

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  if (!m_elementsAtNodesCounted) {
    countElementsAtNodes(MyMesh);
    m_elementsAtNodesCounted = true;
  }

  // --- output of the results in different formats
  std::list<cOutputfile *>::iterator it;
  for (it = m_OutToFile.begin(); it != m_OutToFile.end(); it++) {
    cOutputfile *ptrOutputfile = *it;
    cOutputSTP *ptrOutputSTP = dynamic_cast<cOutputSTP *>(*it);
    if (ptrOutputSTP == NULL) setOutSolution(&solution);
    ptrOutputfile->writeResultsToFile(MyMesh, step, StepVal);
  }

  VecDestroy(&solution);
}

void cOutput::writeAllSolutionsfromStepFiles(cMesh &MyMesh,
                                             PetscInt poststeps) {
  PetscReal StepVal = 0.;
  PetscInt Rows = 0;
  PetscInt Steps = 0;
  std::string FilenameRoot;

  // --- join filename of outputfile
  cOutput::createFilenameRoot(FilenameRoot, &MyMesh);

  // --- read in the .res file
  cOutputfile::readResFile(FilenameRoot, MyMesh, Rows, Steps);

  if (~poststeps < 0) Steps = poststeps;

  Vec solution;
  VecCreate(PETSC_COMM_WORLD, &solution);
  VecSetSizes(solution, PETSC_DECIDE, Rows);
  VecSetFromOptions(solution);

  if (!m_elementsAtNodesCounted) {
    countElementsAtNodes(MyMesh);
    m_elementsAtNodesCounted = true;
  }

  for (int k = 1; k < (Steps + 1); k++) {
    VecZeroEntries(solution);

    // --- read in the .stp file
    cOutputfile::readStpGzFile(FilenameRoot, Rows, k, StepVal, solution);

    VecAssemblyBegin(solution);
    VecAssemblyEnd(solution);

    // --- output of the results in different formats
    std::list<cOutputfile *>::iterator it;
    for (it = m_OutToFile.begin(); it != m_OutToFile.end(); it++) {
      cOutputfile *ptrOutputfile = *it;
      cOutputSTP *ptrOutputSTP = dynamic_cast<cOutputSTP *>(*it);
      if (ptrOutputSTP == NULL) setOutSolution(&solution);
      ptrOutputfile->writeResultsToFile(MyMesh, k, StepVal);
    }
  }
}

void cOutput::countElementsAtNodes(cMesh &MyMesh) {
  // --- loop over elements and count the stress entries
  for (ItMapElements it = MyMesh.getFirstElement();
       it != MyMesh.getLastElement(); it++) {
    cElementFEM *ptrElement = MyMesh.getElement(it->second->getId());
    int nnod = ptrElement->getNumberOfNodes();

    cElementStructure *ptrELStructure =
        dynamic_cast<cElementStructure *>(ptrElement);
    if (ptrELStructure != NULL) {
      for (int k = 0; k < nnod; k++) {
        cNode *nd = ptrELStructure->getNode(k);

        nd->addAdjacentElement(sigma11);
        nd->addAdjacentElement(sigma22);
        nd->addAdjacentElement(sigma33);
        nd->addAdjacentElement(sigma12);
        nd->addAdjacentElement(sigma13);
        nd->addAdjacentElement(sigma23);
      }
    }

    cElementFluid *ptrELFluid = dynamic_cast<cElementFluid *>(ptrElement);
    if (ptrELFluid != NULL) {
      for (int k = 0; k < nnod; k++) {
        cNode *nd = ptrELFluid->getNode(k);

        nd->addAdjacentElement(gradpx);
        nd->addAdjacentElement(gradpy);
        nd->addAdjacentElement(gradpz);
      }
    }
  }
}

std::ostream &cOutput::write(std::ostream &os) {
  os << "outputoptions" << std::endl;
  std::list<cOutputfile *>::iterator it;
  for (it = m_OutToFile.begin(); it != m_OutToFile.end(); it++) {
    cOutputfile *ptrOutputfile = *it;
    ptrOutputfile->write(os);
  }
  os << " want a ak3 file                 : " << createAK3() << std::endl;

  return os;
}

std::istream &cOutput::read(std::istream &is) { return is; }

std::ostream &cOutput::writeXml(std::ostream &os) {
  os << "<Output>" << std::endl;
  std::list<cOutputfile *>::iterator it;
  for (it = m_OutToFile.begin(); it != m_OutToFile.end(); it++) {
    cOutputfile *ptrOutputfile = *it;
    ptrOutputfile->writeXml(os);
  }
  os << "  <ak3>" << createAK3() << "</ak3>" << std::endl;
  os << "</Output>" << std::endl;

  return os;
}
