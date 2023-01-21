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

#include "./analysis.h"

#include "../basics/mpi/mpitools.h"
#include "../element/interface/elementinterfacestructural.h"
#include "../math/mathlibrary.h"

cAnalysis::cAnalysis() {
  ierr = 0;

  m_MakeSymmetric = false;

  m_ComputeStressesCells = false;
  m_ComputeStressesNodes = false;

  m_nl = false;
  m_elementsDistributed = false;

  m_RowsBCDofs = NULL;
  m_BCValues = NULL;

  m_StressesCells = NULL;
  m_StressesCellsSec2 = NULL;
  m_StressesNodes = NULL;
}

cAnalysis::~cAnalysis() {
  // empty
}

void cAnalysis::checkForLinearAndOrNonlinearElementsAndMaterials(
    cMesh &MyMesh) {
  for (ItMapElements it = MyMesh.getFirstElement();
       it != MyMesh.getLastElement(); it++) {
    if (it->second->isNonlinearElement() == true) m_nl = true;
    if (it->second->getMaterial()->isNonlinearElement() == true) m_nl = true;
    if (m_nl == true) {
      it = MyMesh.getLastElement();
      it--;
    }
  }
  std::cout << "[" << PetscGlobalRank << "] lin: " << m_nl << "\n";
}

void cAnalysis::distributeMesh(cMesh &MyMesh) {
  // Set second id which runs from 0 to n, this one is used for stress
  // computation
  if (checkForStressCellsComputation() == true ||
      checkForStressNodesComputation() == true) {
    int i = 0;
    for (ItMapElements it = MyMesh.getFirstElement();
         it != MyMesh.getLastElement(); it++) {
      it->second->setId0n(i);
      i++;
    }
  }

  // --- leave, if only one process exits
  if (PetscGlobalSize == 1) return;

  trace("  distribute elements ... ");

  // -------------------------------------------------------------------------
  //   distribute elements
  // -------------------------------------------------------------------------
  PetscInt *elems = 0;
  PetscInt nelem = MyMesh.getNumberOfElements();

  elems = new PetscInt[PetscGlobalSize];

  PetscInt rest = nelem % PetscGlobalSize;
  PetscInt ganz = (nelem - rest) / PetscGlobalSize;

  for (int k = 0; k < PetscGlobalSize; k++) elems[k] = ganz;

  for (int k = rest; k > 0; k--) elems[PetscGlobalSize - k] += 1;

  for (int k = 0; k < PetscGlobalSize; k++)
    message("    rank = %d : elements = %d\n", k, elems[k]);

  // -- determine the number of elements to stay in local mesh. Purge denotes
  //    the number of elements at the beginning of the mesh which have to be
  //    deleted.
  PetscInt Keep = elems[PetscGlobalRank];
  PetscInt Purge = 0;
  for (int k = 0; k < PetscGlobalRank; k++) Purge += elems[k];

  MyMesh.purgeElements(Purge, Keep);

  // -------------------------------------------------------------------------
  //   distribute interface elements
  // -------------------------------------------------------------------------
  nelem = MyMesh.getNumberOfInterfaceElements();
  if (nelem > 0) {
    rest = nelem % PetscGlobalSize;
    ganz = (nelem - rest) / PetscGlobalSize;

    for (int k = 0; k < PetscGlobalSize; k++) elems[k] = ganz;

    for (int k = rest; k > 0; k--)  // Rest aufteilen
      elems[PetscGlobalSize - k] += 1;

    for (int k = 0; k < PetscGlobalSize; k++)
      message("    rank = %d : interface elements = %d\n", k, elems[k]);

    // -- determine the number of elements to stay in local mesh. Purge denotes
    //    the number of elements at the beginning of the mesh which have to be
    //    deleted.
    Keep = elems[PetscGlobalRank];
    Purge = 0;
    for (int k = 0; k < PetscGlobalRank; k++) Purge += elems[k];

    MyMesh.purgeInterfaceElements(Purge, Keep);
  }

  // -------------------------------------------------------------------------
  //   distribute FE BE coupling elements
  // -------------------------------------------------------------------------
  nelem = MyMesh.getNumberOfCplFemBemElements();
  if (nelem > 0) {
    rest = nelem % PetscGlobalSize;
    ganz = (nelem - rest) / PetscGlobalSize;

    for (int k = 0; k < PetscGlobalSize; k++) elems[k] = ganz;

    for (int k = rest; k > 0; k--)  // distribute remaining elements
      elems[PetscGlobalSize - k] += 1;

    for (int k = 0; k < PetscGlobalSize; k++)
      message("    rank = %d : FE BE coupling elements = %d\n", k, elems[k]);

    // -- determine the number of elements to stay in local mesh. Purge denotes
    //    the number of elements at the beginning of the mesh which have to be
    //    deleted.
    Keep = elems[PetscGlobalRank];
    Purge = 0;
    for (int k = 0; k < PetscGlobalRank; k++) Purge += elems[k];

    MyMesh.purgeCouplingElements(Purge, Keep);
  }

  // --- freeing memory
  delete[] elems;
  m_elementsDistributed = true;
}

void cAnalysis::applyNodalLoads(cMesh &MyMesh, const PetscReal &factor,
                                const PetscInt &step) {
  applyNodalForces(MyMesh, factor, step);
  applyNodalMoments(MyMesh, factor, step);
}

void cAnalysis::applyNodalForces(cMesh &MyMesh, const PetscReal &factor,
                                 const PetscInt &step) {
  if (!getMute()) trace("  inserting nodal forces ...");

  if (MyMesh.getNumberOfNodalForces() == 0) {
    trace("   there are no nodal forces to apply - leaving!");
    return;
  }

  // --- determine the rows located on local process
  PetscInt Istart, Iend;
  MatGetOwnershipRange(m_K, &Istart, &Iend);

  // --- loop across nodes and look for nodal forces
  cNodalForce *ptrForce = NULL;
  for (ItMapNodes it = MyMesh.getFirstNode(); it != MyMesh.getLastNode();
       it++) {
    // --- get the nodal force of current node, returns NULL if no
    //     load is applied to the current node
    ptrForce = it->second->getNodalForce();

    // --- check if we got a valid pointer
    if (ptrForce != NULL) {
      // ---------------------------------------------------------------------
      //   nodal force structure
      // ---------------------------------------------------------------------
      cNodalForceStructure *ptrNodalForceStructure =
          dynamic_cast<cNodalForceStructure *>(ptrForce);

      if (ptrNodalForceStructure != NULL) {
        // --- check x1 direction
        if (it->second->checkIfActive(disp_x1) == true) {
          PetscInt AffectedRow = it->second->getGlobalRow(disp_x1);

          // --- is local process responsible for the affected row?
          if ((AffectedRow >= Istart) && (AffectedRow < Iend)) {
            // --- if we're responsible read the value and insert
            //     it into the vector
            PetscScalar Value = (*ptrNodalForceStructure)[0];
            ierr = VecSetValue(m_F, AffectedRow, factor * Value, ADD_VALUES);
            INFAMCHKERRQ(ierr);
          }
        }

        // --- check x2 direction
        if (it->second->checkIfActive(disp_x2) == true) {
          PetscInt AffectedRow = it->second->getGlobalRow(disp_x2);

          // --- is local process responsible for the affected row?
          if ((AffectedRow >= Istart) && (AffectedRow < Iend)) {
            // --- if we're responsible read the value and insert
            //     it into the vector
            PetscScalar Value = (*ptrNodalForceStructure)[1];
            ierr = VecSetValue(m_F, AffectedRow, factor * Value, ADD_VALUES);
            INFAMCHKERRQ(ierr);
          }
        }

        // --- check x3 direction
        if (it->second->checkIfActive(disp_x3) == true) {
          PetscInt AffectedRow = it->second->getGlobalRow(disp_x3);

          // --- is local process responsible for the affected row?
          if ((AffectedRow >= Istart) && (AffectedRow < Iend)) {
            // --- if we're responsible read the value and insert
            //     it into the vector
            PetscScalar Value = (*ptrNodalForceStructure)[2];
            ierr = VecSetValue(m_F, AffectedRow, factor * Value, ADD_VALUES);
            INFAMCHKERRQ(ierr);
          }
        }
      }

      // ---------------------------------------------------------------------
      //   nodal force fluid
      // ---------------------------------------------------------------------
      cNodalForceFluid *ptrNodalForceFluid =
          dynamic_cast<cNodalForceFluid *>(ptrForce);

      if (ptrNodalForceFluid != NULL) {
        // --- check x1 direction
        if (it->second->checkIfActive(fluid) == true) {
          PetscInt AffectedRow = it->second->getGlobalRow(fluid);

          // --- is local process responsible for the affected row?
          if ((AffectedRow >= Istart) && (AffectedRow < Iend)) {
            // --- if we're responsible read the value and insert
            //     it into the vector
#ifdef PETSC_USE_COMPLEX
            PetscScalar Value =
                4. * M_PI * ptrNodalForceFluid->getSourceValue() * factor;
            ierr = VecSetValue(m_F, AffectedRow, factor * Value, ADD_VALUES);
            INFAMCHKERRQ(ierr);
#endif
          }
        }
      }
    }

    ptrForce = NULL;
  }

  ierr = VecAssemblyBegin(m_F);
  INFAMCHKERRQ(ierr);
  ierr = VecAssemblyEnd(m_F);
  INFAMCHKERRQ(ierr);

  if (!getMute()) trace("  nodal forces inserted ...");
}

void cAnalysis::applyNodalMoments(cMesh &MyMesh, const PetscReal &factor,
                                  const PetscInt &step) {
  if (!getMute()) trace("  inserting nodal moments ...");

  if (MyMesh.getNumberOfNodalMoments() == 0) {
    trace("   there are no nodal moments to apply - leaving!");
    return;
  }

  // --- determine the rows located on local process
  PetscInt Istart, Iend;
  MatGetOwnershipRange(m_K, &Istart, &Iend);

  // --- loop across nodes and look for nodal moments
  cNodalMoment *ptrMoment = NULL;
  for (ItMapNodes it = MyMesh.getFirstNode(); it != MyMesh.getLastNode();
       it++) {
    // --- get the nodal moment of current node, returns NULL if no
    //     load is applied to the current node
    ptrMoment = it->second->getNodalMoment();

    // --- check if we got a valid pointer
    if (ptrMoment != NULL) {
      // ---------------------------------------------------------------------
      //   nodal moment structure
      // ---------------------------------------------------------------------
      cNodalMomentStructure *ptrNodalMomentStructure =
          dynamic_cast<cNodalMomentStructure *>(ptrMoment);

      if (ptrNodalMomentStructure != NULL) {
        // --- check w1 rotation
        if (it->second->checkIfActive(disp_w1) == true) {
          PetscInt AffectedRow = it->second->getGlobalRow(disp_w1);

          // --- is local process responsible for the affected row?
          if ((AffectedRow >= Istart) && (AffectedRow < Iend)) {
            // --- if we're responsible read the value and insert
            //     it into the vector
            PetscScalar Value = (*ptrNodalMomentStructure)[0];
            ierr = VecSetValue(m_F, AffectedRow, factor * Value, ADD_VALUES);
            INFAMCHKERRQ(ierr);
          }
        }

        // --- check w2 rotation
        if (it->second->checkIfActive(disp_w2) == true) {
          PetscInt AffectedRow = it->second->getGlobalRow(disp_w2);

          // --- is local process responsible for the affected row?
          if ((AffectedRow >= Istart) && (AffectedRow < Iend)) {
            // --- if we're responsible read the value and insert
            //     it into the vector
            PetscScalar Value = (*ptrNodalMomentStructure)[1];
            ierr = VecSetValue(m_F, AffectedRow, factor * Value, ADD_VALUES);
            INFAMCHKERRQ(ierr);
          }
        }

        // --- check w3 rotation
        if (it->second->checkIfActive(disp_w3) == true) {
          PetscInt AffectedRow = it->second->getGlobalRow(disp_w3);

          // --- is local process responsible for the affected row?
          if ((AffectedRow >= Istart) && (AffectedRow < Iend)) {
            // --- if we're responsible read the value and insert
            //     it into the vector
            PetscScalar Value = (*ptrNodalMomentStructure)[2];
            ierr = VecSetValue(m_F, AffectedRow, factor * Value, ADD_VALUES);
            INFAMCHKERRQ(ierr);
          }
        }
      }
    }

    ptrMoment = NULL;
  }

  ierr = VecAssemblyBegin(m_F);
  INFAMCHKERRQ(ierr);
  ierr = VecAssemblyEnd(m_F);
  INFAMCHKERRQ(ierr);

  if (!getMute()) trace("  nodal moments inserted ...");
}

void cAnalysis::getGlobalBCs(cMesh &MyMesh, const PetscReal omega) {
  trace("getGlobalBCs ... ");
  // --- determine the rows located on local process
  PetscInt Istart, Iend;
  m_NumFixedRows = 0;  // (ll)

  MatGetOwnershipRange(m_K, &Istart, &Iend);

  // --- loop across nodes and look for boundary conditions
  for (ItMapNodes it = MyMesh.getFirstNode(); it != MyMesh.getLastNode(); it++)
    // --- loop over boundary conditions applied to current node
    for (ItNodeBCs itBCs = it->second->getFirstBC();
         itBCs != it->second->getLastBC(); itBCs++)
      // --- loop over all degrees of freedom
      for (int dof = 0; dof < cstNumberOfKnownDofs; dof++)
        // --- check if the dof exists at the node and the dof is fixed by the
        // boundary condition
        if (it->second->checkIfActive(dof) == true)
          if (itBCs->second->checkIfFixed(dof) == true)
            if (it->second->getGlobalRow(dof) < Iend &&
                it->second->getGlobalRow(dof) >= Istart)
              m_NumFixedRows++;

  if (m_RowsBCDofs != NULL) delete[] m_RowsBCDofs;
  if (m_BCValues != NULL) delete[] m_BCValues;
  m_RowsBCDofs = new PetscInt[m_NumFixedRows];
  m_BCValues = new PetscScalar[m_NumFixedRows];

  PetscInt BCNum = 0;

  // --- loop across nodes and look for boundary conditions
  for (ItMapNodes it = MyMesh.getFirstNode(); it != MyMesh.getLastNode(); it++)
    // --- loop over boundary conditions applied to current node
    for (ItNodeBCs itBCs = it->second->getFirstBC();
         itBCs != it->second->getLastBC(); itBCs++)
      // --- loop over all degrees of freedom
      for (int dof = 0; dof < cstNumberOfKnownDofs; dof++)
        // --- check if the dof exists at the node and the dof is fixed by the
        // boundary condition
        if (it->second->checkIfActive(dof) == true)
          if (itBCs->second->checkIfFixed(dof) == true)
            if (it->second->getGlobalRow(dof) < Iend &&
                it->second->getGlobalRow(dof) >= Istart) {
              m_RowsBCDofs[BCNum] = it->second->getGlobalRow(dof);
              m_BCValues[BCNum] = itBCs->second->getPrescribedValue(
                  cstAllDofs[dof], it->second, omega);
              BCNum++;
            }
}

// ---------------------------------------------------------------------------
//   insert boundary conditions into system matrices M, K and vector F
// ---------------------------------------------------------------------------
bool cAnalysis::insertBoundaryConditions(PetscBool setBCValue2zero,
                                         PetscReal factor,
                                         PetscBool setBCMassMat) {
  PetscScalar bc_value;

  if (!getMute()) trace("  Insert boundary conditions ... ");

  ierr = MatZeroRows(m_K, m_NumFixedRows, m_RowsBCDofs, 1.0, PETSC_NULL,
                     PETSC_NULL);
  INFAMCHKERRQ(ierr);  // ORIGINAL********
  if (setBCMassMat == PETSC_TRUE)
    ierr = MatZeroRows(m_M, m_NumFixedRows, m_RowsBCDofs, 0.0, PETSC_NULL,
                       PETSC_NULL);
  INFAMCHKERRQ(ierr);

  PetscInt Istart, Iend;
  PetscScalar *F;
  ierr = VecGetOwnershipRange(m_F, &Istart, &Iend);
  INFAMCHKERRQ(ierr);
  ierr = VecGetArray(m_F, &F);
  INFAMCHKERRQ(ierr);
  for (PetscInt i = 0; i < m_NumFixedRows; i++) {
    PetscInt idx = m_RowsBCDofs[i] - Istart;
    if (setBCValue2zero == PETSC_TRUE)
      F[idx] = 0.0;
    else
      F[idx] = factor * m_BCValues[i];
    // PetscPrintf(PETSC_COMM_SELF,"BC index %d value %f\n",
    // (int)m_RowsBCDofs[i],(double)m_BCValues[i]);
  }
  ierr = VecRestoreArray(m_F, &F);
  INFAMCHKERRQ(ierr);

  if (!getMute()) trace(" ok");
  return true;
}

void cAnalysis::countElementsPerNode(cMesh &myMesh) {
  if (!getMute()) trace("  count the number of elements per node ... ");
  // initialize with zero
  for (ItMapNodes it = myMesh.getFirstNode(); it != myMesh.getLastNode();
       it++) {
    it->second->setNumEle(0);
  }

  for (ItMapElements it = myMesh.getFirstElement();
       it != myMesh.getLastElement(); it++) {
    for (int k = 0; k < it->second->getNumberOfNodes(); k++) {
      it->second->getNode(k)->countElement();
    }
  }
}

/**
 * @todo Clean the behaviour of the spring elements
 */
void cAnalysis::activateDofsAtNodes(cMesh &myMesh) {
  countElementsPerNode(myMesh);  // testing
  // -- activate degrees of freedom at nodes
  if (!getMute())
    trace("  determine degrees of freedom existing at nodes ... ");
  for (ItMapElements it = myMesh.getFirstElement();
       it != myMesh.getLastElement(); it++) {
    std::vector<eKnownDofs> ElementsDofs(it->second->getDofs());

    for (int k = 0; k < it->second->getNumberOfNodes(); k++) {
      for (int d = 0; d < (int)ElementsDofs.size(); d++)
        it->second->getNode(k)->activateDof((short)ElementsDofs[d]);
    }
  }

  // -------------------------------------------------------------------------
  //  Numbering degrees of freedom of the nodes. This information will be
  //  used to insert the element contribution properly into the global
  //  matrices. Here all degrees of freedom are numbered not only the ones
  //  that exist for the local portion of elements. This is done because
  //  we have to insert also interface elements.
  // -------------------------------------------------------------------------*/
  if (!getMute()) trace("  numbering active degrees of freedom at nodes ...");

  m_NumberOfUnknowns = 0;
  PetscInt NodeCounter = 0;
  for (ItMapNodes it = myMesh.getFirstNode(); it != myMesh.getLastNode();
       it++) {
    for (short i = 0; i < cstNumberOfKnownDofs; i++) {
      int active = 0;
      int activeSum = 0;
      if (it->second->checkIfActive(i) == true) active = 1;
      MPI_Allreduce(&active, &activeSum, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
      if (activeSum > 0) {
        it->second->setGlobalRow(i, m_NumberOfUnknowns);
        m_NumberOfUnknowns++;
      }
    }
    // --- also initialize sequential numbering of the nodes.
    //     This information will be used when the stresses are
    //     computed on element level.
    it->second->setGlobalSeqId(NodeCounter);
    NodeCounter++;
  }
  if (!getMute()) trace(" numbering ok");
}

void cAnalysis::optimizeMatrixPattern(cMesh &myMesh) {
  initializePETScObjects(myMesh);

  trace("optimizing structure of matrix ...");

  // --- size of local matrix portion
  PetscInt FirstRow, LastRow, LocalSize;
  MatSetUp(m_K);

  MatGetOwnershipRange(m_K, &FirstRow, &LastRow);
  LocalSize = LastRow - FirstRow;

  PetscBool found = PETSC_FALSE;
  PetscInt reorder = 0;
  std::vector<PetscInt> Ele2del;

  PetscOptionsGetInt(PETSC_NULL, "", "-reorder", &reorder, &found);
  if (found == PETSC_TRUE) {
    std::vector<std::set<PetscInt> > Pattern;
    switch (reorder) {
      case 1:  // ParMETIS : reduce fill-in
        Pattern.resize(LocalSize);
        getLocalSparsityPattern(myMesh, FirstRow, LastRow, Pattern);
        renumberMatrixUsingParmetis(myMesh, FirstRow, LastRow, Pattern);
        break;

      case 2:  // reverse Cuthill McKee
        Pattern.resize(m_NumberOfUnknowns);
        if (PetscGlobalRank == 0) {
          getLocalSparsityPattern(myMesh, 0, m_NumberOfUnknowns, Pattern);
        }
        renumberMatrixUsingBoost(myMesh, FirstRow, LastRow, Pattern);
        break;

      case 3:  // Reorder using only mesh information
        reorderUsingParMetis(myMesh, 3, Ele2del);
        break;

      case 4:  // Reorder using mesh and co-ordinate information
        reorderUsingParMetis(myMesh, 4, Ele2del);
        break;

      case 5:  // Reorder using mesh and co-ordinate information
        reorderUsingParMetis(myMesh, 5, Ele2del);
        break;

      default:
        message(" *** INVALID VALUE FOR -order : %d\n", reorder);
        trace("     CHECK MANUAL FOR VALID VALUES");
    }
  } else
    trace(" *** no reordering of matrix requested by user");

  std::vector<std::set<PetscInt> > Pattern(LocalSize);
  getLocalSparsityPattern(myMesh, FirstRow, LastRow, Pattern);
  preallocateMatrix(FirstRow, LastRow, Pattern);
  MatGetOwnershipRange(m_K, &FirstRow, &LastRow);
  Pattern.clear();

  if (found == PETSC_TRUE && (reorder >= 3 && reorder <= 5)) {
    distributeParMetis(myMesh, Ele2del);
  } else
    distributeMesh(myMesh);
}

void cAnalysis::renumberMatrixUsingParmetis(
    cMesh &myMesh, const PetscInt &FirstRow, const PetscInt &LastRow,
    std::vector<std::set<PetscInt> > &Pattern) {
#ifndef PETSC_HAVE_PARMETIS
  trace("WARNING: matrix won't be renumbered using ParMETIS");
  trace("         because ParMETIS is not installed on this system");
  return;
#else

  trace("  compute fill-in reducing ordering of matrix using ParMETIS ...");

  // --- setup vtxdist array
  MPI_Status status;
  PetscInt LocalSize = LastRow - FirstRow;

  PetscInt *vtxdist = new PetscInt[PetscGlobalSize + 1];
  for (int k = 0; k < PetscGlobalSize + 1; k++) vtxdist[k] = 0;

  // --- all beside rank 0 receive array from their left neighbour
  //     but they have to wait until their neighbour itself received
  //     the array and added it local number of entries ...
  if (PetscGlobalRank != 0) {
    MPI_Recv(vtxdist, PetscGlobalSize + 1, MPIU_INT, PetscGlobalRank - 1, 1001,
             PETSC_COMM_WORLD, &status);
  }

  vtxdist[PetscGlobalRank + 1] = vtxdist[PetscGlobalRank] + LocalSize;

  // --- send my neighbour to the right the partial array
  if (PetscGlobalRank != PetscGlobalSize - 1)
    MPI_Send(vtxdist, PetscGlobalSize + 1, MPIU_INT, PetscGlobalRank + 1, 1001,
             PETSC_COMM_WORLD);

  // --- last proc sends complete array to all procs
  MPI_Bcast(vtxdist, PetscGlobalSize + 1, MPIU_INT, PetscGlobalSize - 1,
            PETSC_COMM_WORLD);

  PetscInt *xadj = new PetscInt[LocalSize + 1];  // index in adjacency array
  xadj[0] = 0;

  for (int k = 0; k < (int)Pattern.size(); k++) {
    xadj[k + 1] = xadj[k] + Pattern[k].size() -
                  1;  // don't count edge from vertex k to vertex k
  }

  // --- setup adjacency graph
  std::set<PetscInt>::iterator itP;
  PetscInt m = 0;
  PetscInt globNum = 0;
  PetscInt *adjc = new PetscInt[xadj[LocalSize]];
  PetscInt *original = new PetscInt[LocalSize];

  for (int k = 0; k < (int)Pattern.size(); k++) {
    globNum = vtxdist[PetscGlobalRank] + k;
    original[k] = globNum;
    for (itP = Pattern[k].begin(); itP != Pattern[k].end(); itP++) {
      if (*itP != globNum) {
        adjc[m] = *itP;
        m++;
      }
    }
  }

  // --- options[0] == 0 : use default settings, options[2..3] are ignored
  //                   1 : user defined settings, evaluate options[2..3]
  //     options[1] == 0 : no timing results
  //                   1 : show timing results
  //     options[2] == n : random number seed (set to n)
  PetscInt options[3] = {0, 0, 15};
  PetscInt numflag = 0;  // C-style numbering (zero based)
  PetscInt *order = new PetscInt[LocalSize];
  PetscInt *locsizes = new PetscInt[2 * PetscGlobalSize];

  ParMETIS_V3_NodeND(vtxdist, xadj, adjc, &numflag, options, order, locsizes,
                     &PETSC_COMM_WORLD);

  // --- permute the matrix
  AO myAO;
  AOCreateBasic(PETSC_COMM_WORLD, LocalSize, original, order, &myAO);
  permuteMatrixEntries(myMesh, myAO);
  AODestroy(&myAO);

  // --- free memory
  delete[] vtxdist;
  delete[] xadj;
  delete[] adjc;
  delete[] original;
  delete[] order;
  delete[] locsizes;

  trace("  finished !");
#endif
}

void cAnalysis::renumberMatrixUsingBoost(
    cMesh &myMesh, const PetscInt &FirstRow, const PetscInt &LastRow,
    std::vector<std::set<PetscInt> > &Pattern) {
#ifndef FRINK_HAVE_BOOST
  trace(" WARNING: matrix won't be renumbered using BOOST");
  trace("          because BOOST is not installed on this system.");
  return;
#else

  // --- some boost::graph types used in the following lines
  typedef boost::adjacency_list<
      boost::vecS, boost::vecS, boost::undirectedS,
      boost::property<boost::vertex_color_t, boost::default_color_type,
                      boost::property<boost::vertex_degree_t, int> > >
      Graph;
  typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
  typedef boost::graph_traits<Graph>::vertices_size_type size_type;
  typedef std::pair<PetscInt, PetscInt> Pair;

  std::vector<Vertex> inv_perm(m_NumberOfUnknowns);
  PetscInt *permutation = new PetscInt[m_NumberOfUnknowns];
  PetscInt *original = new PetscInt[m_NumberOfUnknowns];

  if (PetscGlobalRank == 0) {
    // --- copy my pattern to boost::graph
    std::set<PetscInt>::iterator itPattern;

    Graph myGraph(Pattern.size());
    for (PetscInt row = 0; row < Pattern.size(); row++) {
      for (itPattern = Pattern[row].begin(); itPattern != Pattern[row].end();
           itPattern++) {
        if (row != *itPattern) {
          Pair p(row, *itPattern);
          boost::add_edge(p.first, p.second, myGraph);
        }
      }
    }

    // --- get bandwidth of original matrix
    boost::graph_traits<Graph>::vertex_iterator ui, ui_end;

    boost::property_map<Graph, boost::vertex_degree_t>::type deg =
        get(boost::vertex_degree, myGraph);
    for (boost::tie(ui, ui_end) = vertices(myGraph); ui != ui_end; ++ui)
      deg[*ui] = degree(*ui, myGraph);

    boost::property_map<Graph, boost::vertex_index_t>::type index_map =
        get(boost::vertex_index, myGraph);

    message("  reorder matrix using Cuthill-McKee ...\n");
    message("    original bandwidth : %d\n", boost::bandwidth(myGraph));

    // --- Reverse Cuthill-McKee-ordering
    boost::cuthill_mckee_ordering(myGraph, inv_perm.rbegin(),
                                  boost::get(boost::vertex_color, myGraph),
                                  boost::make_degree_map(myGraph));
    for (size_type c = 0; c != inv_perm.size(); ++c) {
      original[c] = c;
      permutation[index_map[inv_perm[c]]] = c;
    }

    message(
        "    optimized bandwidth: %d\n",
        boost::bandwidth(myGraph, boost::make_iterator_property_map(
                                      permutation, index_map, permutation[0])));
  } else {
    for (PetscInt k = 0; k < m_NumberOfUnknowns; k++) original[k] = k;
  }

  // --- spread the new permutation computed on rank 0 to all other processes
  MPI_Bcast(permutation, m_NumberOfUnknowns, MPIU_INT, 0, PETSC_COMM_WORLD);

  // --- permute the matrix
  AO myAO;
  AOCreateBasic(PETSC_COMM_WORLD, LastRow - FirstRow, &(original[FirstRow]),
                &(permutation[FirstRow]), &myAO);
  permuteMatrixEntries(myMesh, myAO);

  delete[] permutation;
  delete[] original;
  AODestroy(&myAO);

#endif
}

void cAnalysis::permuteMatrixEntries(cMesh &MyMesh, AO &ordering) {
  PetscInt index;

  trace("  permute matrix entries ...");

  for (ItMapNodes itN = MyMesh.getFirstNode(); itN != MyMesh.getLastNode();
       itN++) {
    for (int dof = 0; dof < cstNumberOfKnownDofs; dof++) {
      index = itN->second->getGlobalRow(dof);
      ierr = AOApplicationToPetsc(ordering, 1, &index);
      INFAMCHKERRQ(ierr);
      // AOPetscToApplication(myAO, 1, &index);
      itN->second->setGlobalRow(dof, index);
    }
  }
}

void cAnalysis::preallocateMatrix(const PetscInt &FirstRow,
                                  const PetscInt &LastRow,
                                  std::vector<std::set<PetscInt> > &Pattern) {
  trace("  preallocate matrix ...");
  PetscInt maxNonZerosperRow = 1;
  std::vector<PetscInt> Ii, J;
  Ii.resize(Pattern.size() + 1);
  Ii[0] = 0;
  for (PetscInt k = 0, n = Pattern.size(); k < n; k++) {
    Pattern[k].insert(k);
    Ii[k + 1] = Ii[k] + Pattern[k].size();
    J.insert(J.end(), Pattern[k].begin(), Pattern[k].end());
    if (int(Pattern[k].size()) > maxNonZerosperRow) {
      maxNonZerosperRow = int(Pattern[k].size());
    }
  }
  MatMPIAIJSetPreallocation(m_K, maxNonZerosperRow, NULL, maxNonZerosperRow,
                            NULL);
  MatSeqAIJSetPreallocation(m_K, maxNonZerosperRow, NULL);
  MatSetOption(m_K, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
  MatSetOption(m_K, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
}

void cAnalysis::getLocalSparsityPattern(
    cMesh &myMesh, const PetscInt &FirstRow, const PetscInt &LastRow,
    std::vector<std::set<PetscInt> > &Pattern) {
  trace("  determine local sparsity pattern of matrix ...");

  // --- determine pattern for normal elements
  trace("    adding general finite elements ...");
  for (ItMapElements itE = myMesh.getFirstElement();
       itE != myMesh.getLastElement(); itE++) {
    std::vector<PetscInt> RowsCols(itE->second->getGlobalPositions());
    for (int z = 0; z < (int)RowsCols.size(); z++) {
      if ((RowsCols[z] >= FirstRow) && (RowsCols[z] < LastRow))
        Pattern[RowsCols[z] - FirstRow].insert(RowsCols.begin(),
                                               RowsCols.end());
    }
  }

  // --- determine pattern for interface elements
  trace("    adding fsi elements ...");
  for (ItMapInterface itF = myMesh.getFirstInterfaceElement();
       itF != myMesh.getLastInterfaceElement(); itF++) {
    std::vector<eKnownDofs> DofsStructure(itF->second->getDofsStructure());

    const int ndofs = (int)DofsStructure.size();
    const int nnod = itF->second->getNumberOfNodes();

    // --- collect the rows and cols where this interface element
    //     will generate entries in the global matrices
    std::vector<PetscInt> Rows(
        ndofs *
        nnod);  // rows in global matrix corresponding to structural dofs
    std::vector<PetscInt> Cols(
        nnod);  // rows in global matrix corresponding to fluid dofs
    for (int k = 0; k < nnod; k++) {
      for (int d = 0; d < (int)DofsStructure.size(); d++) {
        Rows[k * ndofs + d] =
            itF->second->getMatchingNode(k)->getGlobalRow(DofsStructure[d]);
      }
      Cols[k] = itF->second->getNode(k)->getGlobalRow(fluid);
    }

    // --- add the interface element to the matrix pattern
    // --- insert C
    for (int z = 0; z < (int)Rows.size(); z++) {
      if ((Rows[z] >= FirstRow) && (Rows[z] < LastRow)) {
        for (int s = 0; s < (int)Cols.size(); s++) {
          Pattern[Rows[z] - FirstRow].insert(Cols[s]);
        }
      }
    }

    // --- insert C^T
    for (int z = 0; z < (int)Cols.size(); z++) {
      if ((Cols[z] >= FirstRow) && (Cols[z] < LastRow)) {
        for (int s = 0; s < (int)Rows.size(); s++) {
          // if (Rows[s] > 0)
          Pattern[Cols[z] - FirstRow].insert(Rows[s]);
        }
      }
    }
  }
}

void cAnalysis::computeStressesCells(cMesh &MyMesh) {
  trace("  compute stresses ...");
  // --- collect solution on all ranks
  Vec FullVector;
  VecScatter ctx;

  ierr = VecScatterCreateToAll(m_x, &ctx, &FullVector);
  INFAMCHKERRQ(ierr);
  ierr = VecScatterBegin(ctx, m_x, FullVector, INSERT_VALUES, SCATTER_FORWARD);
  INFAMCHKERRQ(ierr);
  ierr = VecScatterEnd(ctx, m_x, FullVector, INSERT_VALUES, SCATTER_FORWARD);
  INFAMCHKERRQ(ierr);
  ierr = VecScatterDestroy(&ctx);
  INFAMCHKERRQ(ierr);

  // --- compute stresses for each element
  //     and add them to global vector
  VecZeroEntries(m_StressesCells);
  VecZeroEntries(m_StressesCellsSec2);
  for (ItMapElements itE = MyMesh.getFirstElement();
       itE != MyMesh.getLastElement(); itE++) {
    itE->second->computeStressesCells(FullVector, m_StressesCells,
                                      m_StressesCellsSec2);
  }

  ierr = VecAssemblyBegin(m_StressesCells);
  INFAMCHKERRQ(ierr);
  ierr = VecAssemblyEnd(m_StressesCells);
  INFAMCHKERRQ(ierr);
  ierr = VecAssemblyBegin(m_StressesCellsSec2);
  INFAMCHKERRQ(ierr);
  ierr = VecAssemblyEnd(m_StressesCellsSec2);
  INFAMCHKERRQ(ierr);
  // --- free memory
  ierr = VecDestroy(&FullVector);
  INFAMCHKERRQ(ierr);
}

void cAnalysis::computeStressesNodes(cMesh &MyMesh) {
  // --- collect solution on all ranks
  Vec FullVector;
  VecScatter ctx;

  ierr = VecScatterCreateToAll(m_x, &ctx, &FullVector);
  INFAMCHKERRQ(ierr);
  ierr = VecScatterBegin(ctx, m_x, FullVector, INSERT_VALUES, SCATTER_FORWARD);
  INFAMCHKERRQ(ierr);
  ierr = VecScatterEnd(ctx, m_x, FullVector, INSERT_VALUES, SCATTER_FORWARD);
  INFAMCHKERRQ(ierr);
  ierr = VecScatterDestroy(&ctx);
  INFAMCHKERRQ(ierr);

  // --- compute stresses for each element at the nodes
  //     and add them to global vector
  VecZeroEntries(m_StressesNodes);
  for (ItMapElements itE = MyMesh.getFirstElement();
       itE != MyMesh.getLastElement(); itE++) {
    itE->second->computeStressesNodes(FullVector, m_StressesNodes);
  }

  ierr = VecAssemblyBegin(m_StressesNodes);
  INFAMCHKERRQ(ierr);
  ierr = VecAssemblyEnd(m_StressesNodes);
  INFAMCHKERRQ(ierr);

  // --- free memory
  ierr = VecDestroy(&FullVector);
  INFAMCHKERRQ(ierr);
}

cElementVector cAnalysis::getStressesVanMieses(cMesh &MyMesh) {
  trace("    cAnalysis::getStressesVanMieses ...");
  Vec FullStresses;
  PetscInt elementId0n;
  PetscInt pos;
  cElementFEM *ptrElement = NULL;

  if (m_StressesCells != NULL) {
    VecScatter ctx;
    VecScatterCreateToAll(m_StressesCells, &ctx, &FullStresses);
    VecScatterBegin(ctx, m_StressesCells, FullStresses, INSERT_VALUES,
                    SCATTER_FORWARD);
    VecScatterEnd(ctx, m_StressesCells, FullStresses, INSERT_VALUES,
                  SCATTER_FORWARD);
    VecScatterDestroy(&ctx);
  }
  //---------------------------------------
  // stress
  //---------------------------------------
  cElementVector sigma(6);  // local stresses
  PetscInt size = 0;
  VecGetSize(FullStresses, &size);

  // create vector to store von-Mises stresses
  cElementVector sigma_v(size / 9);

  // loop over elements
  for (ItMapElements it = MyMesh.getFirstElement();
       it != MyMesh.getLastElement(); it++) {
    ptrElement = it->second;
    elementId0n = ptrElement->getId0n();

    for (int i = 0; i < 6; ++i) {
      pos = (elementId0n)*9 + i;
      VecGetValues(FullStresses, 1, &pos, &sigma[i]);
    }

    sigma_v[elementId0n] = std::sqrt(
        sigma[0] * sigma[0] + sigma[1] * sigma[1] + sigma[2] * sigma[2] -
        sigma[0] * sigma[1] - sigma[0] * sigma[2] - sigma[1] * sigma[2] +
        3. * (sigma[3] * sigma[3] + sigma[4] * sigma[4] + sigma[5] * sigma[5]));
  }
  VecDestroy(&FullStresses);
  ptrElement = NULL;
  return sigma_v;
}

cElementVector cAnalysis::getStressesPlaneStress(cMesh &MyMesh) {
  trace("    cAnalysis::getStressesPlaneStress ...");
  Vec FullStresses;

  if (m_StressesCells != NULL) {
    VecScatter ctx;
    VecScatterCreateToAll(m_StressesCells, &ctx, &FullStresses);
    VecScatterBegin(ctx, m_StressesCells, FullStresses, INSERT_VALUES,
                    SCATTER_FORWARD);
    VecScatterEnd(ctx, m_StressesCells, FullStresses, INSERT_VALUES,
                  SCATTER_FORWARD);
    VecScatterDestroy(&ctx);
  }
  //---------------------------------------
  // stress
  //---------------------------------------
  cElementVector sigma(6);  // local stresses
  PetscInt size = 0;
  VecGetSize(FullStresses, &size);

  // create vector to store plane stresses
  cElementVector sigma_p(size / 9);

#ifdef PETSC_USE_COMPLEX
  // to do...
#else
  PetscInt elementId0n;
  PetscInt pos;
  cElementFEM *ptrElement = NULL;
  PetscReal a, b, c;       // Coefficients of char. polynomial
  PetscReal p, q, phi, r;  // dummy variables
  // PetscReal s2imag,s3imag,u,v;
  PetscReal s1, s2, s3;  // plane stresses

  // loop over elements
  for (ItMapElements it = MyMesh.getFirstElement();
       it != MyMesh.getLastElement(); it++) {
    ptrElement = it->second;
    elementId0n = ptrElement->getId0n();

    for (int i = 0; i < 6; ++i) {
      pos = (elementId0n)*9 + i;
      VecGetValues(FullStresses, 1, &pos, &sigma[i]);
    }

    a = -sigma[0] - sigma[1] - sigma[2];
    b = sigma[0] * sigma[1] + sigma[0] * sigma[2] + sigma[1] * sigma[2] -
        sigma[3] * sigma[3] - sigma[4] * sigma[4] - sigma[5] * sigma[5];
    c = sigma[0] * sigma[5] * sigma[5] + sigma[1] * sigma[4] * sigma[4] +
        sigma[2] * sigma[3] * sigma[3] - sigma[0] * sigma[1] * sigma[2] -
        2. * sigma[3] * sigma[4] * sigma[5];

    p = one_over_three * (b - one_over_three * a * a);
    q = 0.5 * (2. * one_over_27 * a * a * a - one_over_three * a * b + c);

    if (q > 0)
      r = std::sqrt(fabs(p));
    else
      r = -std::sqrt(fabs(p));

    if (p < 0) {
      if ((q * q + p * p * p) <= 0) {
        phi = std::acos(q / (r * r * r));
        s1 = -2. * r * std::cos(one_over_three * phi) - one_over_three * a;
        s2 = 2. * r * std::cos(one_over_three * pi - one_over_three * phi) -
             one_over_three * a;
        s3 = 2. * r * std::cos(one_over_three * pi + one_over_three * phi) -
             one_over_three * a;
      } else {
        trace("cAnalysis::getStressesPlaneStress3 - not defined yet!");
        ExitApp();
      }
    } else {
      trace("cAnalysis::getStressesPlaneStress3 - not defined yet!");
      ExitApp();
    }

    sigma_p[elementId0n] = std::sqrt(s1 * s1 + s2 * s2 + s3 * s3);
  }
  ptrElement = NULL;
#endif
  VecDestroy(&FullStresses);
  return sigma_p;
}

void cAnalysis::initializeSolverObjects(void) {
  if (m_AnalysisType == Eigenvalue) {
    trace("  initializing eigenvalue solver with SLEPc... ");
    ierr = EPSCreate(PETSC_COMM_WORLD, &m_Eps);
    INFAMCHKERRQ(ierr);
    ierr = EPSSetOperators(m_Eps, m_K, m_M);
    INFAMCHKERRQ(ierr);

    EPSWhich epsWhich;
    EPSType epsType;

    /* -------------------------------
     * Hard Coded Settings -----------
     * -------------------------------
     *
     * Default: Finds 6 largest eigen values using solver type EPSKRYLOVSCHUR
     *
     * To be extended for: 1. Efficient solving for smallest eigenvalues -
     * KRYLOV methods failed to find eigenvalues which are close to eachother
     *                     2. Generalization with AK3
     *                     3. TARGET Solving
     *                     4. Solving within an interval
     *                     5. Opportunity to choose a solver
     *
     */

    PetscInt solverSetup = 1;

    PetscBool found = PETSC_FALSE;
    PetscOptionsGetInt(PETSC_NULL, "", "-solver", &solverSetup, &found);

    PetscInt numEigPairs = 6;
    PetscOptionsGetInt(PETSC_NULL, "", "-numEigs", &numEigPairs, &found);

    switch (solverSetup) {
      case 1:
        message(
            "  initializing default solve for %i smallest eigenpairs (default "
            "/ case 1).\n",
            numEigPairs);
        trace("     WARNING!! This solver setting is slower!\n");
        epsType = EPSLAPACK;
        epsWhich = EPS_SMALLEST_MAGNITUDE;
        break;
      case 2:
        message("  initializing solve for %i largest eigenpairs (case 2).",
                numEigPairs);
        epsWhich = EPS_LARGEST_MAGNITUDE;
        epsType = EPSKRYLOVSCHUR;
        break;
      case 3:
        message(
            "  initializing default solve for %i largest eigenpairs (case 3).",
            numEigPairs);
        epsWhich = EPS_SMALLEST_MAGNITUDE;
        epsType = EPSKRYLOVSCHUR;
        break;
      default:
        message(
            "  initializing solve for %i smallest eigenpairs (default / case "
            "1). \n",
            numEigPairs);
        trace("     WARNING!! This solver setting is slower!\n");
        epsType = EPSLAPACK;
        epsWhich = EPS_SMALLEST_MAGNITUDE;
    }

    ierr = EPSSetWhichEigenpairs(m_Eps, epsWhich);
    INFAMCHKERRQ(ierr);
    ierr = EPSSetType(m_Eps, epsType);
    INFAMCHKERRQ(ierr);
    ierr = EPSSetDimensions(m_Eps, numEigPairs, PETSC_DEFAULT, PETSC_DEFAULT);
    INFAMCHKERRQ(ierr);
    ierr = EPSSetFromOptions(m_Eps);
    INFAMCHKERRQ(ierr);

  } else {
    ierr = KSPCreate(PETSC_COMM_WORLD, &m_KSP);
    INFAMCHKERRQ(ierr);
    // ierr = KSPSetOperators(m_KSP, m_K, m_K, DIFFERENT_NONZERO_PATTERN);
    // INFAMCHKERRQ(ierr);
    ierr = KSPGetPC(m_KSP, &m_PC);
    INFAMCHKERRQ(ierr);

    PetscInt solverSetup = 1;
    PetscBool found = PETSC_FALSE;
    PetscOptionsGetInt(PETSC_NULL, "", "-solver", &solverSetup, &found);

    switch (solverSetup) {
      case 1:
        trace(
            "  initializing direct mumps solver (OpenSource default / case "
            "1) ");
        ierr = KSPSetType(m_KSP, KSPPREONLY);
        INFAMCHKERRQ(ierr);
        ierr = PCSetType(m_PC, PCLU);
        INFAMCHKERRQ(ierr);
        ierr = PCFactorSetMatSolverType(m_PC, MATSOLVERMUMPS);
        INFAMCHKERRQ(ierr);
        break;
      case 2:
        trace(
            "  initializing direct mumps solver with ou-of-core option (case "
            "2) ");
        ierr = KSPSetType(m_KSP, KSPPREONLY);
        INFAMCHKERRQ(ierr);
        ierr = PCSetType(m_PC, PCLU);
        INFAMCHKERRQ(ierr);
        ierr = PCFactorSetMatSolverType(m_PC, MATSOLVERMUMPS);
        INFAMCHKERRQ(ierr);
        break;
      case 3:
        trace(
            "  initializing petsc solver by command line options (case 3). If "
            "no options are given, petsc's default lu is taken. ");
        ierr = KSPSetType(m_KSP, KSPPREONLY);
        INFAMCHKERRQ(ierr);
        ierr = PCSetType(m_PC, PCLU);
        INFAMCHKERRQ(ierr);
        ierr = PCSetFromOptions(m_PC);
        INFAMCHKERRQ(ierr);
        ierr = KSPSetFromOptions(m_KSP);
        INFAMCHKERRQ(ierr);
        break;
      case 4:
        trace("  initializing direct MKL PARDISO solver (SMP default/case 4) ");
        ierr = KSPSetType(m_KSP, KSPPREONLY);
        INFAMCHKERRQ(ierr);
        ierr = PCSetType(m_PC, PCLU);
        INFAMCHKERRQ(ierr);
        ierr = PCFactorSetMatSolverType(m_PC, MATSOLVERMKL_PARDISO);
        INFAMCHKERRQ(ierr);
        break;
      case 5:
        trace(
            "  initializing direct MKL *C*PARDISO solver (D-SMP default/case "
            "5) ");
        ierr = KSPSetType(m_KSP, KSPPREONLY);
        INFAMCHKERRQ(ierr);
        ierr = PCSetType(m_PC, PCLU);
        INFAMCHKERRQ(ierr);
        ierr =  PCFactorSetMatSolverType(m_PC, MATSOLVERMKL_CPARDISO );
        INFAMCHKERRQ(ierr);
        break;
      default:
        trace("  initializing direct mumps solver (default / case 1) ");
        ierr = KSPSetType(m_KSP, KSPPREONLY);
        INFAMCHKERRQ(ierr);
        ierr = PCSetType(m_PC, PCLU);
        INFAMCHKERRQ(ierr);
        ierr = PCFactorSetMatSolverType(m_PC, MATSOLVERMUMPS);
        INFAMCHKERRQ(ierr);
        break;
    }
  }
}

void cAnalysis::deleteSolverObjects(void) {
  if (m_AnalysisType == Eigenvalue) {
    ierr = EPSDestroy(&m_Eps);
    INFAMCHKERRQ(ierr);
  } else {
    ierr = KSPDestroy(&m_KSP);
    INFAMCHKERRQ(ierr);
  }
}

void cAnalysis::insertElementMatrix(cElementFEM *ptrElement, cElementMatrix &EM,
                                    Mat &SysMatrix) {
  std::vector<PetscInt> GlobalRows(ptrElement->getGlobalPositions());
  PetscInt n = (PetscInt)GlobalRows.size();

  ierr = MatSetValues(SysMatrix, n, &(GlobalRows[0]), n, &(GlobalRows[0]),
                      &(EM(0, 0)), ADD_VALUES);
  INFAMCHKERRQ(ierr);
}

void cAnalysis::insertElementVector(cElementFEM *ptrElement, cElementVector &V,
                                    Vec &SysVector) {
  std::vector<PetscInt> GlobalRows(ptrElement->getGlobalPositions());
  PetscInt n = (PetscInt)GlobalRows.size();
  ierr = VecSetValues(SysVector, n, &(GlobalRows[0]), &(V[0]), ADD_VALUES);
  INFAMCHKERRQ(ierr);
}

void cAnalysis::dumpMatrixToFile(Mat &SysMatrix, const std::string &Filename,
                                 const std::string &Identifier) {
  PetscViewer viewing;

  message("  writing global matrix to textfile ... ");

  ierr = MatAssemblyBegin(SysMatrix, MAT_FINAL_ASSEMBLY);
  INFAMCHKERRQ(ierr);
  ierr = MatAssemblyEnd(SysMatrix, MAT_FINAL_ASSEMBLY);
  INFAMCHKERRQ(ierr);

  PetscObjectSetName((PetscObject)SysMatrix, Identifier.c_str());
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, Filename.c_str(), &viewing);
  PetscViewerPushFormat(viewing, PETSC_VIEWER_ASCII_MATLAB);
  MatView(SysMatrix, viewing);
  PetscViewerPopFormat(viewing);
  PetscViewerDestroy(&viewing);

  trace("  finished!");
}

void cAnalysis::dumpVectorToFile(Vec &SysVector, const std::string &Filename,
                                 const std::string &Identifier) {
  PetscViewer viewing;

  message("  writing global matrix to textfile ... ");

  PetscObjectSetName((PetscObject)SysVector, Identifier.c_str());
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, Filename.c_str(), &viewing);
  PetscViewerPushFormat(viewing, PETSC_VIEWER_ASCII_MATLAB);
  VecView(SysVector, viewing);
  PetscViewerPopFormat(viewing);
  PetscViewerDestroy(&viewing);

  trace("  ready!");
}

void cAnalysis::solveEquationSystem(Mat &matrix, Vec &solution, Vec &vec,
                                    bool resetPreconditioner) {
  PetscInt its = 0, M, N;

  MatGetSize(matrix, &M, &N);
  if (!getMute()) message("  solving linear system of %d equations ...\n", M);
  ierr = KSPSetOperators(m_KSP, matrix, matrix);
  INFAMCHKERRQ(ierr);

  // This loop is placed here AND in funtion "initializeSolverObjects" as some
  // setup require KSPSetOperators before (e.g. mumps icntl)
  PetscInt solverSetup = 1;
  PetscBool found = PETSC_FALSE;
  PetscOptionsGetInt(PETSC_NULL, "", "-solver", &solverSetup, &found);

  ierr = PCFactorSetUpMatSolverType(m_PC);
  INFAMCHKERRQ(ierr);
  switch (solverSetup) {
    case 1:
      break;
    case 2:
      Mat F;
      ierr = PCFactorGetMatrix(m_PC, &F);
      INFAMCHKERRQ(ierr);
      ierr = MatMumpsSetIcntl(F, 22, 1);
      INFAMCHKERRQ(ierr);
      break;
    case 3:
      break;
    default:
      break;
  }

  ierr = KSPSolve(m_KSP, vec, solution);
  INFAMCHKERRQ(ierr);
  ierr = KSPGetIterationNumber(m_KSP, &its);
  INFAMCHKERRQ(ierr);

  if (!getMute()) message("    number of iterations : %d\n", its);
}
