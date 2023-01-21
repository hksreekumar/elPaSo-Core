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

#include "analysisfrequencybasic.h"

#include "./../../misc/hdf5/inputh5.h"
#include "./../../misc/hdf5/outputh5.h"

cAnalysisFrequencyBasic::cAnalysisFrequencyBasic() {
  m_AnalysisType = Frequency;
  m_Start = 0.;
  m_Steps = 0;
  m_Increment = 0.;
}

cAnalysisFrequencyBasic::~cAnalysisFrequencyBasic() {
  // empty
}

void cAnalysisFrequencyBasic::generateFrequencies(const std::string& Filename) {
  trace("determine frequencies to compute ... ");

  // --- compose current filename
  std::string FrqName(Filename);
  if (FrqName.rfind(".ak3") == (FrqName.length() - 4))
    FrqName.erase(FrqName.length() - 4);
  FrqName += ".frq";

  std::ifstream FrqFile;
  FrqFile.open(FrqName.c_str());

  if (FrqFile.fail()) {
    trace("    no file with frequencies found - using parameters of inputfile");

    for (int step = 0; step < getNumberOfSteps(); step++) {
      const PetscReal Value =
          getStartFrequency() + getIncrement() * (PetscReal)step;
      m_Frequencies.insert(Value);
    }
  } else {
    message("    file with frequencies found - reading data\n");

    PetscReal value;

    while (!FrqFile.eof()) {
      if (FrqFile >> value) {
        // std::cout << "Compute frequency " << value << " Hz " << std::endl;
        // --- values will be sorted automatically (std::set<>)
        m_Frequencies.insert(value);
      }
    }
  }

  FrqFile.close();
  trace("  finished !");
}

void cAnalysisFrequencyBasic::FullRun(cMesh& MyMesh, PetscBool freemem) {
  activateDofsAtNodes(MyMesh);
  optimizeMatrixPattern(MyMesh);

  generateFrequencies(MyMesh.getFilename());
  writeInfoFile(MyMesh, m_Frequencies.size(), m_NumberOfUnknowns);

  initializeSolverObjects();

#ifdef HAVE_HDF5
  int writehdf5 =
      cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(
          "hdf5", "/Analysis");
  std::string _filename = "eGenOutput_" + MyMesh.getFilename();
  cOutputH5 systemoutput_hdf5;
  if (writehdf5) {
    trace("  opening global solutionhdf5 ...");
    systemoutput_hdf5.setGlobalFileName(_filename);

    if (systemoutput_hdf5.exists_hdf5())
      systemoutput_hdf5.openContainer(ELPASO_H5_READWRITE);
    else
      systemoutput_hdf5.openContainer(ELPASO_H5_READWRITE_FORCE);
  }
#endif

  int step = 0;
  for (m_ItFrequencies = m_Frequencies.begin();
       m_ItFrequencies != m_Frequencies.end(); m_ItFrequencies++) {
    message("\n Step = %d\n", step);

    const PetscReal frequency = *m_ItFrequencies;
    const PetscReal omega = 2. * M_PI * frequency;

    message(
        "  compute elementmatrices for f = %5.3g Hz, omega = %5.3g rad/s... \n",
        frequency, omega);
    std::cout << "Solving frequency: " << frequency << " Hz \n";

    for (ItMapMaterials itM = MyMesh.getFirstMaterial();
         itM != MyMesh.getLastMaterial(); itM++) {
      // --- tell the materials as well as the elements the new
      //     angular frequency; also update visco elastic material
      //     parameters
      itM->second->setOmega(omega);
      // --- first set Omega, then do the update! (Otherwise Omega=0 and results
      // give nan.)
      itM->second->updateMaterial();
    }

    // --- compute elementmatrices
    trace("    compute K and M ... ");
    cLogging::startClock(SetupKM);

    assembleGlobalTensors(MyMesh, 0, omega);

    cLogging::stopClock(SetupKM);

    // --- consider couplingterms
    //     boundaryconditions are also evaluated

    trace("    consider FSI coupling ... ");

    for (ItMapInterface it = MyMesh.getFirstInterfaceElement();
         it != MyMesh.getLastInterfaceElement(); it++) {
      it->second->addContributionFrequencyDomain(omega, m_K, m_F);
    }

    PetscReal factor = 1;  // inserted to comply with call to apply nodal loads
    applyNodalLoads(MyMesh, factor, omega);

    trace("  assembling matrix K ...");
    ierr = MatAssemblyBegin(m_K, MAT_FINAL_ASSEMBLY);
    INFAMCHKERRQ(ierr);
    ierr = MatAssemblyEnd(m_K, MAT_FINAL_ASSEMBLY);
    INFAMCHKERRQ(ierr);
    ierr = VecAssemblyBegin(m_F);
    INFAMCHKERRQ(ierr);
    ierr = VecAssemblyEnd(m_F);
    INFAMCHKERRQ(ierr);
    trace("  finished");
    // --- insert boundary conditions
    getGlobalBCs(MyMesh, omega);
    insertBoundaryConditions();
    // insertConstraints(MyMesh);

    // ---- export assembled matrix to H5
#if (defined HAVE_HDF5)
    PetscInt exportMatrix = 0;
    PetscBool found = PETSC_FALSE;
    PetscOptionsGetInt(PETSC_NULL, "", "-exportSysH5", &exportMatrix, &found);
    if (exportMatrix) {
      trace("  exporting system as H5 ...");
#ifdef PETSC_USE_COMPLEX
      MathLibrary::ExportPetscMatToHDF5Complex(
          m_K, "/DynamicStiffness", "eGenSystem_" + MyMesh.getFilename());
      MathLibrary::ExportPetscVecToHDF5Complex(
          m_F, "/SystemMatrices", "/LoadVector", "/vecFemLoad",
          "eGenSystem_" + MyMesh.getFilename());
#endif
      trace("  finished");
      ExitApp();
    }
#endif

    // --- solve linear system of equations
    cLogging::startClock(Solver);
    solveEquationSystem(m_K, m_x, m_F);
    cLogging::stopClock(Solver);

    // ---- export solution to H5
#if (defined HAVE_HDF5)
#ifdef PETSC_USE_COMPLEX
    writehdf5 = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(
        "hdf5", "/Analysis");
    if (writehdf5) {
      trace("  exporting solution as HDF5 ...");

      std::string _namePrimGroup = "/Solution";
      std::string _nameSecoGroup = "/State";
      std::string _nameTerGroup = "/vecFemStep" + std::to_string(step + 1);
      // Initialize
      std::vector<PetscScalar> vec;

      // Retrieve Information of Petsc Vec
      MathLibrary::RetrieveVectorFromVecComplex(m_x, vec);

      // Export to HDF5
      systemoutput_hdf5.createGroup(_namePrimGroup);
      systemoutput_hdf5.createSecoGroup(_namePrimGroup, _nameSecoGroup);

      PetscInt Istart, Iend;
      VecGetOwnershipRange(m_x, &Istart, &Iend);

      systemoutput_hdf5.appendComplexVector(vec, _namePrimGroup, _nameSecoGroup,
                                            _nameTerGroup, Istart, Iend);
      systemoutput_hdf5.writeStringAttributeToGroup(
          "FreqStep" + std::to_string(step + 1), "/Solution/State",
          std::to_string(frequency));
    }
#endif
#endif

    // --- compute stresses at elements etc.
    if (checkForStressCellsComputation() == true) {
      trace("  computing stresses at cells ...");
      cAnalysis::computeStressesCells(MyMesh);
    } else {
      trace("  no stress computation at cells requested ...");
    }

    // --- compute stresses at the nodes etc.
    if (checkForStressNodesComputation() == true) {
      trace("  computing stresses at the nodes ...");
      cAnalysis::computeStressesNodes(MyMesh);
    } else {
      trace("  no stress computation at nodes requested ...");
    }

    // reconstructConstraints(MyMesh);
    writeResults(MyMesh, step + 1, frequency);

    // --- now we go on to the next frequency step
    step++;
    // --- zeroing global matrices and vectors
    MatZeroEntries(m_K);
    VecZeroEntries(m_F);
    VecZeroEntries(m_x);
  }
#ifdef HAVE_HDF5
  writehdf5 = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(
      "hdf5", "/Analysis");
  if (writehdf5) {
    trace("  exporting node-dof map ...");
    // Export node map
    std::vector<PetscInt> nodeDofMap;
    nodeDofMap.reserve(m_NumberOfUnknowns * 2);
    for (ItMapNodes it = MyMesh.getFirstNode(); it != MyMesh.getLastNode();
         it++)
      for (size_t j = 0; j < cstNumberOfKnownDofs; j++)
        if (it->second->checkIfActive(j) == true) {
          nodeDofMap.push_back(it->second->getId());
          nodeDofMap.push_back(j);
        }
    systemoutput_hdf5.createSecoGroup("/Solution", "/Maps");
    systemoutput_hdf5.appendIntegerDenseMatrix(
        nodeDofMap, m_NumberOfUnknowns, 2, "/Solution", "/Maps", "/mtxDofMap");

    trace("  closing global solutionhdf5 ...");
    systemoutput_hdf5.closeContainer();
  }
#endif

  deleteSolverObjects();

  if (freemem) deletePETScObjects();

  PetscLogDouble Ram = 0;
  ierr = PetscMemoryGetMaximumUsage(&Ram);
  INFAMCHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Max process memory %g M\n",
                     Ram / (1024 * 1024));
  INFAMCHKERRQ(ierr);
  message("Memory Usage:\n");
  message("max. Ram                      = %.2f MB\n", Ram / (1024 * 1024));
}

void cAnalysisFrequencyBasic::assembleGlobalTensors(cMesh& MyMesh, int nmodes,
                                                    const PetscReal omega) {
  if (!getMute()) trace("  compute element matrices ... ");

  // --- collect solution on all ranks
  VecZeroEntries(m_F);
  for (ItMapElements it = MyMesh.getFirstElement();
       it != MyMesh.getLastElement(); it++) {
    int n =
        it->second->getNumberOfNodes() * it->second->getNumberOfDofsPerNode();

    cElementMatrix KM(n, n);
    cElementVector F(n);

    it->second->assembleDynamicStiffnessMatrix(omega, KM);  // KM=(K-omega^2*M)
    it->second->assembleDynamicLoadVector(omega, F, KM);

    insertElementMatrix(it->second, KM, m_K);
    insertElementVector(it->second, F, m_F);
  }

  VecAssemblyBegin(m_F);
  VecAssemblyEnd(m_F);
}

void cAnalysisFrequencyBasic::initializePETScObjects(cMesh& MyMesh) {
  trace(
      "  cAnalysisFrequencyBasic::initializePETScObjects -> initialize global "
      "matrices ... ");

  message("    dimension n = %d\n", m_NumberOfUnknowns);

  int nn = m_NumberOfUnknowns;

  ierr = MatCreate(PETSC_COMM_WORLD, &m_K);
  INFAMCHKERRQ(ierr);
  ierr = MatSetSizes(m_K, PETSC_DECIDE, PETSC_DECIDE, nn, nn);
  INFAMCHKERRQ(ierr);

  ierr = MatSetType(m_K, MATAIJ);
  INFAMCHKERRQ(ierr);
  trace("    using MUMPS matrix ... ");

  // --- check for command line options
  // ierr = MatSetFromOptions(m_K); INFAMCHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)m_K, "K");
  INFAMCHKERRQ(ierr);
  trace("  initialize global vectors ... ");

  // --- solution vector
  ierr = VecCreate(PETSC_COMM_WORLD, &m_x);
  INFAMCHKERRQ(ierr);
  ierr = VecSetSizes(m_x, PETSC_DECIDE, m_NumberOfUnknowns);
  INFAMCHKERRQ(ierr);

  ierr = VecSetFromOptions(m_x);
  INFAMCHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)m_x, "x");
  INFAMCHKERRQ(ierr);

  VecDuplicate(m_x, &m_F);

  ierr = PetscObjectSetName((PetscObject)m_F, "F");
  INFAMCHKERRQ(ierr);

  // --- global stress vector
  if (checkForStressCellsComputation() == true) {
    int nElements = MyMesh.getNumberOfElements();
    int sumElements = 0;
    MPI_Allreduce(&nElements, &sumElements, 1, MPI_INT, MPI_SUM,
                  PETSC_COMM_WORLD);

    const PetscInt nst =
        (cstNumberOfStressDofs)*sumElements;  // sigma
                                              // _11,_22,_33,_12,_13,_23,gradpx,gradpy,gradpz

    ierr = VecCreate(PETSC_COMM_WORLD, &m_StressesCells);
    INFAMCHKERRQ(ierr);
    ierr = VecSetSizes(m_StressesCells, PETSC_DECIDE, nst);
    INFAMCHKERRQ(ierr);
    ierr = VecSetFromOptions(m_StressesCells);
    INFAMCHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)m_StressesCells, "sigmaCells");
    INFAMCHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &m_StressesCellsSec2);
    INFAMCHKERRQ(ierr);
    ierr = VecSetSizes(m_StressesCellsSec2, PETSC_DECIDE, nst);
    INFAMCHKERRQ(ierr);
    ierr = VecSetFromOptions(m_StressesCellsSec2);
    INFAMCHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)m_StressesCellsSec2, "sigmaCells");
    INFAMCHKERRQ(ierr);
  }
  // --- global stress vector at the nodes
  if (checkForStressNodesComputation() == true) {
    int nElements = MyMesh.getNumberOfNodes();
    int sumNodes = 0;
    MPI_Allreduce(&nElements, &sumNodes, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);

    const PetscInt nst =
        (cstNumberOfStressDofs)*sumNodes;  // sigma
                                           // _11,_22,_33,_12,_13,_23,gradpx,gradpy,gradpz

    ierr = VecCreate(PETSC_COMM_WORLD, &m_StressesNodes);
    INFAMCHKERRQ(ierr);
    ierr = VecSetSizes(m_StressesNodes, PETSC_DECIDE, nst);
    INFAMCHKERRQ(ierr);
    ierr = VecSetFromOptions(m_StressesNodes);
    INFAMCHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)m_StressesNodes, "sigmaNodes");
    INFAMCHKERRQ(ierr);
  }

  setOutSolution(&m_x);
  if (checkForStressCellsComputation() == true)
    setOutStressesCells(&m_StressesCells);
  if (checkForStressCellsComputation() == true)
    setOutStressesCellsSec2(&m_StressesCellsSec2);
  if (checkForStressNodesComputation() == true)
    setOutStressesNodes(&m_StressesNodes);
}

void cAnalysisFrequencyBasic::deletePETScObjects(void) {
  trace("  cAnalysisFrequencyBasic deleting PETSc objects ... ");

  ierr = VecDestroy(&m_x);
  INFAMCHKERRQ(ierr);
  ierr = VecDestroy(&m_F);
  INFAMCHKERRQ(ierr);
  ierr = MatDestroy(&m_K);
  INFAMCHKERRQ(ierr);

  if (checkForStressCellsComputation() == true) {
    ierr = VecDestroy(&m_StressesCells);
    INFAMCHKERRQ(ierr);
    ierr = VecDestroy(&m_StressesCellsSec2);
    INFAMCHKERRQ(ierr);
  }
  if (checkForStressNodesComputation() == true) {
    ierr = VecDestroy(&m_StressesNodes);
    INFAMCHKERRQ(ierr);
  }

  trace("  finished!");
}

std::istream& cAnalysisFrequencyBasic::read(std::istream& is) {
  is >> m_Start >> m_Steps >> m_Increment;

  return is;
}

std::ostream& cAnalysisFrequencyBasic::write(std::ostream& os) const {
  return os;
}

std::ostream& cAnalysisFrequencyBasic::writeXml(std::ostream& os) const {
  os << "<Analysis Type=\"frequency\">";
  os << "<start>" << std::showpoint << getStartFrequency() << "</start>";
  os << "<steps>" << getNumberOfSteps() << "</steps>";
  os << "<delta>" << std::showpoint << getIncrement() << "</delta>";
  os << "</Analysis>";

  return os;
}
