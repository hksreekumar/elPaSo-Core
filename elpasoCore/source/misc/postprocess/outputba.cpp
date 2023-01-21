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

#include "outputba.h"

#include "../../element/fluid/elementfluid3d.h"
#include "outputvtk.h"


cOutputBA::cOutputBA() {
  this->m_RevMode = 0;
  this->m_RevLoadIDPs = 0;
  this->m_RevLoadIDPin = 0;
  this->m_RevVolume = -1.;
  this->m_RevFileName.clear();
  this->m_RevTimeFileName.clear();
  this->m_RevDeviceArea = 0.;
  this->m_RevRoomNr = 1;

  this->m_GlobalNumRowsPerRoom[0] = 0;
  this->m_GlobalNumRowsPerRoom[1] = 0;

  this->RevTimeFrequ.resize(0);
  this->RevTimeVal.resize(0);

  this->Press2ValRoomA.resize(0);
  this->Press2ValRoomB.resize(0);
  this->Pin_act_terzsum = 0.;
  this->lastMiddleTerz = 0.;
  this->nr_lastMiddleTerz = 0;

  //  this->m_LevelFileName.clear();
}

cOutputBA::cOutputBA(cOutputBA &other) : cOutputfile(other) {
  this->m_RevMode = other.m_RevMode;
  this->m_RevLoadIDPs = other.m_RevLoadIDPs;
  this->m_RevLoadIDPin = other.m_RevLoadIDPin;
  this->m_RevVolume = other.m_RevVolume;
  this->m_RevFileName.clear();
  this->m_RevTimeFileName.clear();
  this->m_RevDeviceArea = other.m_RevDeviceArea;
  this->m_RevRoomNr = other.m_RevRoomNr;

  //  this->m_RevFileName		= std::string(other.m_RevFileName);
  this->m_GlobalNumRowsPerRoom[0] = other.m_GlobalNumRowsPerRoom[0];
  this->m_GlobalNumRowsPerRoom[1] = other.m_GlobalNumRowsPerRoom[1];

  this->RevTimeFrequ = other.RevTimeFrequ;
  this->RevTimeVal = other.RevTimeVal;

  this->Press2ValRoomA = other.Press2ValRoomA;
  this->Press2ValRoomB = other.Press2ValRoomB;
  this->Pin_act_terzsum = other.Pin_act_terzsum;
  this->lastMiddleTerz = other.lastMiddleTerz;
  this->nr_lastMiddleTerz = 0;

  //  this->m_LevelFileName.clear();
  // this->m_LevelFileName	= std::string(other.m_LevelFileName);
}

cOutputBA::~cOutputBA() {
  if (wantOutput() > 0) {
    m_RevFile.close();
    //		 m_LevelFile.close();
  }
}

void cOutputBA::computeBuildingAcoustics(cMesh &MyMesh, const int &NrStep,
                                         const double CurrentFrequency,
                                         Vec &solution) {
  trace("  computing calculations for Building Acoustics ... ");

  writeHeader(MyMesh);

  // getRevMode = 1: reverberation time sound power method
  // 2: sound pressure level with sound power through sample
  // 3: energy density and intensity in all fluid points
  // 5: reverberation time normal
  // 8: reverberation time trz
  // 9: sound pressure level difference averaging trz
  if (getRevMode() == 1 || getRevMode() == 2 || getRevMode() == 5 ||
      getRevMode() == 8)
    writeSoundPowerStep(MyMesh, CurrentFrequency, solution);

  else if (getRevMode() == 3) {
    // energy density and intensity in all fluid points
    writeEDandIntNodes(MyMesh, NrStep, CurrentFrequency, solution);
  } else if (getRevMode() == 4) {
    // computes the phase shift of pressure and velocity at the nodes
    writePhaseShiftNodes(MyMesh, NrStep, CurrentFrequency, solution);
  } else if (getRevMode() == 6) {
    // computes the SRI with RTC for all frequ
    computePressureLevelRTC(CurrentFrequency, solution, MyMesh);
  } else if (getRevMode() == 7) {
    // computes the SRI with RTC averaging terz frequ
    computePressuresRTCAtNodesTerz(CurrentFrequency, solution, MyMesh);
  } else if (getRevMode() == 9) {
    // computes the SRI without RTC averaging terz frequ
    computePressuresAtNodesTerz(CurrentFrequency, solution, MyMesh);
  } else {
    trace(
        "  !!! Error: this mode for bulding acoustic calculation is not valid "
        "!!! ");
  }

  trace("  Building Acoustics computed! ");
}

void cOutputBA::writeHeader(cMesh &MyMesh) {
  // --- if the stream LevelFile is closed, this function is called
  //     for the first time. So we have to adjust the local range of
  //     rows and open the file.
  if (m_RevFileName.empty()) {
    if (getRevMode() == 1) {
      // -------------------------------------------------------------------------
      //   create file name
      // -------------------------------------------------------------------------
      std::string text = MyMesh.getFilename();
      if (text.rfind(cstInputFileSuffix) ==
          (text.length() - cstInputFileSuffix.length()))
        text.erase(text.length() - cstInputFileSuffix.length());
      m_RevFileName = text + ".rev";

      // Only the first process writes the header
      if (PetscGlobalRank == 0) {
        m_RevFile.open(m_RevFileName.c_str());
        if (m_RevFile.fail()) {
          throw cException("Aborting! Unable to open " + m_RevFileName,
                           __FILE__, __LINE__);
        }

        // --- writing header of file
        m_RevFile << "#--------------------------------------------------------"
                     "-------------------------------------------"
                  << std::endl;
        //			  m_RevFile << "#    f [Hz]       Pin_act[dB]
        //Pin_react[dB]  Ps_act[dB]  Ps_react[dB]  Tr_act[s]   Tr_react[s] "  <<
        //std::endl;
        m_RevFile << "#    f [Hz]       Pin_act[dB]  Ps_act[dB]  Tr_act[s]  "
                  << std::endl;
        m_RevFile << "#--------------------------------------------------------"
                     "-------------------------------------------"
                  << std::endl;

        m_RevFile.setf(std::ios::scientific);
      }
    } else if (getRevMode() == 2) {
      // -------------------------------------------------------------------------
      //   create file name
      // -------------------------------------------------------------------------
      std::string text = MyMesh.getFilename();
      if (text.rfind(cstInputFileSuffix) ==
          (text.length() - cstInputFileSuffix.length()))
        text.erase(text.length() - cstInputFileSuffix.length());
      m_RevFileName = text + ".levP";

      // Only the first process writes the header
      if (PetscGlobalRank == 0) {
        m_RevFile.open(m_RevFileName.c_str());
        if (m_RevFile.fail()) {
          throw cException("Aborting! Unable to open " + m_RevFileName,
                           __FILE__, __LINE__);
        }

        // --- writing header of file
        m_RevFile << "#--------------------------------------------------------"
                     "-------------------------------"
                  << std::endl;
        //			  m_RevFile << "#    f [Hz]       Pin_act[dB]
        //Pin_react[dB]  Ps_act[dB]  Ps_react[dB]  L_act[dB]   	" << std::endl;
        m_RevFile << "#    f [Hz]       Pin_act[dB] Ps_act[dB]  L_act[dB] "
                  << std::endl;
        m_RevFile << "#--------------------------------------------------------"
                     "-------------------------------"
                  << std::endl;

        m_RevFile.setf(std::ios::scientific);
      }

    } else if (getRevMode() == 3) {
    } else if (getRevMode() == 4) {
    } else if (getRevMode() == 5) {
      parseNdsFile(MyMesh);
      checkForSides(MyMesh);

      // -------------------------------------------------------------------------
      //   create file name
      // -------------------------------------------------------------------------
      std::string text = MyMesh.getFilename();
      if (text.rfind(cstInputFileSuffix) ==
          (text.length() - cstInputFileSuffix.length()))
        text.erase(text.length() - cstInputFileSuffix.length());
      m_RevFileName = text + ".revPressure";

      if (PetscGlobalRank == 0) {
        m_RevFile.open(m_RevFileName.c_str());
        if (m_RevFile.fail()) {
          throw cException("Aborting! Unable to open " + m_RevFileName,
                           __FILE__, __LINE__);
        }

        // --- writing header of file
        m_RevFile << "#--------------------------------------------------------"
                  << std::endl;
        m_RevFile << "#    f [Hz]        Pin           p2[dB]         Tr[s]    "
                  << std::endl;
        m_RevFile << "#                  Source        Raum 2                  "
                  << std::endl;
        m_RevFile << "#--------------------------------------------------------"
                  << std::endl;

        if (m_GlobalNumRowsPerRoom[1] <= 0)
          m_RevFile << "#  no nodes specified for Room B !!" << std::endl;

        m_RevFile.setf(std::ios::scientific);
      }
    } else if (getRevMode() == 6) {
      parseNdsFile(MyMesh);
      checkForSides(MyMesh);

      // -------------------------------------------------------------------------
      //   create file name
      // -------------------------------------------------------------------------
      std::string text = MyMesh.getFilename();
      if (text.rfind(cstInputFileSuffix) ==
          (text.length() - cstInputFileSuffix.length()))
        text.erase(text.length() - cstInputFileSuffix.length());
      m_RevFileName = text + ".levRTC";

      if (PetscGlobalRank == 0) {
        m_RevFile.open(m_RevFileName.c_str());
        if (m_RevFile.fail()) {
          throw cException("Aborting! Unable to open " + m_RevFileName,
                           __FILE__, __LINE__);
        }

        // --- writing header of file
        m_RevFile << "#--------------------------------------------------------"
                     "----------------------	"
                  << std::endl;
        m_RevFile << "#    f [Hz]        L1[dB]        L2[dB]      L1-L2[dB]   "
                     "  Tr[s]     RTC[dB]    R[dB] "
                  << std::endl;
        m_RevFile << "#                  Raum 1        Raum 2		"
                     "				                              "
                  << std::endl;
        m_RevFile << "#--------------------------------------------------------"
                     "---------------------- "
                  << std::endl;

        if (m_GlobalNumRowsPerRoom[0] <= 0)
          m_RevFile << "#  no nodes specified for Room A !!" << std::endl;

        if (m_GlobalNumRowsPerRoom[1] <= 0)
          m_RevFile << "#  no nodes specified for Room B !!" << std::endl;

        m_RevFile.setf(std::ios::scientific);
      }
    } else if (getRevMode() == 7) {
      parseNdsFile(MyMesh);
      checkForSides(MyMesh);

      // -------------------------------------------------------------------------
      //   create file name
      // -------------------------------------------------------------------------
      std::string text = MyMesh.getFilename();
      if (text.rfind(cstInputFileSuffix) ==
          (text.length() - cstInputFileSuffix.length()))
        text.erase(text.length() - cstInputFileSuffix.length());
      m_RevFileName = text + ".levRTCtrz";

      if (PetscGlobalRank == 0) {
        m_RevFile.open(m_RevFileName.c_str());
        if (m_RevFile.fail()) {
          throw cException("Aborting! Unable to open " + m_RevFileName,
                           __FILE__, __LINE__);
        }

        // --- writing header of file
        m_RevFile << "#--------------------------------------------------------"
                     "----------------------	"
                  << std::endl;
        m_RevFile << "#    f [Hz]        L1[dB]        L2[dB]      L1-L2[dB]   "
                     " Tr[s]     RTC[dB]    R[dB] "
                  << std::endl;
        m_RevFile << "#                  Raum 1        Raum 2		"
                     "	                                    "
                  << std::endl;
        m_RevFile << "#--------------------------------------------------------"
                     "---------------------- "
                  << std::endl;

        if (m_GlobalNumRowsPerRoom[0] <= 0)
          m_RevFile << "#  no nodes specified for Room A !!" << std::endl;

        if (m_GlobalNumRowsPerRoom[1] <= 0)
          m_RevFile << "#  no nodes specified for Room B !!" << std::endl;

        m_RevFile.setf(std::ios::scientific);
      }

      Press2ValRoomA.resize(m_GlobalNumRowsPerRoom[0]);
      Press2ValRoomB.resize(m_GlobalNumRowsPerRoom[1]);

    } else if (getRevMode() == 8) {
      parseNdsFile(MyMesh);
      checkForSides(MyMesh);

      // -------------------------------------------------------------------------
      //   create file name
      // -------------------------------------------------------------------------
      std::string text = MyMesh.getFilename();
      if (text.rfind(cstInputFileSuffix) ==
          (text.length() - cstInputFileSuffix.length()))
        text.erase(text.length() - cstInputFileSuffix.length());
      m_RevFileName = text + ".revPressuretrz";

      if (PetscGlobalRank == 0) {
        m_RevFile.open(m_RevFileName.c_str());
        if (m_RevFile.fail()) {
          throw cException("Aborting! Unable to open " + m_RevFileName,
                           __FILE__, __LINE__);
        }

        // --- writing header of file
        m_RevFile << "#--------------------------------------------------------"
                  << std::endl;
        m_RevFile << "#    f [Hz]        Pin           p2[dB]         Tr[s]    "
                  << std::endl;
        m_RevFile << "#                  Source        Raum 2                  "
                  << std::endl;
        m_RevFile << "#--------------------------------------------------------"
                  << std::endl;

        if (m_GlobalNumRowsPerRoom[1] <= 0)
          m_RevFile << "#  no nodes specified for Room B !!" << std::endl;

        m_RevFile.setf(std::ios::scientific);

        Press2ValRoomA.resize(m_GlobalNumRowsPerRoom[0]);
        Press2ValRoomB.resize(m_GlobalNumRowsPerRoom[1]);
      }
    } else if (getRevMode() == 9) {
      parseNdsFile(MyMesh);
      checkForSides(MyMesh);

      // -------------------------------------------------------------------------
      //   create file name
      // -------------------------------------------------------------------------
      std::string text = MyMesh.getFilename();
      if (text.rfind(cstInputFileSuffix) ==
          (text.length() - cstInputFileSuffix.length()))
        text.erase(text.length() - cstInputFileSuffix.length());
      m_RevFileName = text + ".levtrz";

      if (PetscGlobalRank == 0) {
        m_RevFile.open(m_RevFileName.c_str());
        if (m_RevFile.fail()) {
          throw cException("Aborting! Unable to open " + m_RevFileName,
                           __FILE__, __LINE__);
        }

        // --- writing header of file
        m_RevFile << "#--------------------------------------------------------"
                     "----------------------	"
                  << std::endl;
        m_RevFile << "#    f [Hz]        L1[dB]        L2[dB]      L1-L2[dB]   "
                     "                       "
                  << std::endl;
        m_RevFile << "#                  Raum 1        Raum 2		"
                     "	                                    "
                  << std::endl;
        m_RevFile << "#--------------------------------------------------------"
                     "---------------------- "
                  << std::endl;

        if (m_GlobalNumRowsPerRoom[0] <= 0)
          m_RevFile << "#  no nodes specified for Room A !!" << std::endl;

        if (m_GlobalNumRowsPerRoom[1] <= 0)
          m_RevFile << "#  no nodes specified for Room B !!" << std::endl;

        m_RevFile.setf(std::ios::scientific);
      }

      Press2ValRoomA.resize(m_GlobalNumRowsPerRoom[0]);
      Press2ValRoomB.resize(m_GlobalNumRowsPerRoom[1]);

    } else {
      trace(
          "  !!! Error: this mode for bulding acoustic calculation is not "
          "valid !!! ");
    }
  }
}

void cOutputBA::writeResultsToFile(cMesh &MyMesh, const int &NrStep,
                                   const double &CurrentFrequency) {
  if (wantOutput() > 0) {
    // --- tell the materials as well as the elements the new
    //     angular frequency; also update visco elastic material
    //     parameters
    cMaterial::setOmega(2. * M_PI * CurrentFrequency);
    for (ItMapMaterials itM = MyMesh.getFirstMaterial();
         itM != MyMesh.getLastMaterial(); itM++) {
      itM->second->updateMaterial();
    }

    Vec FullVector;
    VecScatter ctx;

    ierr = VecScatterCreateToAll(m_OutSolution[0], &ctx, &FullVector);
    INFAMCHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, m_OutSolution[0], FullVector, INSERT_VALUES,
                           SCATTER_FORWARD);
    INFAMCHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, m_OutSolution[0], FullVector, INSERT_VALUES,
                         SCATTER_FORWARD);
    INFAMCHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx);
    INFAMCHKERRQ(ierr);

    computeBuildingAcoustics(MyMesh, NrStep, CurrentFrequency, FullVector);

    // --- free memory
    ierr = VecDestroy(&FullVector);
    INFAMCHKERRQ(ierr);
  }
}

void cOutputBA::writeEDandIntNodes(cMesh &MyMesh, const int &NrStep,
                                   const double CurrentFrequency,
                                   Vec &solution) {
  trace("  calculate energy density and intensity ... ");

  // EnergyDensity vector
  Vec EnergyDensity;
  // --- dublicate solution vector
  VecDuplicate(solution, &EnergyDensity);

  // Zeros to the EnergyDensity vector
  VecZeroEntries(EnergyDensity);

  // intensity vector
  Vec intensity;

  int nElements = MyMesh.getNumberOfNodes();
  int sumNodes = 0;
  MPI_Allreduce(&nElements, &sumNodes, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);

  const PetscInt nst =
      (cstNumberOfStressDofs)*sumNodes;  // sigma
                                         // _11,_22,_33,_12,_13,_23,gradpx,gradpy,gradpz

  ierr = VecCreate(PETSC_COMM_WORLD, &intensity);
  INFAMCHKERRQ(ierr);
  ierr = VecSetSizes(intensity, PETSC_DECIDE, nst);
  INFAMCHKERRQ(ierr);
  ierr = VecSetFromOptions(intensity);
  INFAMCHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)intensity, "intensityNodes");
  INFAMCHKERRQ(ierr);

  // Zeros to the intensity vector
  VecZeroEntries(intensity);

  for (ItMapElements itE = MyMesh.getFirstElement();
       itE != MyMesh.getLastElement(); itE++) {
    // The EnergyDensity calculation is now only available for Fluid3d elements
    cElementFluid3d *ptrE = dynamic_cast<cElementFluid3d *>(itE->second);
    if (ptrE != NULL) {
      // compute the energy density in each element and write it to the
      // corresponding row of the EnergyDensity vector
      ptrE->computeEDandIntNodes(solution, EnergyDensity, intensity);
    }
  }
  ierr = VecAssemblyBegin(EnergyDensity);
  INFAMCHKERRQ(ierr);
  ierr = VecAssemblyEnd(EnergyDensity);
  INFAMCHKERRQ(ierr);
  ierr = VecAssemblyBegin(intensity);
  INFAMCHKERRQ(ierr);
  ierr = VecAssemblyEnd(intensity);
  INFAMCHKERRQ(ierr);

  // std::cout <<" EnergyDensity " << EnergyDensity << std::endl;
  // VecView(EnergyDensity,PETSC_VIEWER_STDOUT_WORLD);

  // VTK - Output of energy density
  cOutputVTK *outp_vtk = new cOutputVTK();
  outp_vtk->setWantOutput(true);

  outp_vtk->setFilenameRootadd("_ED_int");
  outp_vtk->setModifyFilenameRoot(true);

  /// function missused //Marco
  // outp_vtk->writeResultsToFile(MyMesh, NrStep, CurrentFrequency,
  // EnergyDensity, &intensity); workaround: set neu vectors for output! -> ask
  // Marco
  trace(
      "ERROR: cOutputBA::writeEDandIntNodes() function missused! Check the "
      "Code!");
  ExitApp();

  // delete the output object
  delete outp_vtk;
}

void cOutputBA::writePhaseShiftNodes(cMesh &MyMesh, const int &NrStep,
                                     const double CurrentFrequency,
                                     Vec &solution) {
  trace("  calculate phase shift at the nodes ... ");

  // phaseShiftAverage vector
  Vec phaseShiftAverage;
  // --- dublicate solution vector
  VecDuplicate(solution, &phaseShiftAverage);

  // Zeros to the phaseShiftAverage vector
  VecZeroEntries(phaseShiftAverage);

  // intensity vector
  Vec phaseShift;

  int nElements = MyMesh.getNumberOfNodes();
  int sumNodes = 0;
  MPI_Allreduce(&nElements, &sumNodes, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);

  const PetscInt nst =
      (cstNumberOfStressDofs)*sumNodes;  // sigma
                                         // _11,_22,_33,_12,_13,_23,gradpx,gradpy,gradpz

  ierr = VecCreate(PETSC_COMM_WORLD, &phaseShift);
  INFAMCHKERRQ(ierr);
  ierr = VecSetSizes(phaseShift, PETSC_DECIDE, nst);
  INFAMCHKERRQ(ierr);
  ierr = VecSetFromOptions(phaseShift);
  INFAMCHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)phaseShift, "phaseShiftNodes");
  INFAMCHKERRQ(ierr);

  // Zeros to the intensity vector
  VecZeroEntries(phaseShift);

  for (ItMapElements itE = MyMesh.getFirstElement();
       itE != MyMesh.getLastElement(); itE++) {
    // The EnergyDensity calculation is now only available for Fluid3d elements
    cElementFluid3d *ptrE = dynamic_cast<cElementFluid3d *>(itE->second);
    if (ptrE != NULL) {
      // compute the energy density in each element and write it to the
      // corresponding row of the EnergyDensity vector
      ptrE->computePhaseShiftNodes(solution, phaseShiftAverage, phaseShift);
    }
  }
  ierr = VecAssemblyBegin(phaseShiftAverage);
  INFAMCHKERRQ(ierr);
  ierr = VecAssemblyEnd(phaseShiftAverage);
  INFAMCHKERRQ(ierr);
  ierr = VecAssemblyBegin(phaseShift);
  INFAMCHKERRQ(ierr);
  ierr = VecAssemblyEnd(phaseShift);
  INFAMCHKERRQ(ierr);

  // std::cout <<" EnergyDensity " << EnergyDensity << std::endl;
  // VecView(EnergyDensity,PETSC_VIEWER_STDOUT_WORLD);

  // VTK - Output of energy density
  cOutputVTK *outp_vtk = new cOutputVTK();
  outp_vtk->setWantOutput(true);

  outp_vtk->setFilenameRootadd("_PS");
  outp_vtk->setModifyFilenameRoot(true);

  /// function missused //Marco
  // outp_vtk->writeResultsToFile(MyMesh, NrStep, CurrentFrequency,
  // phaseShiftAverage, NULL, &phaseShift); workaround: set neu vectors for
  // output! -> ask Marco
  trace(
      "ERROR: cOutputBA::writePhaseShiftNodes() function missused! Check the "
      "Code!");
  ExitApp();

  // delete the output object
  delete outp_vtk;
}

void cOutputBA::writeSoundPowerStep(cMesh &MyMesh,
                                    const double CurrentFrequency,
                                    Vec &solution) {
  // --- if the file is opened successfully, we'll now compute
  //     the reverberation time - first the sound power at the boundaries P_s
  //     ...
  trace("    computing local sound power P_s/P_in ... ");
  PetscReal Ps_res_act =
      0.;  // sound power at the boundaries based on the active sound intensity
  std::vector<PetscReal> Ps_act_sum;
  //	PetscReal Ps_res_react = 0.;	// sound power at the boundaries based
  //on the reactive sound intensity

  PetscReal Pin_res_act =
      0.;  // sound power at the source based on the active sound intensity
  std::vector<PetscReal> Pin_act_sum;
  //	PetscReal Pin_res_react = 0.;	// sound power at the source based on
  //the reactive sound intensity

  PetscReal vn_in = 0.;           // sound velocity of source
  PetscReal boundary_area = 0.;   // area of the boundaries
  PetscReal speedofsound = 343.;  // speed of sound
  PetscReal rho0 = 1.21;          // fluid density

  // Loop over all Elements
  for (ItMapElements itE = MyMesh.getFirstElement();
       itE != MyMesh.getLastElement(); itE++) {
    // The sound power calculation is now only available for Fluid3d elements
    cElementFluid3d *ptrE = dynamic_cast<cElementFluid3d *>(itE->second);
    if (ptrE != NULL) {
      // Loop over all Loads
      for (ItLoadsOnElementMap it = itE->second->getFirstElementLoad();
           it != itE->second->getLastElementLoad(); it++) {
        // If the LoadIDPs is found
        if (it->second->getId() == m_RevLoadIDPs) {
          PetscReal area = ptrE->getFaceArea(it->first - 1);
          PetscReal Ps_face_act = 0;
          //					PetscReal Ps_face_react = 0;

          // For reverberation time calculation the source velocity is allready
          // known
          if (getRevMode() == 5 || getRevMode() == 8) {
            /*cElementLoadVn *ptrELvn = dynamic_cast<cElementLoadVn
            *>(it->second); if (ptrELvn != NULL) {
              //PetscReal GivenVn = 0.1;
              vn_in = ptrELvn->getValue();
            }
            else
              trace("error !!!  load on element for rev-time computation must be
            of type cElementLoadVn ");
            */
            ptrE->calculateActiveSoundPowerSurfaceIntegralGivenVn(
                solution, it->first - 1, Ps_face_act, vn_in);
          } else {
            ptrE->calculateActiveSoundPowerSurfaceIntegral(
                solution, it->first - 1, Ps_face_act);
            //					ptrE->calculateReactiveSoundPowerSurfaceIntegral(solution,
            //it->first - 1, Ps_face_react);
            //					ptrE->calculateActiveSoundPowerSurfaceAverage(solution,
            //it->first - 1, Ps_face_act);
            //					ptrE->calculateReactiveSoundPowerSurfaceAverage(solution,
            //it->first - 1, Ps_face_react);
          }

          Ps_act_sum.push_back(Ps_face_act);
          //					Ps_res_act   += Ps_face_act;
          //					Ps_res_react += Ps_face_react;
          boundary_area += area;
        }
        // If the LoadIDPin is found
        if (it->second->getId() == m_RevLoadIDPin) {
          // PetscReal area = ptrE->getFaceArea(it->first - 1); // Comment due
          // to "unused variable" PetscReal Pin_face_andereFormel = 0; // Comment
          // due to "unused variable" cMaterialFluid *ptrF =
          // dynamic_cast<cMaterialFluid *>(ptrE->getMaterial());

          PetscReal Pin_face_act = 0;
          //					PetscReal Pin_face_react = 0;

          // For reverberation time calculation the source velocity is allready
          // known
          if (getRevMode() == 5 || getRevMode() == 8) {
            /*cElementLoadVn *ptrELvn = dynamic_cast<cElementLoadVn
            *>(it->second); if (ptrELvn != NULL) {
              //PetscReal GivenVn = 0.1;
              vn_in = ptrELvn->getValue();
            }
            else
              trace("error !!!  load on element for rev-time computation must be
            of type cElementLoadVn ");
            */
            ptrE->calculateActiveSoundPowerSurfaceIntegralGivenVn(
                solution, it->first - 1, Pin_face_act, vn_in);
          } else {
            ptrE->calculateActiveSoundPowerSurfaceIntegral(
                solution, it->first - 1, Pin_face_act);
            //					ptrE->calculateReactiveSoundPowerSurfaceIntegral(solution,
            //it->first - 1, Pin_face_react);
            //					ptrE->calculateActiveSoundPowerSurfaceAverage(solution,
            //it->first - 1, Pin_face_act);
            //					ptrE->calculateReactiveSoundPowerSurfaceAverage(solution,
            //it->first - 1, Pin_face_react);
          }

          Pin_act_sum.push_back(Pin_face_act);
          //					Pin_res_act   += Pin_face_act;
          //					Pin_res_react += Pin_face_react;
        }
      }
    }
  }

  std::sort(Ps_act_sum.begin(), Ps_act_sum.end(), abs_order);
  std::sort(Pin_act_sum.begin(), Pin_act_sum.end(), abs_order);

  for (int i = Pin_act_sum.size() - 1; i >= 0; i--)
    Pin_res_act += Pin_act_sum[i];

  for (int i = Ps_act_sum.size() - 1; i >= 0; i--) Ps_res_act += Ps_act_sum[i];

  // --- now we compute the sum over all processes and gather the result on rank
  // 0
  trace("    computing global sound power ... ");
  PetscReal globalPin_act = 0.;
  //	  PetscReal globalPin_react		= 0.;
  PetscReal globalPs_act = 0.;
  //	  PetscReal globalPs_react		= 0.;
  PetscReal globalboundary_area = 0.;
  MPI_Reduce(&Pin_res_act, &globalPin_act, 1, MPIU_REAL, MPI_SUM, 0,
             PETSC_COMM_WORLD);
  //	  MPI_Reduce(&Pin_res_react,	&globalPin_react,	1, MPIU_REAL,
  //MPI_SUM, 0, PETSC_COMM_WORLD);
  MPI_Reduce(&Ps_res_act, &globalPs_act, 1, MPIU_REAL, MPI_SUM, 0,
             PETSC_COMM_WORLD);
  //	  MPI_Reduce(&Ps_res_react	,	&globalPs_react,	1, MPIU_REAL,
  //MPI_SUM, 0, PETSC_COMM_WORLD);
  MPI_Reduce(&boundary_area, &globalboundary_area, 1, MPIU_REAL, MPI_SUM, 0,
             PETSC_COMM_WORLD);

  // --- output of the computed reverberation times
  if (PetscGlobalRank == 0) {
    if (getRevMode() == 1) {
      PetscReal Tr_act = 0.;  // reverberation time with active mehthod
      //			PetscReal Tr_react  = 0.;     // reverberation
      //time with reactive method Schallleistung

      // The facenormal is in wrong direction, so the Power has the opposite
      // sign
      //			Tr_act   = 55.262 * m_RevVolume / (
      //globalPin_act * globalboundary_area * speedofsound) * globalPs_act;
      Tr_act = 55.262 * m_RevVolume /
               (-1 * globalPin_act * globalboundary_area * speedofsound) *
               globalPs_act;
      //			Tr_react = 55.262 * m_RevVolume / ( -1*
      //globalPin_react * globalboundary_area * speedofsound) * globalPs_react;

      m_RevFile << std::setw(15) << CurrentFrequency;
      m_RevFile << std::setw(15) << globalPin_act;
      //			m_RevFile << std::setw(15)	<<
      //globalPin_react;
      m_RevFile << std::setw(15) << globalPs_act;
      //			m_RevFile << std::setw(15)	<<
      //globalPs_react;
      m_RevFile << std::setw(15) << Tr_act;
      //			m_RevFile << std::setw(15)	<< Tr_react;
      m_RevFile << std::endl << std::flush;
    } else if (getRevMode() == 2) {
      PetscReal L_act = 0.;  // transmission loss with active mehthod
      //			PetscReal L_react	= 0.;     //
      //transmission loss with reactive method

      // The facenormal is in wrong direction, so the Power Ps has the opposite
      // sign
      L_act = 10 * std::log10(globalPin_act / (-1 * globalPs_act));
      //			L_react = 10 * std::log10(globalPin_react / (-1
      //* globalPs_react));

      // std::cout <<" CurrentFrequency  " << CurrentFrequency << std::endl;

      m_RevFile << std::setw(15) << CurrentFrequency;
      m_RevFile << std::setw(15) << globalPin_act;
      //			m_RevFile << std::setw(15)	<<
      //globalPin_react;
      m_RevFile << std::setw(15) << -1 * globalPs_act;
      //			m_RevFile << std::setw(15)	<< -1 *
      //globalPs_react;
      m_RevFile << std::setw(15) << L_act;
      //			m_RevFile << std::setw(15)	<< L_react;
      m_RevFile << std::endl << std::flush;
    } else if (getRevMode() == 5) {
      PetscReal Tr_act = 0.;  // reverberation time
      PetscReal globalSumRoomA = 0.;
      PetscReal globalSumRoomB = 0.;

      computePressures_in_Room(CurrentFrequency, solution, MyMesh,
                               globalSumRoomA, globalSumRoomB);

      // computation on LUDWIG cluster leads to opposite sign in Pin calculation
      if (globalPin_act < 0.0) globalPin_act = globalPin_act * (-1);

      if (getRevRoomNr() == 0)
        Tr_act = 55.262 * m_RevVolume /
                 (globalPin_act * speedofsound * speedofsound * 4 * rho0) *
                 globalSumRoomA;
      else
        Tr_act = 55.262 * m_RevVolume /
                 (globalPin_act * speedofsound * speedofsound * 4 * rho0) *
                 globalSumRoomB;

      m_RevFile << std::setw(15) << CurrentFrequency;
      m_RevFile << std::setw(15) << globalPin_act;
      if (getRevRoomNr() == 0)
        m_RevFile << std::setw(15) << globalSumRoomA;
      else
        m_RevFile << std::setw(15) << globalSumRoomB;
      m_RevFile << std::setw(15) << Tr_act;
      m_RevFile << std::endl << std::flush;
    }
  }

  if (getRevMode() == 8) {
    std::vector<PetscReal> TerzMitte;
    std::vector<PetscReal> TerzOben;
    std::vector<PetscReal> TerzUnten;
    getTerzFrequencies(TerzMitte, TerzOben, TerzUnten);

    PetscReal currentMiddleTerz;
    const PetscInt numNodesRoomA = m_GlobalNumRowsPerRoom[0];
    const PetscInt numNodesRoomB = m_GlobalNumRowsPerRoom[1];

    // computation on LUDWIG cluster leads to opposite sign in Pin calculation
    if (globalPin_act < 0.0) globalPin_act = globalPin_act * (-1);

    // -------------------------------------------------------------------------
    //   First Averaging the Reverberation time data into thirds,
    //   than averaging over the room
    // -------------------------------------------------------------------------

    for (int k = 0; k < (int)TerzMitte.size(); k++) {
      if ((TerzUnten[k] <= CurrentFrequency) &&
          (CurrentFrequency < TerzOben[k])) {
        if (lastMiddleTerz == 0.0) lastMiddleTerz = TerzMitte[k];

        currentMiddleTerz = TerzMitte[k];

        if (std::abs(currentMiddleTerz - lastMiddleTerz) < cstGeomEps) {
          computePressuresAtNodesFrequ_Room(solution, MyMesh, Press2ValRoomA,
                                            Press2ValRoomB);
          Pin_act_terzsum = Pin_act_terzsum + globalPin_act;

          nr_lastMiddleTerz = nr_lastMiddleTerz + 1;
          lastMiddleTerz = currentMiddleTerz;
        } else {
          // --- now we compute the sum over all processes and gather the result
          // on rank 0
          trace("    computing global sum for Terz... ");
          PetscReal pEff2SumRoomA = 0.;
          PetscReal pEff2SumRoomB = 0.;

          for (int n = 0; n < (int)Press2ValRoomA.size(); n++) {
            pEff2SumRoomA = pEff2SumRoomA + Press2ValRoomA[n];
          }
          for (int n = 0; n < (int)Press2ValRoomB.size(); n++) {
            pEff2SumRoomB = pEff2SumRoomB + Press2ValRoomB[n];
          }

          PetscReal globalSumRoomA = 0.;
          PetscReal globalSumRoomB = 0.;

          MPI_Reduce(&pEff2SumRoomA, &globalSumRoomA, 1, MPIU_REAL, MPI_SUM, 0,
                     PETSC_COMM_WORLD);
          MPI_Reduce(&pEff2SumRoomB, &globalSumRoomB, 1, MPIU_REAL, MPI_SUM, 0,
                     PETSC_COMM_WORLD);

          if (PetscGlobalRank == 0) {
            PetscReal Tr_act = 0;

            if (getRevRoomNr() == 0)
              Tr_act =
                  55.262 * m_RevVolume /
                  (Pin_act_terzsum * speedofsound * speedofsound * 4 * rho0) *
                  globalSumRoomA / ((PetscReal)numNodesRoomA);
            else
              Tr_act =
                  55.262 * m_RevVolume /
                  (Pin_act_terzsum * speedofsound * speedofsound * 4 * rho0) *
                  globalSumRoomB / ((PetscReal)numNodesRoomB);

            m_RevFile << std::setw(15) << lastMiddleTerz;
            m_RevFile << std::setw(15) << Pin_act_terzsum;
            if (getRevRoomNr() == 0)
              m_RevFile << std::setw(15) << globalSumRoomA;
            else
              m_RevFile << std::setw(15) << globalSumRoomB;
            m_RevFile << std::setw(15) << Tr_act;
            m_RevFile << std::endl << std::flush;
          }

          trace("  reverberation time for Terz computed! ");

          Press2ValRoomA.clear();
          Press2ValRoomB.clear();
          Pin_act_terzsum = globalPin_act;

          Press2ValRoomA.resize(m_GlobalNumRowsPerRoom[0]);
          Press2ValRoomB.resize(m_GlobalNumRowsPerRoom[1]);

          computePressuresAtNodesFrequ_Room(solution, MyMesh, Press2ValRoomA,
                                            Press2ValRoomB);
          nr_lastMiddleTerz = 1;
          lastMiddleTerz = currentMiddleTerz;
        }
      }
    }
  }
}

void cOutputBA::computePressureLevelRTC(const PetscReal &CurrentFrequency,
                                        Vec &solution, cMesh &MyMesh) {
  const PetscReal pRef = 20.0e-6;  // reference pressure in air

  PetscReal globalSumRoomA = 0.;
  PetscReal globalSumRoomB = 0.;

  // --- if the file is opened successfully, we'll now compute
  //     the sound pressure levels in each room ...
  trace("    computing local sum ... ");

  computePressures_in_Room(CurrentFrequency, solution, MyMesh, globalSumRoomA,
                           globalSumRoomB);

  // --- output of the computed sound pressure level and the difference
  //     between them
  if (PetscGlobalRank == 0) {
    PetscReal L_A = 0.;  // sound pressure level room A
    PetscReal L_B = 0.;  // sound pressure level room B
    PetscReal Tr = 0.;   // reverberation time
    PetscReal RTC = 0.;  // correction factor

    if (globalSumRoomA > 0.)
      L_A = 10. * std::log10(globalSumRoomA / (pRef * pRef));
    else
      trace("please check : LA <= 0.0");

    if (globalSumRoomB > 0.)
      L_B = 10. * std::log10(globalSumRoomB / (pRef * pRef));
    else
      trace("please check : LB <= 0.0");

    Tr = getInterpolatedValue(CurrentFrequency, RevTimeFrequ, RevTimeVal);
    RTC = 10 * std::log10(m_RevDeviceArea * Tr / (0.163 * m_RevVolume));

    m_RevFile << std::setw(15) << CurrentFrequency;
    m_RevFile << std::setw(15) << L_A;
    m_RevFile << std::setw(15) << L_B;
    m_RevFile << std::setw(15) << L_A - L_B;
    m_RevFile << std::setw(15) << Tr;   // reverberation time
    m_RevFile << std::setw(15) << RTC;  // correction factor
    m_RevFile << std::setw(15)
              << L_A - L_B + RTC;  // corrected transmission loss
    m_RevFile << std::endl << std::flush;
  }

  trace("  sound pressure level with reverberation time correction computed! ");
}

PetscReal cOutputBA::getInterpolatedValue(const PetscReal &CurrentFrequency,
                                          std::vector<PetscReal> &frequencies,
                                          std::vector<PetscReal> &values) {
  PetscReal value = 0.;

  // -- check for lower and upper bound; else interpolate
  if (CurrentFrequency <= frequencies[0]) {
    return values[0];
  } else if (CurrentFrequency >= frequencies[frequencies.size() - 1]) {
    return values[frequencies.size() - 1];
  } else {
    for (int k = 0; k < (int)frequencies.size() - 1; k++) {
      if (frequencies[k] == CurrentFrequency) return values[k];

      if ((frequencies[k] <= CurrentFrequency) &&
          (CurrentFrequency < frequencies[k + 1])) {
        value = (values[k + 1] - values[k]) /
                    (frequencies[k + 1] - frequencies[k]) *
                    (CurrentFrequency - frequencies[k]) +
                values[k];
        return value;
      }
    }

    if (frequencies[(int)frequencies.size() - 1] <= CurrentFrequency)
      return values[frequencies.size() - 1];
  }

  trace(
      "Error in cOutputBA::getInterpolatedValue!!! Please check the "
      "Reverberation File");

  // --- should not get here !
  return value;
}

void cOutputBA::getTerzFrequencies(std::vector<PetscReal> &TerzMitte,
                                   std::vector<PetscReal> &TerzOben,
                                   std::vector<PetscReal> &TerzUnten) {
  int nr = 22;
  const int N = nr + 11;
  TerzMitte.resize(N);  // Werte der Terzmitten-Frequenzen
  TerzOben.resize(N);   // Obergrenze der i-ten Terz
  TerzUnten.resize(N);  // Untergrenze der i-ten Terz

  TerzMitte[0] = 50;
  TerzMitte[1] = 63;
  TerzMitte[2] = 80;
  TerzMitte[3] = 100;
  TerzMitte[4] = 125;
  TerzMitte[5] = 160;
  TerzMitte[6] = 200;
  TerzMitte[7] = 250;
  TerzMitte[8] = 315;
  TerzMitte[9] = 400;
  TerzMitte[10] = 500;

  TerzUnten[0] = 45;
  TerzUnten[1] = 56;
  TerzUnten[2] = 71;
  TerzUnten[3] = 90;
  TerzUnten[4] = 112;
  TerzUnten[5] = 140;
  TerzUnten[6] = 180;
  TerzUnten[7] = 224;
  TerzUnten[8] = 280;
  TerzUnten[9] = 355;
  TerzUnten[10] = 450;

  TerzOben[0] = 56;
  TerzOben[1] = 71;
  TerzOben[2] = 90;
  TerzOben[3] = 112;
  TerzOben[4] = 140;
  TerzOben[5] = 180;
  TerzOben[6] = 224;
  TerzOben[7] = 280;
  TerzOben[8] = 355;
  TerzOben[9] = 450;
  TerzOben[10] = 560;

  for (int k = 0; k < nr; k++) {
    TerzMitte[k + 11] = TerzMitte[k] * 10;
    TerzUnten[k + 11] = TerzUnten[k] * 10;
    TerzOben[k + 11] = TerzOben[k] * 10;
  }
}

void cOutputBA::writeMeshToFile(cMesh &MyMesh) {
  if (wantOutput() > 0) {
    trace("no mesh output possible for BA-Output ...");
  }
}

void cOutputBA::writeMeshToFileGroup(cMesh &MyMesh) {
  if (wantOutput() > 0) {
    trace("no mesh output for groups possible for BA-Output ...");
  }
}

//! Creates a copy of current Object
cOutputfile *cOutputBA::copy() {
  cOutputfile *copy_object = new cOutputBA(*this);

  return copy_object;
}

//! set the filename of the file from where the reverberation time values has to
//! be read
void cOutputBA::setRevTimeFileName(std::string filename) {
  m_RevTimeFileName = filename;

  std::ifstream input;  // inputstream
  PetscReal frequ, Tr;
  std::string dummy;

  message("trying to read '%s' ...\n", m_RevTimeFileName.c_str());
  input.open(m_RevTimeFileName.c_str());

  if (input.fail()) {
    message("  Unable to open data file '%s'\n", m_RevTimeFileName.c_str());
    input.close();
  } else {
    trace("  file opened successfully ...");
  }

  // Read in the header
  std::getline(input, dummy);
  std::getline(input, dummy);
  std::getline(input, dummy);
  std::getline(input, dummy);

  // --------------------------------------------------------------------------
  //   read everything
  // --------------------------------------------------------------------------
  while (input >> frequ) {
    RevTimeFrequ.push_back(frequ);
    input >> dummy;
    input >> dummy;
    input >> Tr;
    RevTimeVal.push_back(Tr);
  }
  input.close();
}

std::ostream &cOutputBA::write(std::ostream &os) const {
  os << " want to compute Build.Acoustics : " << wantOutput() << std::endl;
  if (wantOutput() > 0) {
    os << " mode						  : "
       << getRevMode() << std::endl;
    os << " LoadIDPs					  : "
       << getRevLoadIDPs() << std::endl;
    os << " LoadIDPin					  : "
       << getRevLoadIDPin() << std::endl;
    os << " Volume						  : "
       << getRevVolume() << std::endl;
    os << " DeviceArea          : " << getRevDeviceArea() << std::endl;
    os << " RevFileName         : " << getRevTimeFileName() << std::endl;
    os << " RevRoomNr           : " << getRevRoomNr() << std::endl;
  }

  return os;
}

std::istream &cOutputBA::read(std::istream &is) { return is; }

std::ostream &cOutputBA::writeXml(std::ostream &os) const {
  os << "  <outba><active>" << wantOutput() << "</active>";
  if (wantOutput() > 0) {
    os << "<mode>" << getRevMode() << "</mode>";
    os << "<LoadIDPs>" << getRevLoadIDPs() << "</LoadIDPs>";
    os << "<LoadIDPin>" << getRevLoadIDPin() << "</LoadIDPin>";
    os << "<Volume>" << getRevVolume() << "</Volume>";
    os << "<DeviceArea>" << getRevDeviceArea() << "</DeviceArea>";
    os << "<RevFileName>" << getRevTimeFileName() << "</RevFileName>";
    os << "<RevRoomNr>" << getRevRoomNr() << "</RevRoomNr>";
  }
  os << "</outba>" << std::endl;

  return os;
}

void cOutputBA::computePressuresRTCAtNodesTerz(
    const PetscReal &CurrentFrequency, Vec &solution, cMesh &MyMesh) {
  const PetscReal pRef = 20.0e-6;  // reference pressure in air

  const PetscInt numNodesRoomA = m_GlobalNumRowsPerRoom[0];
  const PetscInt numNodesRoomB = m_GlobalNumRowsPerRoom[1];

  std::vector<PetscReal> TerzMitte;
  std::vector<PetscReal> TerzOben;
  std::vector<PetscReal> TerzUnten;
  getTerzFrequencies(TerzMitte, TerzOben, TerzUnten);

  PetscReal currentMiddleTerz;

  // -------------------------------------------------------------------------
  //   First Averaging the Reverberation time data into thirds,
  //   than averaging over the room
  // -------------------------------------------------------------------------

  for (int k = 0; k < (int)TerzMitte.size(); k++) {
    if ((TerzUnten[k] <= CurrentFrequency) &&
        (CurrentFrequency < TerzOben[k])) {
      if (lastMiddleTerz == 0.0) lastMiddleTerz = TerzMitte[k];

      currentMiddleTerz = TerzMitte[k];

      if (std::abs(currentMiddleTerz - lastMiddleTerz) < cstGeomEps) {
        computePressuresAtNodesFrequ_Room(solution, MyMesh, Press2ValRoomA,
                                          Press2ValRoomB);
        nr_lastMiddleTerz = nr_lastMiddleTerz + 1;
        lastMiddleTerz = currentMiddleTerz;
      } else {
        // --- now we compute the sum over all processes and gather the result
        // on rank 0
        trace("    computing global sum for Terz... ");
        PetscReal pEff2SumRoomA = 0.;
        PetscReal pEff2SumRoomB = 0.;

        for (int n = 0; n < (int)Press2ValRoomA.size(); n++) {
          pEff2SumRoomA = pEff2SumRoomA + Press2ValRoomA[n];
        }
        for (int n = 0; n < (int)Press2ValRoomB.size(); n++) {
          pEff2SumRoomB = pEff2SumRoomB + Press2ValRoomB[n];
        }

        nr_lastMiddleTerz = 0;

        PetscReal globalSumRoomA = 0.;
        PetscReal globalSumRoomB = 0.;

        MPI_Reduce(&pEff2SumRoomA, &globalSumRoomA, 1, MPIU_REAL, MPI_SUM, 0,
                   PETSC_COMM_WORLD);
        MPI_Reduce(&pEff2SumRoomB, &globalSumRoomB, 1, MPIU_REAL, MPI_SUM, 0,
                   PETSC_COMM_WORLD);

        if (PetscGlobalRank == 0) {
          PetscReal L_A = 0.;  // sound pressure level room A
          PetscReal L_B = 0.;  // sound pressure level room B
          PetscReal Tr = 0.;   // reverberation time
          PetscReal RTC = 0.;  // correction factor

          Tr = getInterpolatedValue(lastMiddleTerz, RevTimeFrequ, RevTimeVal);
          RTC = 10 * std::log10(m_RevDeviceArea * Tr / (0.163 * m_RevVolume));

          if (globalSumRoomA > 0.)
            L_A = 10. * std::log10(globalSumRoomA / ((PetscReal)numNodesRoomA) /
                                   (pRef * pRef));

          if (globalSumRoomB > 0.)
            L_B = 10. * std::log10(globalSumRoomB / ((PetscReal)numNodesRoomB) /
                                   (pRef * pRef));

          m_RevFile << std::setw(15) << lastMiddleTerz;
          m_RevFile << std::setw(15) << L_A;
          m_RevFile << std::setw(15) << L_B;
          m_RevFile << std::setw(15) << L_A - L_B;
          m_RevFile << std::setw(15) << Tr;
          m_RevFile << std::setw(15) << RTC;
          m_RevFile << std::setw(15) << L_A - L_B + RTC;
          m_RevFile << std::endl << std::flush;
        }

        trace("  sound pressure level for Terz computed! ");

        Press2ValRoomA.clear();
        Press2ValRoomB.clear();

        Press2ValRoomA.resize(m_GlobalNumRowsPerRoom[0]);
        Press2ValRoomB.resize(m_GlobalNumRowsPerRoom[1]);

        computePressuresAtNodesFrequ_Room(solution, MyMesh, Press2ValRoomA,
                                          Press2ValRoomB);
        nr_lastMiddleTerz = nr_lastMiddleTerz + 1;
        lastMiddleTerz = currentMiddleTerz;
      }
    }
  }
}

void cOutputBA::computePressuresAtNodesFrequ_Room(
    Vec &solution, cMesh &MyMesh, std::vector<PetscReal> &Press2ValRoomA,
    std::vector<PetscReal> &Press2ValRoomB) {
  PetscInt FirstRow = 0;  // first row of local portion of solution vector
  PetscInt LastRow = 0;   // last row of local portion of solution vector

  //  trace("  computing sound pressure level ... ");

  ierr = VecGetOwnershipRange(solution, &FirstRow, &LastRow);
  INFAMCHKERRQ(ierr);

  // --- if the file is opened successfully, we'll now compute
  //     the sound pressure levels in each room ...
  trace("    computing local sum ... ");
  PetscScalar value = 0.;

  PetscReal pEff2;
  PetscInt pNrRoomA = 0;
  PetscInt pNrRoomB = 0;

  short room;
  PetscInt row;
  for (ItMapNodes itN = MyMesh.getFirstNode(); itN != MyMesh.getLastNode();
       itN++) {
    row = itN->second->getGlobalRow(fluid);

    if ((FirstRow <= row) && (row < LastRow)) {
      room = itN->second->getRoom();
      VecGetValues(solution, 1, &row, &value);

#ifdef PETSC_USE_COMPLEX
      pEff2 = (value.real() * value.real() + value.imag() * value.imag()) / 2;
#else
      pEff2 = std::abs(value / M_SQRT2) * std::abs(value / M_SQRT2);
#endif

      if (room == 1) {
        Press2ValRoomA[pNrRoomA] = Press2ValRoomA[pNrRoomA] + pEff2;
        pNrRoomA = pNrRoomA + 1;
      } else if (room == 2) {
        Press2ValRoomB[pNrRoomB] = Press2ValRoomB[pNrRoomB] + pEff2;
        pNrRoomB = pNrRoomB + 1;
      }
    }
  }
}

void cOutputBA::computePressures_in_Room(const PetscReal &CurrentFrequency,
                                         Vec &solution, cMesh &MyMesh,
                                         PetscReal &globalSumRoomA,
                                         PetscReal &globalSumRoomB) {
  PetscInt FirstRow = 0;  // first row of local portion of solution vector
  PetscInt LastRow = 0;   // last row of local portion of solution vector

  //  trace("  computing sound pressure level ... ");

  ierr = VecGetOwnershipRange(solution, &FirstRow, &LastRow);
  INFAMCHKERRQ(ierr);

  // --- if the file is opened successfully, we'll now compute
  //     the sound pressure levels in each room ...
  trace("    computing local sum ... ");
  PetscScalar value = 0.;

  PetscReal pEff2;  // current effective value of the sound pressure
  PetscReal pEff2SumRoomA = 0.;  // sum of effective values in room A
  PetscReal pEff2SumRoomB = 0.;  // sum of effective values in room B

  //  const PetscReal pRef           = 20.0e-6; // reference pressure in air
  const PetscReal numNodesRoomA = (PetscReal)m_GlobalNumRowsPerRoom[0];
  const PetscReal numNodesRoomB = (PetscReal)m_GlobalNumRowsPerRoom[1];

  short room;
  PetscInt row;
  for (ItMapNodes itN = MyMesh.getFirstNode(); itN != MyMesh.getLastNode();
       itN++) {
    row = itN->second->getGlobalRow(fluid);

    if ((FirstRow <= row) && (row < LastRow)) {
      room = itN->second->getRoom();
      VecGetValues(solution, 1, &row, &value);

#ifdef PETSC_USE_COMPLEX
      pEff2 = (value.real() * value.real() + value.imag() * value.imag()) / 2;
#else
      pEff2 = std::abs(value / M_SQRT2) * std::abs(value / M_SQRT2);
#endif

      if (room == 1)
        pEff2SumRoomA += pEff2 / (numNodesRoomA);
      else if (room == 2)
        pEff2SumRoomB += pEff2 / (numNodesRoomB);
    }
  }

  // --- now we compute the sum over all processes and gather the result on rank
  // 0
  trace("    computing global sum ... ");
  globalSumRoomA = 0.;
  globalSumRoomB = 0.;
  MPI_Reduce(&pEff2SumRoomA, &globalSumRoomA, 1, MPIU_REAL, MPI_SUM, 0,
             PETSC_COMM_WORLD);
  MPI_Reduce(&pEff2SumRoomB, &globalSumRoomB, 1, MPIU_REAL, MPI_SUM, 0,
             PETSC_COMM_WORLD);

  trace("  sound pressure level computed! ");
}

void cOutputBA::checkForSides(cMesh &MyMesh) {
  for (ItMapNodes it = MyMesh.getFirstNode(); it != MyMesh.getLastNode();
       it++) {
    if (it->second->getRoom() == 1)
      m_GlobalNumRowsPerRoom[0]++;
    else if (it->second->getRoom() == 2)
      m_GlobalNumRowsPerRoom[1]++;
  }

  message("  total number : %6d\n",
          m_GlobalNumRowsPerRoom[0] + m_GlobalNumRowsPerRoom[1]);
  message("        room A : %6d\n", m_GlobalNumRowsPerRoom[0]);
  message("        room B : %6d\n", m_GlobalNumRowsPerRoom[1]);

  trace(" finished!");
}

void cOutputBA::parseNdsFile(cMesh &MyMesh) {
  std::ifstream rein;  // inputstream
  std::string text;

  short seite = 1;  // side [1,2]

  // -------------------------------------------------------------------------
  //   opening nds file
  //   used to determine the nodes belonging to the two adjacent rooms
  // -------------------------------------------------------------------------
  text = MyMesh.getFilename();
  if (text.rfind(cstInputFileSuffix) ==
      (text.length() - cstInputFileSuffix.length()))
    text.erase(text.length() - cstInputFileSuffix.length());
  text += ".nds";

  rein.open(text.c_str());
  message("trying to read '%s' ...\n", text.c_str());

  if (rein.fail()) {
    throw cException("not able to open " + text, __FILE__, __LINE__);
    ExitApp();
  } else {
    trace("  file opened successfully ...");
  }

  // --------------------------------------------------------------------------
  //   read everything
  // --------------------------------------------------------------------------
  while (rein >> text) {
    // ------------------------------------------------------------------------
    //   look which room we're processing (1 or 2)
    // ------------------------------------------------------------------------
    if (text == "Seite") {
      rein >> seite;
      if ((seite != 1) && (seite != 2)) {
        message("invalid value for 'seite' : %d\n", seite);
        throw cException("cMesh::parseNds() : invalid value for 'Seite'",
                         __FILE__, __LINE__);
      }
    }
    // ------------------------------------------------------------------------
    //   ... nodes we have to evaluate ?
    // ------------------------------------------------------------------------
    else {
      // ------------------------------------------------------------------------
      //   looking for ':'
      // ------------------------------------------------------------------------
      int pos1 = (int)text.find(":");
      int pos2 = (int)text.rfind(":");

      // ------------------------------------------------------------------------
      //   ':' found ?
      // ------------------------------------------------------------------------
      if (pos1 >= 0) {
        // ----------------------------------------------------------------------
        //   there is one ':' (pos1 == pos2) or are there two ':'
        // ----------------------------------------------------------------------
        if (pos1 == pos2) {
          std::string oben(text);
          std::string unten(text);

          oben.erase(0, pos1 + 1);
          unten.erase(pos1);

          for (int i = std::atoi(unten.c_str()); i <= std::atoi(oben.c_str());
               i++)
            MyMesh.getNode(i)->setRoom(seite);
        }
        // ----------------------------------------------------------------------
        //   more than one ':' (pos1 != pos2)
        // ----------------------------------------------------------------------
        else {
          std::string oben(text);
          std::string unten(text);
          std::string step(text);

          oben.erase(0, pos1 + 1);
          oben.erase(oben.find(":"));
          unten.erase(pos1);
          step.erase(0, pos2 + 1);

          for (int i = std::atoi(unten.c_str()); i <= std::atoi(oben.c_str());
               i += std::atoi(step.c_str()))
            MyMesh.getNode(i)->setRoom(seite);
        }
      }
      // ------------------------------------------------------------------------
      //   no ':' found
      // ------------------------------------------------------------------------
      else {
        MyMesh.getNode(std::atoi(text.c_str()))->setRoom(seite);
      }
    }
  }

  // --- close nds file
  rein.close();

  /*
  // -------------------------------------------------------------------------
  //   tell the fluid-structure coupling elements to which room they belong
  //   This information will be used later on to compute the sound intensity
  //   on each room's wall
  // -------------------------------------------------------------------------
  MyMesh.associateInterfaceElementsToRoom();
  */
  trace("finished.");
}

void cOutputBA::computePressuresAtNodesTerz(const PetscReal &CurrentFrequency,
                                            Vec &solution, cMesh &MyMesh) {
  const PetscReal pRef = 20.0e-6;  // reference pressure in air

  const PetscInt numNodesRoomA = m_GlobalNumRowsPerRoom[0];
  const PetscInt numNodesRoomB = m_GlobalNumRowsPerRoom[1];

  std::vector<PetscReal> TerzMitte;
  std::vector<PetscReal> TerzOben;
  std::vector<PetscReal> TerzUnten;
  getTerzFrequencies(TerzMitte, TerzOben, TerzUnten);

  PetscReal currentMiddleTerz;

  // -------------------------------------------------------------------------
  //   First Averaging the Reverberation time data into thirds,
  //   than averaging over the room
  // -------------------------------------------------------------------------

  for (int k = 0; k < (int)TerzMitte.size(); k++) {
    if ((TerzUnten[k] <= CurrentFrequency) &&
        (CurrentFrequency < TerzOben[k])) {
      if (lastMiddleTerz == 0.0) lastMiddleTerz = TerzMitte[k];

      currentMiddleTerz = TerzMitte[k];

      if (std::abs(currentMiddleTerz - lastMiddleTerz) < cstGeomEps) {
        computePressuresAtNodesFrequ_Room(solution, MyMesh, Press2ValRoomA,
                                          Press2ValRoomB);
        nr_lastMiddleTerz = nr_lastMiddleTerz + 1;
        lastMiddleTerz = currentMiddleTerz;
      } else {
        // --- now we compute the sum over all processes and gather the result
        // on rank 0
        trace("    computing global sum for Terz... ");
        PetscReal pEff2SumRoomA = 0.;
        PetscReal pEff2SumRoomB = 0.;

        for (int n = 0; n < (int)Press2ValRoomA.size(); n++) {
          pEff2SumRoomA = pEff2SumRoomA + Press2ValRoomA[n];
        }
        for (int n = 0; n < (int)Press2ValRoomB.size(); n++) {
          pEff2SumRoomB = pEff2SumRoomB + Press2ValRoomB[n];
        }

        nr_lastMiddleTerz = 0;

        PetscReal globalSumRoomA = 0.;
        PetscReal globalSumRoomB = 0.;

        MPI_Reduce(&pEff2SumRoomA, &globalSumRoomA, 1, MPIU_REAL, MPI_SUM, 0,
                   PETSC_COMM_WORLD);
        MPI_Reduce(&pEff2SumRoomB, &globalSumRoomB, 1, MPIU_REAL, MPI_SUM, 0,
                   PETSC_COMM_WORLD);

        if (PetscGlobalRank == 0) {
          PetscReal L_A = 0.;  // sound pressure level room A
          PetscReal L_B = 0.;  // sound pressure level room B

          if (globalSumRoomA > 0.)
            L_A = 10. * std::log10(globalSumRoomA / ((PetscReal)numNodesRoomA) /
                                   (pRef * pRef));

          if (globalSumRoomB > 0.)
            L_B = 10. * std::log10(globalSumRoomB / ((PetscReal)numNodesRoomB) /
                                   (pRef * pRef));

          m_RevFile << std::setw(15) << lastMiddleTerz;
          m_RevFile << std::setw(15) << L_A;
          m_RevFile << std::setw(15) << L_B;
          m_RevFile << std::setw(15) << L_A - L_B;
          m_RevFile << std::endl << std::flush;
        }

        trace("  sound pressure level for Terz computed! ");

        Press2ValRoomA.clear();
        Press2ValRoomB.clear();

        Press2ValRoomA.resize(m_GlobalNumRowsPerRoom[0]);
        Press2ValRoomB.resize(m_GlobalNumRowsPerRoom[1]);

        computePressuresAtNodesFrequ_Room(solution, MyMesh, Press2ValRoomA,
                                          Press2ValRoomB);
        nr_lastMiddleTerz = nr_lastMiddleTerz + 1;
        lastMiddleTerz = currentMiddleTerz;
      }
    }
  }
}
