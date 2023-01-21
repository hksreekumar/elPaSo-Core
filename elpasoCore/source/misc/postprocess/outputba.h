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

#ifndef INFAM_OUTPUTBA_H
#define INFAM_OUTPUTBA_H

#include "outputfile.h"

//! @brief accessing the results
//! @author Meike Wulkau
//! @date 27.10.2010
//!
//! cOutputBA writes the results of the computation to an vtk outputfile.
class cOutputBA : public cOutputfile {
 private:
  std::string m_RevFileName;  ///< name of the file for the results of the
                              ///< building acoustics calculation
  std::ofstream m_RevFile;  ///< file to which computed building acoustic values
                            ///< are written

  PetscInt m_RevLoadIDPs;   ///<  number of the first  Load ID Pabs
  PetscInt m_RevLoadIDPin;  ///<  number of the second Load ID Pin
  PetscInt m_RevMode;       ///<  Mode of the building acoustics calculation
  ///<  1: reverberation time from surface sound power 2: Transmission loss with
  ///<  sound power through sample
  ///   3: energy density and intensity vectors in all fluid points
  ///   4: phase shift of pressure and velocity 5: reverberation time with mean
  ///   sound pressure 6: Transmission loss with mean sound power level and
  ///   reverberation time correction
  PetscReal m_RevVolume;      ///<  Volume of the reverberation room
  PetscReal m_RevDeviceArea;  ///<  Area of the test device
  PetscInt m_RevRoomNr;       ///<  Nr of the room for the reverberation time

  std::string m_RevTimeFileName;  ///<  name of the file of the precalculated
                                  ///<  reverberation time of room 2
  std::vector<PetscReal>
      RevTimeFrequ;  ///< vector where the reverberation time frequencies can be
                     ///< read in from file
  std::vector<PetscReal> RevTimeVal;  ///< vector where the reverberation time
                                      ///< values can be read in from file

  PetscInt m_GlobalNumRowsPerRoom[2];  ///< global number of nodes located
                                       ///< within sending/receiving room

  std::vector<PetscReal> Press2ValRoomA;
  std::vector<PetscReal> Press2ValRoomB;
  PetscReal Pin_act_terzsum;
  PetscReal lastMiddleTerz;
  PetscInt nr_lastMiddleTerz;

  //  std::string   m_LevelFileName;  ///< name of the file with computed sound
  //  pressures std::ofstream m_LevelFile;      ///< file to which computed
  //  sound pressures are written

  //! reads a file with suffix .nds that tells us which nodes
  //! belongs to room A and room B.
  //! @param Filename
  void checkForSides(cMesh &MyMesh);

  //! parse the file that defines which node is used either for
  //! the sending room or the receiving room during the computation
  //! of the sound pressure level
  void parseNdsFile(cMesh &MyMesh);

  PetscReal getInterpolatedValue(const PetscReal &CurrentFrequency,
                                 std::vector<PetscReal> &frequencies,
                                 std::vector<PetscReal> &values);
  void getTerzFrequencies(std::vector<PetscReal> &TerzMitte,
                          std::vector<PetscReal> &TerzOben,
                          std::vector<PetscReal> &TerzUnten);

 protected:
  //! computes the reverberation time of one room
  //! @param CurrentFrequency current frequency - used as x-axis
  //! @param solution the current solution vector of frequency response analysis
  //! @param MyMesh used to prepare the output
  void computeBuildingAcoustics(cMesh &MyMesh, const int &NrStep,
                                const double CurrentFrequency, Vec &solution);

  void writeHeader(cMesh &MyMesh);

  void writeSoundPowerStep(cMesh &MyMesh, const double CurrentFrequency,
                           Vec &solution);

  void writeEDandIntNodes(cMesh &MyMesh, const int &NrStep,
                          const double CurrentFrequency, Vec &solution);

  void writePhaseShiftNodes(cMesh &MyMesh, const int &NrStep,
                            const double CurrentFrequency, Vec &solution);

  void computePressures_in_Room(const PetscReal &CurrentFrequency,
                                Vec &solution, cMesh &MyMesh,
                                PetscReal &pEffSumRoomA,
                                PetscReal &pEffSumRoomB);

  void computePressuresAtNodesTerz(const PetscReal &CurrentFrequency,
                                   Vec &solution, cMesh &MyMesh);
  void computePressuresRTCAtNodesTerz(const PetscReal &CurrentFrequency,
                                      Vec &solution, cMesh &MyMesh);
  void computePressuresAtNodesFrequ_Room(
      Vec &solution, cMesh &MyMesh, std::vector<PetscReal> &Press2ValRoomA,
      std::vector<PetscReal> &Press2ValRoomB);

  void computePressureLevelRTC(const PetscReal &CurrentFrequency, Vec &solution,
                               cMesh &MyMesh);

 public:
  cOutputBA();
  cOutputBA(cOutputBA &other);
  ~cOutputBA();

  //! set the number of the first Load ID Pabs
  void setRevLoadIDPs(PetscInt id) { m_RevLoadIDPs = id; }

  //! get the number of the first Load ID Pabs
  inline PetscInt getRevLoadIDPs(void) const { return m_RevLoadIDPs; }

  //! set the number of the second Load ID Pin
  void setRevLoadIDPin(PetscInt id) { m_RevLoadIDPin = id; }

  //! get the number of the second Load ID Pin
  inline PetscInt getRevLoadIDPin(void) const { return m_RevLoadIDPin; }

  //! set the Volume of the reverberation room
  void setRevVolume(PetscReal volume) { m_RevVolume = volume; }

  //! get the Volume of the reverberation room
  inline PetscReal getRevVolume(void) const { return m_RevVolume; }

  //! set the area of the test device
  void setRevDeviceArea(PetscReal area) { m_RevDeviceArea = area; }

  //! get the area of the test device
  inline PetscReal getRevDeviceArea(void) const { return m_RevDeviceArea; }

  //! set the filename of the file from where the reverberation time values has
  //! to be read
  void setRevTimeFileName(std::string filename);

  //! get the filename of the file from where the reverberation time values has
  //! to be read
  inline std::string getRevTimeFileName(void) const {
    return m_RevTimeFileName;
  }

  //! set the room number where the reverberation time values have to be
  //! calculated
  void setRevRoomNr(PetscInt nr) { m_RevRoomNr = nr; }

  //! get the room number where the reverberation time values have to be
  //! calculated
  inline PetscInt getRevRoomNr(void) const { return m_RevRoomNr; }

  //! set the Mode for the building acoustics calculation
  void setRevMode(PetscInt mode) { m_RevMode = mode; }

  //! get the Mode for the building acoustics calculation
  inline PetscInt getRevMode(void) const { return m_RevMode; }

  //! write the results of a single computation step to a file.
  //! @param MyMesh given mesh
  //! @param NrStep number of the current computational step
  //! @param Step frequency or time step
  //! @param Solution solution vector
  //! @param stressesCells solution stress vector of the cells
  //! @param stressesNodes solution stress vector of the nodes
  void writeResultsToFile(cMesh &MyMesh, const int &NrStep, const double &Step);

  //! write the Mesh to a file.
  //! @param MyMesh given mesh
  void writeMeshToFile(cMesh &MyMesh);

  //! write the Mesh to a file separeted by groups.
  //! @param MyMesh given mesh
  void writeMeshToFileGroup(cMesh &MyMesh);

  //! Creates a copy of current Object
  cOutputfile *copy(void);

  //! write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &write(std::ostream &os) const;

  //! reads a single object of a stream
  //! @param is inputstream
  //! @return modified outputstream
  std::istream &read(std::istream &is);

  //! write this object in XML
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &writeXml(std::ostream &os) const;
};

#endif
