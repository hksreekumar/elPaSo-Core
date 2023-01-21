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

#ifndef INFAM_OUTPUTLEV_H
#define INFAM_OUTPUTLEV_H

#include "outputfile.h"

//! @brief accessing the results
//! @author Meike Wulkau
//! @date 27.10.2010
//!
//! cOutputLEV writes the results of the computation to an lev outputfile.
class cOutputLEV : public cOutputfile {
 private:
  PetscInt m_GlobalNumRowsPerRoom[2];  ///< global number of nodes located
                                       ///< within sending/receiving room

  std::string
      m_LevelFileName;  ///< name of the file with computed sound pressures
  std::ofstream
      m_LevelFile;  ///< file to which computed sound pressures are written

  //! reads a file with suffix .nds that tells us which nodes
  //! belongs to room A and room B.
  //! @param Filename
  void checkForSides(cMesh &MyMesh);

  //! parse the file that defines which node is used either for
  //! the sending room or the receiving room during the computation
  //! of the sound pressure level
  void parseNdsFile(cMesh &MyMesh);

 protected:
  //! computes the soundpressurelevel of the two rooms described by
  //! RowsRoomA and RowsRoomB.
  //! @param CurrentFrequency current frequency - used as x-axis
  //! @param MyMesh used to prepare the output
  void computePressureLevel(const PetscReal &CurrentFrequency, cMesh &MyMesh);

 public:
  cOutputLEV();
  cOutputLEV(cOutputLEV &other);
  ~cOutputLEV();

  //! write the results of a single computation step to a file.
  //! @param MyMesh given mesh
  //! @param NrStep number of the current computational step
  //! @param Step frequency or time step
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
