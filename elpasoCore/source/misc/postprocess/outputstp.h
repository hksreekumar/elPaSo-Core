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

#ifndef INFAM_OUTPUTSTP_H
#define INFAM_OUTPUTSTP_H

#include "../filter.h"
#include "outputfile.h"

//! @brief accessing the results
//! @author Meike Wulkau
//! @date 27.10.2010
//!
//! cOutputAVS writes the results of the computation to an avs outputfile.
class cOutputSTP : public cOutputfile, public cFilter {
 private:
  //! write the results of a single computation step to a file.
  //! For each compuation step a new file is created.
  //! The values of each degree of freedom are written to a
  //! auto-compressing stream. To view this file use zmore.
  //! @param Step number of the current computational step
  //! @param Frequenz frequency or time step
  //! @param Filename name of outputfile (scheme: <Filename>.<Step>.stp.gz)
  void writeResultSTP(const int &step, const double &frequenz,
                      const std::string &Filename);

 public:
  cOutputSTP();
  cOutputSTP(cOutputSTP &other);
  ~cOutputSTP();

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
