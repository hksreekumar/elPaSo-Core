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

#ifndef INFAM_OUTPUTVTK_H
#define INFAM_OUTPUTVTK_H

#include "outputfile.h"

//! @brief accessing the results
//! @author Meike Wulkau
//! @date 27.10.2010
//!
//! cOutputVTK writes the results of the computation to an vtk outputfile.
class cOutputVTK : public cOutputfile {
 private:
  std::ofstream VtkOutput;      ///< vtk file (vtk)
  std::ofstream nodeConOutput;  ///< node connectivity file

  //! exports mesh to vtk in ASCII
  //! @param myMesh mesh to export
  //! @param byGroup mesh is ordered by groups defined in Patran
  void exportMeshToVtk(cMesh *myMesh);

  //! exports the solution to a file in vkt format (ASCII)
  //! @param myMesh underlying mesh
  //! @param NumStep sequential step number
  //! @param Step frequency/time step where NumStep corresponds to
  void exportSolutionToVtk(cMesh *myMesh, const PetscInt &NumStep,
                           const PetscReal &Step);

  //! export to vtk file
  //! @param MyMesh underlying mesh
  //! @param output outputstream (vtk-file)
  void exportToVtk(cMesh &MyMesh, std::ostream &output,
                   std::ostream &nodeConOutput = std::cout,
                   const PetscReal &Step = 0.0);

  void nodesNeeded(cElementFEM *ptrElement, int &NodesNeeded,
                   int &VtkElementType);

  void collectOutputStrings(std::ostream &output, std::stringstream &string,
                            std::string tag1, std::string tag2);

 public:
  cOutputVTK();
  cOutputVTK(cOutputVTK &other);
  ~cOutputVTK();

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
