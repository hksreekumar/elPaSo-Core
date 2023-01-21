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

#ifndef INFAM_OUTPUTFILE_H
#define INFAM_OUTPUTFILE_H

#include "../../fedata/mesh.h"
#include "../log/logging.h"
#include "../mytypes.h"
#include "gzstream.h"

//! @brief accessing the results
//! @author Meike Wulkau
//! @date 27.10.2010
//!
//! cOutputfile is the base class of all file output objects.
class cOutputfile : public virtual cLogging {
 private:
  int m_WantOutput;  ///< tells us if an output is requested and if yes, what
                     ///< step period is defined

  bool modify_FilenameRoot;
  std::string FilenameRoot_add;

 protected:
  bool m_InitPostProcess;          // becomes true if Output is initialized
  std::vector<bool> m_DispVector;  // entity becomes true if dof exists
  int m_numberofDofs;              // number of unknowns
  int m_dispX, m_dispW;            // rows os disp_x and disp_w vector

  Vec *m_OutSolution;
  Vec *m_OutVelocity;
  Vec *m_OutAcceleration;
  Vec *m_OutStressesCells;
  Vec *m_OutStressesCellsSec2;
  Vec *m_OutStressesNodes;

  //! Prescribed displacements will be checked and stored in displacement vector
  void initPostProcess();

  //! create the name of inputfile without any suffix
  //! @param FileNameRoot name of inputfile without any suffix
  //! @param MyMesh mesh definition
  void createFilenameRoot(std::string &FilenameRoot, cMesh *myMesh);

 public:
  cOutputfile();
  cOutputfile(cOutputfile &other);
  virtual ~cOutputfile();

  void setOutSolution(Vec *vec) { m_OutSolution = vec; }
  Vec *getOutSolution() { return m_OutSolution; }
  void setOutVelocity(Vec *vec) { m_OutVelocity = vec; }
  Vec *getOutVelocity() { return m_OutVelocity; }
  void setOutAcceleration(Vec *vec) { m_OutAcceleration = vec; }
  Vec *getOutAcceleration() { return m_OutAcceleration; }
  void setOutStressesCells(Vec *vec) { m_OutStressesCells = vec; }
  Vec *getOutStressesCells() { return m_OutStressesCells; }
  void setOutStressesCellsSec2(Vec *vec) { m_OutStressesCellsSec2 = vec; }
  Vec *getOutStressesCellsSec2() { return m_OutStressesCellsSec2; }
  void setOutStressesNodes(Vec *vec) { m_OutStressesNodes = vec; }
  Vec *getOutStressesNodes() { return m_OutStressesNodes; }

  //! check, if an output is requested and if yes, what step period is defined
  inline int wantOutput(void) const { return m_WantOutput; }

  //! tell the program if an output is requested and if yes, what step period is
  //! defined
  void setWantOutput(int flag) { m_WantOutput = flag; }

  inline bool getModifyFilenameRoot(void) const { return modify_FilenameRoot; }

  void setModifyFilenameRoot(bool flag) { modify_FilenameRoot = flag; }

  inline std::string getFilenameRootadd(void) const { return FilenameRoot_add; }

  void setFilenameRootadd(std::string addstr) { FilenameRoot_add = addstr; }

  //! write the results of a single computation step to a file.
  //! @param MyMesh given mesh
  //! @param NrStep number of the current computational step
  //! @param Step frequency or time step
  virtual void writeResultsToFile(cMesh &MyMesh, const int &NrStep,
                                  const double &Step) = 0;

  //! write the Mesh to a file in the different formats.
  //! @param MyMesh given mesh
  virtual void writeMeshToFile(cMesh &MyMesh) = 0;

  //! write the Mesh to a file in the different formats separeted by groups.
  //! @param MyMesh given mesh
  virtual void writeMeshToFileGroup(cMesh &MyMesh) = 0;

  //! Creates a copy of current Object!
  virtual cOutputfile *copy() = 0;

  //! reads the mesh information of the .res file
  static void readResFile(const std::string &FileNameRoot, cMesh &MyMesh,
                          PetscInt &Rows, PetscInt &Steps);

  //! reads the solution of the .stp file
  static void readStpGzFile(const std::string &FileNameRoot, PetscInt Rows,
                            int step, PetscReal &StepVal, Vec &solution);

  //! write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &write(std::ostream &os) const = 0;

  //! reads a single object of a stream
  //! @param is inputstream
  //! @return modified outputstream
  virtual std::istream &read(std::istream &is) = 0;

  //! write this object in XML
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &writeXml(std::ostream &os) const = 0;
};

#endif
