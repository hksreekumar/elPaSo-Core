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

#ifndef INFAM_OUTPUT_H
#define INFAM_OUTPUT_H

#include "../../fedata/dof.h"
#include "../../fedata/mesh.h"
#include "../filter.h"
#include "../log/logging.h"
#include "../mytypes.h"
#include "gzstream.h"
#include "outputfile.h"

//! @brief accessing the results
//! @author Dirk Clasen
//! @date 30.05.2005
//! cOutput writes the results of the computation to an outputfile.
class cOutput : public virtual cLogging {
 private:
  std::string
      m_BaseName;  ///< prefix of the outputfile: <m_BaseName>.<step>.stp.gz

  std::list<cOutputfile *> m_OutToFile;

  bool m_CreateAK3;  ///< tells us if the users wants a ak3 file of current mesh

  bool m_LogEnergy;    ///< tells us if the users wants an energy-output-file of
                       ///< each step's solution
  bool m_LogVelocity;  ///< tells us if the users wants an velocity-output-file
                       ///< of each step's solution

  //! check, if we have to create an info-file that contains the
  //! mesh information
  bool wantResFile(void);

  bool m_elementsAtNodesCounted;

 protected:
  void countElementsAtNodes(cMesh &MyMesh);

 public:
  cOutput();
  ~cOutput();

  void setOutSolution(Vec *vec);
  void setOutVelocity(Vec *vec);
  void setOutAcceleration(Vec *vec);
  void setOutStressesCells(Vec *vec);
  void setOutStressesCellsSec2(Vec *vec);
  void setOutStressesNodes(Vec *vec);

  //! inserst a new pointer to an output-to-file-object
  void insertOutput(cOutputfile *m_out);

  //! iterator pointing to first output-to-file-object
  std::list<cOutputfile *>::iterator getFirstOutputElement(void) {
    return m_OutToFile.begin();
  }

  //! iterator pointing to last output-to-file-object
  std::list<cOutputfile *>::iterator getLastOutputElement(void) {
    return m_OutToFile.end();
  }

  /**
   * Generates a XML-file which contains the current mesh
   */
  void createXmlOutputOfCurrentMesh(cMesh *myMesh, std::ostream &myAnalysis,
                                    std::ostream &myOutput, PetscInt step);

  /**
   * create the name of inputfile without any suffix
   * @param FileNameRoot name of inputfile without any suffix
   * @param MyMesh mesh definition
   */
  static void createFilenameRoot(std::string &FilenameRoot, cMesh *myMesh);

  //! returns the prefix of the outputfiles
  inline std::string getBaseName(void) const { return m_BaseName; }

  //! sets a new prefix for the outputfiles
  inline void setBaseName(const std::string &BaseName) {
    m_BaseName = BaseName;
  }

  //! check, if we have to create a ak3-file that contains the
  //! current mesh
  inline bool createAK3(void) const { return m_CreateAK3; }

  //! check, if we have to create a energy-output-file that contains the
  //! solution of each step
  inline bool logEnergy(void) const { return m_LogEnergy; }

  //! check, if we have to create a velocity-output-file that contains the
  //! solution of each step
  inline bool logVelocity(void) const { return m_LogVelocity; }

  //! tell the program if we want to get vtk files of the result
  void setCreateAK3(bool flag) { m_CreateAK3 = flag; }

  //! tell the program if we want to get vtk files of the result
  void setLogEnergy(bool flag) { m_LogEnergy = flag; }

  //! tell the program if we want to get the velocity written to file
  void setLogVelocity(bool flag) { m_LogVelocity = flag; }

  void writeResults(cMesh &MyMesh, const int &NrStep, const double &Step);

  //! write the Mesh to a file in the different formats.
  //! @param MyMesh given mesh
  void writeMesh(cMesh &MyMesh);

  //! write the Mesh to a file in the different formats separeted by groups.
  //! @param MyMesh given mesh
  void writeMeshGroup(cMesh &MyMesh);

  //! write the solution to a file in the different formats.
  //! @param MyMesh given mesh
  //! @param step step to write
  void writeSolutionfromStepFile(cMesh &MyMesh, PetscInt &step);

  //! write all solutions to a file in the different formats.
  //! @param MyMesh given mesh
  void writeAllSolutionsfromStepFiles(cMesh &MyMesh, PetscInt steps);

  //! output of eigenvalues
  //! @param RealParts real part of eigenvalues
  //! @param ImagParts imaginary part of eigenvalues
  //! @param Mesh discretisation used to compute
  void writeResults(std::list<PetscReal> &RealParts,
                    std::list<PetscReal> &ImagParts, cMesh &Mesh);
  void writeResults(std::list<PetscReal> &RealParts,
                    std::list<PetscReal> &ImagParts,
                    std::list<PetscReal> &RealParts2,
                    std::list<PetscReal> &ImagParts2, cMesh &Mesh);

  //! write .res file
  //! @param Mesh used discretisation
  //! @param Steps number of frequency/time steps computed
  //! @param NumberOfUnknowns number of degrees of freedom in system matrix
  void writeInfoFile(cMesh &Mesh, const int &Steps,
                     const int &NumberOfUnknowns);

  //! write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &write(std::ostream &os);

  //! reads a single object of a stream
  //! @param is inputstream
  //! @return modified outputstream
  std::istream &read(std::istream &is);

  //! write this object in XML
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &writeXml(std::ostream &os);
};

//! overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os, cOutput &other) {
  return other.write(os);
}

//! overloaded inputoperator
inline std::istream &operator>>(std::istream &is, cOutput &other) {
  return other.read(is);
}

#endif
