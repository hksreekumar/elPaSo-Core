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

#ifndef INFAM_ELEMENTINTERFACE_H
#define INFAM_ELEMENTINTERFACE_H

#include "../../material/fluid/materialfluid.h"
#include "../element.h"

//! @brief simple structure that holds all information used to setup an
//! interfaceelement
//! @author Dirk Clasen
//! @date 19.12.2005
//! The program wrote trash to the files after searching for a large number of
//! interfaceelements on the cluster and opening the new inputfile several
//! times. Therefore, i decided to send the interfaceelements by making use of
//! MPI's user defined datatypes. This seems to be much more robust than re-open
//! the files. The MPI datatype is constructed in the constructor of
//! cPreprocessor.
typedef struct {
  PetscInt Id;        ///< Id of Element
  PetscInt N;         ///< number of nodes
  PetscInt MatFluId;  ///< Id of fluid's material
  PetscInt MatStrId;  ///< Id of structure's material
  PetscInt order;  ///< order of interpolation in thickness direction (for poro
                   ///< elements)
  PetscInt typ;    ///< type of structural element
  PetscReal ori;   ///< orientation of normal vector
  PetscInt NS[9];  ///< nodeids of plate-/shellelement
  PetscInt NF[9];  ///< nodeids of fluidelement
  PetscInt ElFluId;  ///< Id of fluid element
  PetscInt ElStrId;  ///< Id of structural element
} sElementInterface;

inline std::ostream &operator<<(std::ostream &os,
                                const sElementInterface &other) {
  std::string strTyp;

  if (other.typ == KirchhoffPlate)
    strTyp = "Kirch";
  else if ((other.typ == MindlinPlate) || (other.typ == PlaneShell) ||
           (other.typ == PlaneShellDrilling))
    strTyp = "Mindlin";
  else if (other.typ == TimoshenkoBeam)
    strTyp = "Beam";
  else if (other.typ == PoroTrash)
    strTyp = "Poro";
  else if (other.typ == MindlinPoro3dUP)
    strTyp = "Poro3dUP";
  else if (other.typ == MindlinEquiporo)
    strTyp = "MindlinEquiporo";
  else if (other.typ == SwebemBnd)
    strTyp = "Swebem";
  else if (other.typ == FFPoro3dUP)
    strTyp = "FFPoro3dUP";
  else if (other.typ == FFShell)
    strTyp = "FFShell";
  else if (other.typ == PlaneShellPoro2dUP)
    strTyp = "PlaneShellPoro2dUP";
  else {
    PetscPrintf(PETSC_COMM_SELF,
                "UNABLE TO EXPORT THIS TYPE OF INTERFACE ELEMENT (%d)\n",
                other.typ);
    ExitApp();
  }

  os << "<Interface" << strTyp << " N=\"" << other.N << "\"";

  if ((other.typ == PoroTrash) || (other.typ == MindlinPoro3dUP) ||
      (other.typ == MindlinEquiporo) || (other.typ == PlaneShellPoro2dUP))
    os << " mats=\"" << other.MatStrId << "\"";

  os << ">" << std::endl;
  os << "<Id>" << other.Id << "</Id>";

  if (other.typ != SwebemBnd && other.typ != FFShell)
    os << "<MatF>" << other.MatFluId << "</MatF>";

  os << "<ori>" << std::showpoint << other.ori << "</ori>" << std::endl;
  for (int k = 0; k < other.N; k++) os << "<NS>" << other.NS[k] << "</NS>";
  os << std::endl;
  for (int k = 0; k < other.N; k++) os << "<NF>" << other.NF[k] << "</NF>";
  if (other.typ == FFPoro3dUP)
    os << "<StrId>" << other.ElStrId << "</StrId><FluId>" << other.ElFluId
       << "</FluId>";
  os << std::endl;
  os << "</Interface" << strTyp << ">";

  return os;
}

//! @brief base class for coupling acoustic fluid domains and structural domains
//! @author Dirk Clasen
//! @date 16.08.2006
//! m_Nodes[] holds the nodes of the fluid domain, the nodes of the structural
//! domain that correspond to the fluid domain nodes are stored in a different
//! vector.
//! @todo incorporating non-homogenous boundary conditions
class cElementInterface : public cElement, private tCounter<cElementInterface> {
 private:
  short m_Room;  ///< room number element belongs to

  /**
   * orientation of normalvectors of both domains
   *  \f$ n_f \f$ : normalvector of fluid's surface
   *  \f$ n_s \f$ : normalvector of structural domain
   *  +1 : \f$ n_f = -n_s \f$
   *  -1 : \f$ n_f =  n_s \f$
   */
  PetscReal m_Orientation;

 protected:
  //! @brief this function performs the computation of C and C^T, multiplies
  //! both by the factors C and C^T and inserts them into the system matrices A
  //! and B by calling insertMatrices(). A and B are both the stiffnessmatrix
  //! for frequency domain analysis. In time domain, A is the stiffness matrix
  //! and B the massmatrix, respectively.
  virtual void addContribution(Mat &A, Mat &B, const PetscScalar &factorC,
                               const PetscScalar &factorCT) = 0;

  //! @brief apply the nodal boundary conditions to the coupling terms
  virtual void applyBoundaryConditions(cElementMatrix &C,
                                       cElementMatrix &CT) = 0;

  //! @brief insert the coupling matrices to global system matrices. C is added
  //! to A and C^T is added to B. In frequency domain analysis, A and B are the
  //! same matrix whereas in time domain A is the stiffnessmatrix and B the
  //! massmatrix, respectively.
  virtual void insertMatrices(Mat &A, Mat &B, cElementMatrix &C,
                              cElementMatrix &CT) = 0;

  //! @brief node ids of structural element that
  //! are attached to the fluid nodes.
  //! The locations of m_Nodes[i] and m_MatchingNodes[i] are
  //! the same.
  std::vector<cNode *> m_MatchingNodes;

  int StrId;  //! Id of structural element
  int FluId;  //! Id of fluid element

  //! @brief compute the coupling matrix
  virtual void computeCouplingMatrix(cElementMatrix &C) = 0;

  //! @brief return the density of the fluid. This is needed by
  //! addContributionFrequencyDomain() and casting to fluid material
  //! does not work for coupling of Mindlin plate and 3d poroelastic
  //! elements
  virtual PetscScalar getRhoF(void) const = 0;

  //! @brief returns the element type
  virtual eElementType getElementType(void) const = 0;

 public:
  //! @brief Constructor
  cElementInterface(short NumberOfNodes);
  //! @brief Copy constructor
  cElementInterface(const cElementInterface &other);
  //! @brief Destructor
  ~cElementInterface();

  //! @brief denotes the shape of this element
  virtual eElementShape getElementShape(void) const = 0;

  //! @brief returns the physical meaning of this element type
  ePhysicsType getPhysicsType(void) const { return Undefined; }

  //! @brief returns the number of instances of this objecttype
  static size_t howMany(void) { return tCounter<cElementInterface>::howMany(); }

  //! @brief assing a new room id to this element
  void setRoomId(const short &id) { m_Room = id; }

  //! @brief return the number of the room this element belongs to
  short getRoomId(void) const { return m_Room; }

  //! @brief return the orientation of the two normals
  PetscReal getOrientation(void) const { return m_Orientation; }

  //! @brief set a new orientation of the two normals
  void setOrientation(const PetscReal &newori) { m_Orientation = newori; }

  //! @brief return a vector that tells us which structural domains of
  //! freedom are coupled to the acoustic/porous domain.
  virtual std::vector<eKnownDofs> getDofsStructure(void) const = 0;

  //! @brief return a coupling node
  //! @param Index index in incidence table
  //! @return pointer to a node
  cNode *getMatchingNode(int Index) const;

  //! @brief insert a coupling node
  //! @param Index index in incidence table
  //! @param ptrNode pointer to a node
  void setMatchingNode(int Index, cNode *ptrNode);

  //! @brief assign the material parameters to one of the coupled domains
  virtual void setMaterial(cMaterial *ptrMaterial, const int &domain = 0) = 0;

  //! @brief return one of the materials of the coupled domains
  virtual cMaterial *getMaterial(const int &domain = 0) const = 0;

  //! @brief following functions needed for fluid flow coupling
  inline void setStrId(const int Index) { StrId = Index; }
  inline void setFluId(const int Index) { FluId = Index; }

  //! @brief following functions needed for fluid flow coupling
  inline int getStrId(void) { return StrId; }
  inline int getFluId(void) { return FluId; }

  //! @brief add the contribution of the interfaceelement to the global
  //! matrix and vector. This is the version for the frequency domain analysis.
  //! @param omega angular frequency
  //! @param K global stiffness matrix
  //! @param F global load vector (used when inhomogenous boundary conditions
  //! are applied)
  //! @todo extend this routine to non-homogenous boundary conditions
  virtual void addContributionFrequencyDomain(const PetscReal &omega, Mat &K,
                                              Vec &F) = 0;

  //! @brief add the contribution of the interfaceelement to the global
  //! matrix and vector. This is the version for the frequency domain analysis
  //! with velocity potential formulation
  //! @param omega angular frequency
  //! @param K global stiffness matrix
  //! @param F global load vector (used when inhomogenous boundary conditions
  //! are applied)
  //! @author Harikrishnan Sreekumar
  virtual void addContributionFrequencyDomainVelocityPotential(
      const PetscReal &omega, Mat &A, Mat &B, Vec &F) = 0;

  //! @brief add the contribution of the interfaceelement to the global
  //! matrix and vector. This is the version for the time domain/eigenvalue
  //! analysis.
  //! @param K global stiffness matrix
  //! @param M global mass matrix
  //! @todo extend this routine to non-homogenous boundary conditions
  virtual void addContributionTimeDomain(Mat &K, Mat &M) = 0;

  //! @brief write this object to a stream
  //! @param os outputstream
  //! @return der modified outputstream
  std::ostream &write(std::ostream &os) const;

  //! @brief write this object in XML format to a stream
  //! @param os outputstream
  //! @return der modified outputstream
  virtual std::ostream &writeXml(std::ostream &os) const = 0;
};

//! @brief overloaded outputoperator
inline std::ostream &write(std::ostream &os, const cElementInterface &other) {
  return other.write(os);
}

#endif
