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

#ifndef INFAM_ANALYSIS_H
#define INFAM_ANALYSIS_H

#include <set>

#ifdef FRINK_HAVE_BOOST
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/minimum_degree_ordering.hpp>
#include <boost/graph/properties.hpp>

#endif

#include "../element/structure/linear/spring/elementstructurespring.h"
#include "../element/structure/linear/spring/elementstructurespringbc.h"
#include "../fedata/mesh.h"
#include "../misc/nodalforce/nodalforcefluid.h"
#include "../misc/nodalforce/nodalforcestructure.h"
#include "../misc/nodalmoment/nodalmomentstructure.h"
#include "../misc/postprocess/output.h"
#include "slepceps.h"

#ifdef PETSC_HAVE_PARMETIS
#include "parmetis.h"
#endif

//! @brief flag to differ between different analysis types
enum eAnalysisType {
  Static = 1,
  Time = 2,
  Frequency = 3,
  Eigenvalue = 4,
  GeoOpt = 5,
  ModRed = 6,
  ParaModRed = 7
};

//! @brief holder for parsed analysis entities
struct sAnalysisEntities {
  // basic, frequency, time
  double start = 0., delta = 0., alpha = 0., beta = 0., factor = 0., gamma = 0.,
         raya = 0., rayb = 0., alpham = 0., end = 0.;
  int steps = 1;
  char swebem[50] = "none";
  char sbfem[50] = "none";
  char similarv4[50] = "none";
  char tbem[50] = "none";
  // eigen
  int computeEigenvectors;
  // geoopt
  int geoOptStep = 0;
  double geoOptRR0 = 0., geoOptER = 0., geoOptMaxRR = 0.;
  char geoOptType[50] = "none";
  char geoOptStressType[50] = "none";
  // mor-offline
  std::string morMethod, morSetting, morSystem, morEstimator, errorQuantity;
  double sigmaTol, errorTol;
  int maxOrder, pointLoading, constantLoading, freqDepLoading;
  int morSteps = 1;
  double morStart, morEnd;
  std::vector<double> inputnodes, outputnodes;
  std::vector<PetscInt> activedofs;
  // pmor-offline
  int numFoms, numRoms;
  std::vector<std::string> filenameFomsHdf5, filenameRomsHdf5;
};

//! @brief Virtual class, foundation of all analysis types
//! @author Dirk Clasen
//! @date 09.08.2005
class cAnalysis : public cOutput {
 private:
  ///! we want to have a symmetric FSI formulation (does not work for all cases)
  bool m_MakeSymmetric;
  ///! compute stresses in element
  bool m_ComputeStressesCells;
  ///! compute stresses at nodes
  bool m_ComputeStressesNodes;

  //! @brief applies a given permutation to the numbering of matrix' entries
  //! @param MyMesh mesh to whom the permutation has to be applied
  //! @param ordering permutation computed by boost::graph or PETSc
  void permuteMatrixEntries(cMesh &MyMesh, AO &ordering);

  //! @brief determine the local sparsity pattern of the system matrix.
  //! @param MyMesh mesh to whom the permutation has to be applied
  void getLocalSparsityPattern(cMesh &myMesh, const PetscInt &FirstRow,
                               const PetscInt &LastRow,
                               std::vector<std::set<PetscInt> > &Pattern);

  //! @brief apply matrix reordering scheme to the system matrix.
  //! use ParMETIS to compute fill-in reducing ordering
  void renumberMatrixUsingParmetis(cMesh &myMesh, const PetscInt &FirstRow,
                                   const PetscInt &LastRow,
                                   std::vector<std::set<PetscInt> > &Pattern);

  //! @brief apply matrix reordering scheme to the system matrix.
  //! use boost's rcm algorithm to reduce bandwidth
  void renumberMatrixUsingBoost(cMesh &myMesh, const PetscInt &FirstRow,
                                const PetscInt &LastRow,
                                std::vector<std::set<PetscInt> > &Pattern);

  //! @brief Distributes a mesh among the processes of PETSC_COMM_WORLD.
  //! First, each process reads the whole mesh. Afterwards, a part of the mesh
  //! is determined that has to be evaluated by a single process. The other
  //! elements for which the current process is not responsible are removed
  //! of the cMesh object. This is done for normal elements as well as
  //! interface elements.
  //! @param MyMesh Mesh to be distributed
  void distributeMesh(cMesh &MyMesh);
  void distributeMeshRowdependent(cMesh &MyMesh);

  //! @brief use PETSc preallocation in order to speed-up matrix assembly
  void preallocateMatrix(const PetscInt &FirstRow, const PetscInt &LastRow,
                         std::vector<std::set<PetscInt> > &Pattern);

  //! @brief parMetisInterface (Prashanth Sheshappa - Studienarbeit)
  void reorderUsingParMetis(cMesh &myMesh, PetscInt optionFlag,
                            std::vector<PetscInt> &Ele2del);
  void distributeParMetis(cMesh &myMesh, std::vector<PetscInt> &Ele2del);

  //! @brief Insert nodal forces into the loadvector
  void applyNodalForces(cMesh &MyMesh, const PetscReal &factor,
                        const PetscInt &step = 0);

  //! @brief Insert nodal moments into the loadvector
  void applyNodalMoments(cMesh &MyMesh, const PetscReal &factor,
                         const PetscInt &step = 0);

 protected:
  /// specifies analysis type
  eAnalysisType m_AnalysisType;

  /// number of rows/columns of system matrices and vectors
  PetscInt m_NumberOfUnknowns;

  /// becomes true if problem contains linear elements or materials
  bool m_nl;
  bool m_elementsDistributed;

  Mat m_K;        ///< global stiffnessmatrix
  Mat m_C;        ///< global dampingmatrix
  Mat m_M;        ///< global massmatrix
  Vec m_x;        ///< solution of linear system of equations A*x=F
  Vec m_delta_x;  ///< solution of nonlinear system
  Vec m_F;        ///< right-hand-side / loadvector
  EPS m_Eps;      ///< solves for eigenvalues
  KSP m_KSP;      ///< PETSc's linear equation solver
  PC m_PC;        ///< preconditioner-

  Vec m_StressesCells;      ///< vector storing stresses at elements
  Vec m_StressesCellsSec2;  ///< vector storing stresses at elements (second
                            ///< section for shell elements)
  Vec m_StressesNodes;      ///< vector storing stresses at nodes

  //! @brief Looking for linear or nonlinear behavior of our problem.
  //! Sets m_linear true if problem contains linear elements or materials.
  //! Sets m_nonlinear true if problem contains nonlinear elements or materials.
  //! @param MyMesh Mesh to be distributed
  void checkForLinearAndOrNonlinearElementsAndMaterials(cMesh &MyMesh);

  //! @brief function to assemble global tensors (stiffness matrix, mass matrix,
  //! load vector)
  //! @param myMesh used mesh
  //! @param nmodes used by similar sv4
  //! @param omega angular freq.
  virtual void assembleGlobalTensors(cMesh &MyMesh, int nmodes = 0,
                                     const PetscReal omega = 0) = 0;

  //! @brief Marks degrees of freedom as active for the nodes of the
  //! discretization. Afterwards, these degrees of freedom are numbered
  //! sequentially.
  //! @param myMesh used mesh
  void activateDofsAtNodes(cMesh &myMesh);

  //! @brief Function to count all elements which are attached to one node
  //! @param myMesh used mesh
  void countElementsPerNode(cMesh &myMesh);

  //! @brief optimizes the pattern of the systemmatrix according to the nodal
  //! degrees of freedom. You can use several reordering schemes to reduce
  //! the bandwidth of the matrix as well as fill-in of matrix factorization.
  void optimizeMatrixPattern(cMesh &myMesh);

  //! @brief allocates memory and initializes the solvers.
  virtual void initializePETScObjects(cMesh &myMesh) = 0;

  //! @brief frees memory and destroys the solver objects.
  virtual void deletePETScObjects(void) = 0;

  //! @brief initialize the PETSc KSP solver
  void initializeSolverObjects(void);

  //! @brief inserts a single element matrix into a global matrix
  //! @param ptrElement pointer to the element that shall be inserted
  //! @param EM element matrix
  //! @param SysMatrix systemmatrix to which EM has to be added
  void insertElementMatrix(cElementFEM *ptrElement, cElementMatrix &EM,
                           Mat &SysMatrix);

  //! @brief Adds the contribution of an elementvector V of element ptrElement
  //! to a global vector SysVector
  void insertElementVector(cElementFEM *ptrElement, cElementVector &V,
                           Vec &SysVector);

  //! @brief Writes a global matrix to a text file in Matlab's format.
  //! @param SysMatrix matrix to write
  //! @param Filename name of outputfile
  //! @param Identifier name of the matrix within Matlab
  void dumpMatrixToFile(Mat &SysMatrix, const std::string &Filename,
                        const std::string &Identifier);

  //! @brief Writes a global vector to a text file in Matlab's format.
  //! @param SysVector vector to write
  //! @param Filename name of outputfile
  //! @param Identifier name of the matrix within Matlab
  void dumpVectorToFile(Vec &SysVector, const std::string &Filename,
                        const std::string &Identifier);

  //! @brief Solves a linear equation system e.g. m_K * m_x = m_F
  void solveEquationSystem(Mat &matrix, Vec &solution, Vec &vec,
                           bool resetPreconitioner = true);

  //! @brief Insert nodal forces and moments into the loadvector
  void applyNodalLoads(cMesh &MyMesh, const PetscReal &factor,
                       const PetscInt &step = 0);

  //! @brief function to compute stresses in the element
  //! @param MyMesh mesh used for computation
  void computeStressesCells(cMesh &MyMesh);

  //! @brief function to compute stresses at the nodes
  //! @param MyMesh mesh used for computation
  void computeStressesNodes(cMesh &MyMesh);

  //! @brief function to compute stresses
  //! sigma_v = sqrt( sigma_11^2 + sigma_22^2 + sigma_33^2
  //!               - sigma_11*sigma_22 - sigma_11*sigma_33 - sigma_22*sigma_33
  //!               + 3*(sigma_12^2 + sigma_13^2 + sigma_23^2)
  cElementVector getStressesVanMieses(cMesh &MyMesh);

  //! @brief function to compute plane stresses
  //!                |s_11-s s_12   s_13  |
  //! s_{1,2,3} = det|s_21   s_22-s s_23  |=0
  //!                |s_31   s_32   s-33-s|
  cElementVector getStressesPlaneStress(cMesh &MyMesh);

  //! @brief Insert boundary conditions into global matrices and load vector
  //! @param setBCValue2zero if it is true, BC value will be set to zero. The
  //! default is flase.
  //! @param factor is used to rate the BC value. The default is 1.0
  //! @param setBCMassMat if it true, the BC is introduced to mass matrix as
  //! well - used for eigenvalue analysis
  bool insertBoundaryConditions(PetscBool setBCValue2zero = PETSC_FALSE,
                                PetscReal factor = 1.0,
                                PetscBool setBCMassMat = PETSC_FALSE);

  //@brief this is used to apply boundary conditions
  void getGlobalBCs(cMesh &MyMesh, const PetscReal omega = 0);
  PetscInt m_NumFixedRows;  ///< number of fixed rows
  PetscInt *m_RowsBCDofs;   ///< list of rows fixed
  PetscScalar *m_BCValues;  ///< respective bc values

  // --- this is used to apply multi point constraints
  Mat m_Gamma;  ///< global transformation matrix

 public:
  //! @brief Constructor
  cAnalysis();
  //! @brief Destructor
  virtual ~cAnalysis();

  //! @brief Returns the total number of unknowns
  //! @return number of total unknowns
  PetscInt getNumberOfUnknowns() { return m_NumberOfUnknowns; }

  //! @brief set flag that tells us if we want a symmetric formulation or not
  void setFlagSymmetricFormulation(bool flag) { m_MakeSymmetric = flag; }

  //! @brief check if symmetric fsi coupling formulation is requested
  bool checkForSymmetricFormulation(void) const { return m_MakeSymmetric; }

  //! @brief set flag that tells us if stresses in the elements shall be
  //! computed
  void setFlagStressCellsComputation(bool flag) {
    m_ComputeStressesCells = flag;
  }

  //! @brief set flag that tells us if stresses at the nodes shall be computed
  void setFlagStressNodesComputation(bool flag) {
    m_ComputeStressesNodes = flag;
  }

  //! @brief check if stress computation in the elements is requested
  bool checkForStressCellsComputation(void) const {
    return m_ComputeStressesCells;
  }

  //! @brief check if stress computation at the nodes is requested
  bool checkForStressNodesComputation(void) const {
    return m_ComputeStressesNodes;
  }

  //! @brief Perform a full computation on the given mesh.
  //! @param MyMesh mesh used for computation
  //! @param freemem tells us whether the Petsc Objects are going to be deletet
  //! at the and of FullRun. The default is TRUE
  virtual void FullRun(cMesh &MyMesh, PetscBool freemem = PETSC_TRUE) = 0;

  //! @brief Write this object to a stream.
  //! @param os outputstream
  //! @return der modified outputstream
  virtual std::ostream &write(std::ostream &os) const = 0;

  //! @brief Writes this element in XML format to a stream.
  //! Used for writing the mesh to a new inputfile after searching
  //! for the interface elements.
  //! @param os outputstream
  //! @return der modified outputstream
  virtual std::ostream &writeXml(std::ostream &os) const = 0;

  //! @brief Read data for this object from a stream.
  //! @param is inputstream
  //! @return der modified inputstream
  virtual std::istream &read(std::istream &is) = 0;

  //! @brief displacement vector
  Vec *getSolution() { return &m_x; }

  //! @brief stress vectors
  Vec *getStressesCells() { return &m_StressesCells; }
  Vec *getStressesCellsSec2() { return &m_StressesCellsSec2; }
  Vec *getStressesNodes() { return &m_StressesNodes; }

  //! @brief shut down the PETSc KSP context
  void deleteSolverObjects(void);
};

//! @brief overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os, const cAnalysis &other) {
  return other.write(os);
}

//! @brief overloaded inputoperator
inline std::istream &operator>>(std::istream &is, cAnalysis &other) {
  return other.read(is);
}

#endif
