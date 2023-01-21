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

#ifndef ELPASO_FEMPARSERINTERFACE_H
#define ELPASO_FEMPARSERINTERFACE_H

#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include "../../analysis/analysisfactory.h"
#include "../../basics/exceptions/exceptions.h"
#include "../../bc/boundaryconditionfactory.h"
#include "../../element/elementfactory.h"
#include "../../element/interface/elementinterfacefactory.h"
#include "../../element/load/elementloadfactory.h"
#include "../../element/ncinterface/ncelementinterfacefactory.h"
#include "../../fedata/entities.h"
#include "../../fedata/problem.h"
#include "../../material/materialfactory.h"
#include "../../misc/nodalforce/nodalforcefactory.h"
#include "../../misc/nodalmoment/nodalmomentfactory.h"
#include "../log/logging.h"

class sAnalysisEntities;

//! @brief holder for output data
struct ParserOutputData {
  int writestp = 0;
  int writestp2 = 0;
  int writelev = 0;
  int writetec = 0;
  int writeavs = 0;
  int writevtk = 0;
  int writeak3 = 0;
  int writeenergy = 0;
};

//! @brief holder for nodal data
struct ParserNodeData {
  std::vector<double> nodalCoordinates;
  std::vector<int> nodalIds;
};

//! @brief maximum material tags
enum {
  num_mat_tags = 97 /*!< number of possible material parameter tags (listed in
                       eMaterialTags) */
};

//! @brief holder for material data
struct ParserMaterialData {
  std::string matType;
  std::string matName = "unknown";
  std::string matParameters[num_mat_tags];
};

//! @brief holder for element blocks
struct ParserElementBlockData {
  PetscInt blockNumElements;
  PetscInt blockId;
  PetscInt blockMaterial;
  std::string blockName;
  std::string blockOrientation;
  std::string blockOrientationFilename;
  std::string blockElementType;
  std::vector<double> blockConnectiviy;  // matrix with ids and nodes
  PetscInt blockConnectiviyNumRows;
  PetscInt blockConnectiviyNumCols;
};

//! @brief holder for node constraints - ids
struct ParserNodeConstraintIdsData {
  PetscInt nconstraintNodeId;
  PetscInt nconstraintBcId;
};

//! @brief holder for node constraints - structre
struct ParserNodeConstraintStructureData {
  std::string nconstraintName;
  PetscInt nconstraintId;
  const PetscInt tagcount = 20;
  std::vector<int> nconstraintFlag;
  std::vector<double> nconstraintVals;
};

//! @brief holder for node constraints - structre
struct ParserNodeConstraintAcousticData {
  std::string nconstraintName;
  PetscInt nconstraintId;
  double nconstraintValue;
};

//! @brief holder for element loads
struct ParserElementLoadsData {
  std::string eloadType;
  std::string eloadDataset;
  PetscInt eloadId;
  PetscInt eloadElementId;
  double eloadVelocity;
  double eloadFace;
};

//! @brief holder for node loads
struct ParserNodeLoadsData {
  std::string nloadType;
  PetscInt nloadId;
  PetscInt nloadNodeId;
  PetscInt nloadSteps;
  double nloadDeltaT;
  double nloadFx;
  double nloadFy;
  double nloadFz;
  double nloadPf;
};

//! @brief holder for node moments
struct ParserNodeMomentsData {
  std::string nmomentType;
  PetscInt nmomentId;
  PetscInt nmomentNodeId;
  PetscInt nmomentSteps;
  double nmomentDeltaT;
  double nmomentMx;
  double nmomentMy;
  double nmomentMz;
};

//! @brief holder for interface elements
struct ParserInterfaceElementsData {
  std::string Type;
  std::vector<int> groupdata;
  int groupdata_m;
  int groupdata_n;
  int MatF, MatS, StrId, FluId;
  int nnod_NF, nnod_NS;
};

//! @brief holder for NC interface elements
struct ParserNCInterfaceElementsData {
  std::string Type;
  std::vector<int> groupdata;
  std::vector<double> matvec_coord;
  std::vector<int> vec_ids;
  int groupdata_m;
  int groupdata_n;
  int MatF, MatS, StrId, FluId;
  int nnod_NF, nnod_NS;
};

/* BEGIN_NO_COVERAGE*/
//! @brief struct for holding interface element data
struct AppInterfaceElements {
  int* dataspace;
  int MatF;
  int MatS;
  int dataspace_size;
  int nnod;
  std::vector<double> _matvec_coord_elem;
  std::string dataset;
};
/* END_NO_COVERAGE*/

//! @brief number of expected entries for a specific part of the inputfile
//! @date 09.11.2005
//! @author Dirk Clasen
struct EstimatedNumberOfEntities {
  bool AttributesExistInInputFile; /*!< says us, if the attribute "number"
                                      exists */
  int Interfaces;                  /*!< number of interfaces */
  int Nodes;                       /*!< number of nodes */
  int Materials;                   /*!< number of materials */
  int Elements;                    /*!< number of elements */
  int NodalForces;                 /*!< number of nodal forces */
  int LoadedNodes;                 /*!< number of nodes with associated load */
  int NodalMoments;                /*!< number of nodal moments */
  int MomentedNodes;       /*!< number of nodes with associated moment */
  int NodalPressures;      /*!< number of nodal pressures */
  int NodalValues;         /*!< number of nodal values (from fluid flow) */
  int PressuredNodes;      /*!< number of nodes with associated pressure */
  int NodeBCs;             /*!< number of nodal boundaryconditions */
  int FixedNodes;          /*!< number of fixed nodes */
  int ElementLoads;        /*!< number of elementloads */
  int LoadedElements;      /*!< number of loaded elements */
  int InterfaceElements;   /*!< number of interface elements used to describe
                              interaction btw. domains */
  int NCInterfaceElements; /*!< number of non-conforming interface elements used
                              to describe interaction btw. domains */
  int ElementsFF;          /*!< number of fluid flow elements */
  int Constraints;         /*!< number of constraints */
};

//! @brief enumeration of all possible tags processed by XML parser
enum eMaterialTags {
  Id,      /*!< material id */
  M,       /*!< pointmass */
  E,       /*!< Young's modulus */
  Ex,      /*!< Young's modulus */
  Ey,      /*!< Young's modulus */
  Ez,      /*!< Young's modulus */
  ExMem,   /*!< Young's modulus for membrane part */
  EyMem,   /*!< Young's modulus for membrane part */
  nu,      /*!< Poisson's ratio */
  nuxy,    /*!< Poisson's ratio */
  nuxz,    /*!< Poisson's ratio */
  nuyz,    /*!< Poisson's ratio */
  nuxyMem, /*!< Poisson's ratio for membrane part */
  Gxy,     /*!< shear modulus */
  Gxz,     /*!< shear modulus */
  Gyz,     /*!< shear modulus */
  GxyMem,  /*!< shear modulus for membrane part */
  A,       /*!< cross section */
  Ix,      /*!< geometrical moment of inertia */
  Iy,
  Iz,
  rho,   /*!< density */
  rhof,  /*!< fluid density */
  c,     /*!< speed of sound */
  cf,    /*!< speed of sound of fluid part in multiphase medium (equiv fluid) */
  t,     /*!< thickness */
  Fi,    /*!< initial force to prestress element */
  preX,  /*!< initial stress in x to prestress DSG9pre and PlShell9pre elements
          */
  preY,  /*!< initial stress in x to prestress DSG9pre and PlShell9pre elements
          */
  mtyp,  /*!< type (viscoelastic materials) */
  eta,   /*!< damping constant (viscoelastic materials) */
  G,     /*!< shear modulus */
  alpha, /*!< Biot factor */
  alinf, /*!< tortuosity alpha_inf */
  phi,   /*!< porosity */
  R,     /*!< flow resistivity */
  K,     /*!< compressibility modulus */
  Ks,    /*!< compressibility modulus structure (poro)*/
  Kf,    /*!< compressibility modulus fluid (poro)*/
  kappa, /*!< */
  ks,    /*!< porous material tortuosity */
  RL,    /*!< porous medium flow resistivity */
  st,    /*!< thermal pore shape factor */
  sv,    /*!< viscous pore shape factor */
  at,    /*!< */
  sigma, /*!< */
  lambda,  /*!< */
  lambdat, /*!< */
  factorE, /*!< factor for Young's modulus in bi-linear hardening material law
            */

  /**
   * Material parameters for the cap model. (plane cap)
   * The notation of the parameters was adopted from the book "inelastic
   * analysis of solids and structures, Kojic/Bathe, 2005", and additionally the
   * letters "_cm" (for cap model) were added. The description and determination
   * of the parameters can also readed in the book "constitutive laws for
   * engineering materials, Desai/Siriwardane, 1984".
   */
  phi_cm, /*!< angle of friction, slope of the Drucker-Prager II function (cap
             model) */
  c_cm, /*!< kohesion, the value of the intercept of the f1 function (cap model)
         */
  alpha_cm, /*!< material constant*/
  k_cm,     /*!< material constant*/
  A_cm,     /*!< material constant*/
  D_cm,     /*!< the volumetric strain rate */
  W_cm,     /*!< the maximum plastic volumetric strain */
  X_cm,     /*!< initial cap position */

  /**
   * Material parameters for anisotropic cloaking material around a sphere
   * the formulation is according to Cummer et al.: Scattering Theory Derivation
   * of a 3D Acoustic Cloaking Shell, Physical Review Letters, Volume 100,
   * 024301 (2008) [0031-9007/08/100(2)/024301(4)] Letters "_cl" added for
   * cloaking material
   */

  R_i_cl, /*!< inner radius of the cloaking material around the sphere */
  R_o_cl, /*!< outer radius of the cloaking material around the sphere */
  CSx_cl, /*!< center of the sphere the sphere [x]*/
  CSy_cl, /*!< center of the sphere the sphere [y]*/
  CSz_cl, /*!< center of the sphere the sphere [z]*/

  /**
   * Material parameters for spring elements
   */

  Cx,  /*!< stiffness trans x for springs */
  Cy,  /*!< stiffness trans y for springs */
  Cz,  /*!< stiffness trans z for springs */
  Crx, /*!< stiffness rot x for springs */
  Cry, /*!< stiffness rot y for springs */
  Crz, /*!< stiffness rot z for springs */

  /**
   * Material parameters for equivfluiddirect (direct input of complex c and
   * rho)
   */

  creal,   /*!< real part of speed of sound */
  cimag,   /*!< imag part of speed of sound */
  rhoreal, /*!< real part of density */
  rhoimag, /*!< imag part of density */

  /**
   * Material parameters for viscofreqparam and for multiple layer models
   */

  E0,   /*!< basic E0 for formula according to Ereal = E0 + E1*f + E2*f^2 */
  E1,   /*!< E1 for above formula and layer index 1 */
  E2,   /*!< E2 for above formula and layer index 2 */
  E3,   /*!< E3 for layer index 3 */
  Eta0, /*!< Eta0 parameter for damping */
  eta1, /*!< eta1 parameter for damping in layer 1 */
  eta2, /*!< eta2 parameter for damping in layer 2*/
  eta3, /*!< eta3 parameter for damping in layer 3*/
  beta,
  /*!< beta parameter for damping */ /*!< alpha is already defined and used for
                                        the viscofreqparam model also */

  /**
   * Material parameters used in constrained layer damping (CLD) models
   */

  B,    /*!< bending stiffness of a plate element */
  Btyp, /*!< bending stiffness type (beam or plate) */
  nu1,  /*!< Poisson's ratio */
  nu2,  /*!< Poisson's ratio */
  nu3,  /*!< Poisson's ratio */
  G2,   /*!< Shear modulus of viscoelastic middle layer (layer 2) */
  H1,   /*!< thickness of base plate */
  H2,   /*!< thickness of viscoelastic middle layer */
  H3,   /*!< thickness of restraining covering layer */
  rho1, /*!< density rho'x' of layer 'x' */
  rho2,
  rho3,

  wE_RKU,      /*!< weight for RKU homogenized Young's modulus */
  wE_Laminate, /*!< weight for homogenized Young's modulus according to
                  classical laminate theory*/
  wNu_RKU,     /*!< weight for RKU homogenized Poisson ratio */
  wNu_Laminate /*!< weight for homogenized Poisson ratio according to classical
                  laminate theory*/

  /**
   * If you add some parameters to this list - increase the parameter
   * num_mat_tags above !!!
   */
};

/**
 *  @brief Parser Interface to read in required data for elPaSo run
 *  @author Harikrishnan Sreekumar
 *  @date 19.10.2020
 */
class cFemParserInterface : public virtual cLogging {
 public:
  const int m_MAX_CHAR_TAG_LENGTH =
      200;  // define maximum number of characters in an xml tag content string
  int m_Revision;  ///< Revision of input file's format
  PetscInt m_SectionsMaterial;
  PetscInt m_lines_in_orient_file;
  std::string m_SectionsOrientation;
  std::string m_SectionsOrientationFilename;
  cMatrix m_SectionsOrientationVectors;

  //! Object factories
  cAnalysisFactory* m_AnalysisFactory;
  cMaterialFactory* m_MaterialFactory;
  cElementFactory* m_ElementFactory;
  cElementLoadFactory* m_ElementLoadFactory;
  cBoundaryConditionFactory* m_BoundaryConditionFactory;
  cNodalForceFactory* m_NodalForceFactory;
  cNodalMomentFactory* m_NodalMomentFactory;
  cElementInterfaceFactory* m_ElementInterfaceFactory;
  cNCElementInterfaceFactory* m_NCElementInterfaceFactory;

  //! compares the number of read entities with the number of expected entities
  //! given by the section's attribute "number".
  //! @param ptrMesh pointer to the mesh object
  void checkCounters(cMesh* ptrMesh);

  //! reads the orientation file
  virtual void readOrientationFromFile(void) = 0;

  //! stores the number of available entities
  EstimatedNumberOfEntities mEstimatedNumberOfEntities;

 public:
  //! @brief constructor
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  cFemParserInterface();

  //! @brief destructor
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  ~cFemParserInterface();

  //! @brief parse the entire input file
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  void parseInputFile(cProblem& myProblem);

  //! @brief opens the file handle
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual void openInputFile(std::string _filename) = 0;

  //! @brief closes the file handle
  //! @date 07.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual void closeInputFile() = 0;

  //! @brief parse section Description
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual void parseDescriptionEntity(cProblem& userData) = 0;

  //! @brief parse section Revision
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual void parseRevisionEntity() = 0;

  //! @brief parse section Analysis
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  void parseAnalysisEntity(cProblem& userData);

  //! @brief parse section Node
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  void parseNodeEntity(cProblem& userData);

  //! @brief parse section Materials
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  void parseMaterialsEntity(cProblem& userData);

  //! @brief parse section Elements
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  void parseElementsEntity(cProblem& userData);

  //! @brief parse section ElementLoads
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  void parseElementLoadsEntity(cProblem& userData);

  //! @brief parse section NodeConstraints
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  void parseNodeConstraintsEntity(cProblem& userData);

  //! @brief parse section NodeLoads
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  void parseNodeLoadsEntity(cProblem& userData);

  //! @brief parse section NodeMoments
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  void parseNodeMomentsEntity(cProblem& userData);

  //! @brief parse section InterfaceElements
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  void parseInterfaceElementsEntity(cProblem& userData);

  //! @brief parse section NCInterfaceElements
  //! @date 05.01.2023
  //! @author Harikrishnan K. Sreekumar
  void parseNCInterfaceElementsEntity(cProblem& userData);

  //! @brief Returns the analysis type
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual std::string getAnalysisType() = 0;

  //! @brief Returns the analysis data - frequency
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual sAnalysisEntities getAnalysisData() = 0;

  //! @brief Returns the output data
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual ParserOutputData getOutputData() = 0;

  //! @brief Returns the noda data
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual ParserNodeData getNodalData() = 0;

  //! @brief Returns the number of materials
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual int getNumberOfMaterials() = 0;

  //! @brief Returns the specific material data
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual ParserMaterialData getMaterialData(int _materialId) = 0;

  //! @brief Returns the number of element blocks
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual int getNumberOfElementBlocks() = 0;

  //! @brief Returns the specific element block data
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual ParserElementBlockData getElementBlockData(int _blockId) = 0;

  //! @brief Parses the element
  //! @date 22.10.2020
  //! @author Harikrishnan K. Sreekumar
  void processElement(cElementFEM* ptrElement, cMesh* ptrMesh);

  //! @brief Returns the number of node constraints
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual int getNumberOfNodeConstraints() = 0;

  //! @brief Returns the type of node constraints
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual std::string getNodeConstraintsType(int _id) = 0;

  //! @brief Returns node constraints - structure data
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual ParserNodeConstraintStructureData getNodeConstraintsStructureData(
      int _id) = 0;

  //! @brief Returns node constraints - acoustic data
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual ParserNodeConstraintAcousticData getNodeConstraintsAcousticData(
      int _id) = 0;

  //! @brief Returns node constraints - ids data
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual ParserNodeConstraintIdsData getNodeConstraintsIdsData(int _id) = 0;

  //! @brief Returns the number of element loads
  //! @date 06.01.2022
  //! @author Harikrishnan K. Sreekumar
  virtual int getNumberOfElementLoads() = 0;

  //! @brief Returns element loads data
  //! @date 06.01.2022
  //! @author Harikrishnan K. Sreekumar
  virtual ParserElementLoadsData getElementLoadsData(int _id) = 0;

  //! @brief Returns the number of node loads
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual int getNumberOfNodeLoads() = 0;

  //! @brief Returns node loads data
  //! @date 07.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual ParserNodeLoadsData getNodeLoadsData(int _id) = 0;

  //! @brief Returns the number of node moments
  //! @date 07.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual int getNumberOfNodeMoments() = 0;

  //! @brief Returns node moments data
  //! @date 07.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual ParserNodeMomentsData getNodeMomentsData(int _id) = 0;

  //! @brief Returns the number of interface elements
  //! @date 07.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual int getNumberOfInterfaceElements() = 0;

  //! @brief Returns the number of interfaces
  //! @date 07.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual int getNumberOfInterfaces() = 0;

  //! @brief Returns interface element data
  //! @date 07.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual ParserInterfaceElementsData getInterfaceElementsData(int _id) = 0;

  //! @brief Returns the number of NC interface elements
  //! @date 05.01.2023
  //! @author Harikrishnan K. Sreekumar
  virtual int getNumberOfNCInterfaceElements() = 0;

  //! @brief Returns the number of NC interfaces
  //! @date 05.01.2023
  //! @author Harikrishnan K. Sreekumar
  virtual int getNumberOfNCInterfaces() = 0;

  //! @brief Returns NC interface element data
  //! @date 05.01.2023
  //! @author Harikrishnan K. Sreekumar
  virtual ParserNCInterfaceElementsData getNCInterfaceElementsData(int _id) = 0;
};

#endif