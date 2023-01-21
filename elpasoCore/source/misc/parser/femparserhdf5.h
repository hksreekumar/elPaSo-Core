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

#ifndef ELPASO_FEMPARSERHDF5_H
#define ELPASO_FEMPARSERHDF5_H

#include "./femparserinterface.h"

namespace HDF5 {
class Handler;
class ElementHandler;
class MaterialHandler;
}  // namespace HDF5

//! @brief HDF5 parser that reads the given discretisation
//! @author Harikrishnan Sreekumar
//! @date 21.10.2020
class cFemParserHDF5 : public virtual cLogging, public cFemParserInterface {
 public:
  //! @brief Constructor
  //! @author Harikrishnan Sreekumar
  //! @date 21.10.2020
  cFemParserHDF5();

  //! @brief Destructor
  //! @author Harikrishnan Sreekumar
  //! @date 21.10.2020
  ~cFemParserHDF5();

  //! @brief Function to open the desired input file
  //! @author Harikrishnan Sreekumar
  //! @date 06.12.2022
  void openInputFile(std::string _filename);

  //! @brief closes the file handle
  //! @date 07.11.2022
  //! @author Harikrishnan K. Sreekumar
  void closeInputFile();

  //! @brief parse section Analysis for HDF5
  //! @date 06.12.2022
  //! @author Harikrishnan K. Sreekumar
  void parseAnalysisEntity(cProblem& userData);

  //! @brief Returns the analysis type
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  virtual std::string getAnalysisType();

  //! @brief Returns the analysis data
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  sAnalysisEntities getAnalysisData();

  //! @brief Returns the output data
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  ParserOutputData getOutputData();

  //! @brief reads the orientation file
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  void readOrientationFromFile(void);

  //! @brief parse section Description
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  void parseDescriptionEntity(cProblem& userData);

  //! @brief parse section Revision
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  void parseRevisionEntity();

  //! @brief Returns the noda data
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  ParserNodeData getNodalData();

  //! @brief Returns the number of materials
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  int getNumberOfMaterials();

  //! @brief Returns the specific material data
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  ParserMaterialData getMaterialData(int _materialId);

  //! @brief parse section material for HDF5
  //! @date 26.04.2021
  //! @author Harikrishnan K. Sreekumar
  static void sMaterialAttributeHandler(HDF5::MaterialHandler* _handlers,
                                        int _size);

  //! @brief Returns the number of element blocks
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  int getNumberOfElementBlocks();

  //! @brief Returns the specific element block data
  //! @param _blockId Block ID in sequential order
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  ParserElementBlockData getElementBlockData(int _blockId);

  //! @brief Returns the number of node constraints
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  int getNumberOfNodeConstraints();

  //! @brief Returns the type of node constraints
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  std::string getNodeConstraintsType(int _id);

  //! @brief Returns node constraints - structure data
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  ParserNodeConstraintStructureData getNodeConstraintsStructureData(int _id);

  //! @brief Returns node constraints - acoustic data
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  ParserNodeConstraintAcousticData getNodeConstraintsAcousticData(int _id);

  //! @brief Returns node constraints - ids data
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  ParserNodeConstraintIdsData getNodeConstraintsIdsData(int _id);

  //! @brief Returns the number of element loads
  //! @date 06.01.2022
  //! @author Harikrishnan K. Sreekumar
  int getNumberOfElementLoads();

  //! @brief Returns element loads data
  //! @date 06.01.2022
  //! @author Harikrishnan K. Sreekumar
  ParserElementLoadsData getElementLoadsData(int _id);

  //! @brief Returns the number of node loads
  //! @date 06.11.2022
  //! @author Harikrishnan K. Sreekumar
  int getNumberOfNodeLoads();

  //! @brief Returns node loads data
  //! @date 07.11.2022
  //! @author Harikrishnan K. Sreekumar
  ParserNodeLoadsData getNodeLoadsData(int _id);

  //! @brief Returns the number of node moments
  //! @date 07.11.2022
  //! @author Harikrishnan K. Sreekumar
  int getNumberOfNodeMoments();

  //! @brief Returns node moments data
  //! @date 07.11.2022
  //! @author Harikrishnan K. Sreekumar
  ParserNodeMomentsData getNodeMomentsData(int _id);

  //! @brief Returns the number of interface elements
  //! @date 07.11.2022
  //! @author Harikrishnan K. Sreekumar
  int getNumberOfInterfaceElements();

  //! @brief Returns the number of interfaces
  //! @date 07.11.2022
  //! @author Harikrishnan K. Sreekumar
  int getNumberOfInterfaces();

  //! @brief Returns interface element data
  //! @date 07.11.2022
  //! @author Harikrishnan K. Sreekumar
  ParserInterfaceElementsData getInterfaceElementsData(int _id);

  //! @brief Returns the number of NC interface elements
  //! @date 05.01.2023
  //! @author Harikrishnan K. Sreekumar
  int getNumberOfNCInterfaceElements();

  //! @brief Returns the number of NC interfaces
  //! @date 05.01.2023
  //! @author Harikrishnan K. Sreekumar
  int getNumberOfNCInterfaces();

  //! @brief Returns NC interface element data
  //! @date 05.01.2023
  //! @author Harikrishnan K. Sreekumar
  ParserNCInterfaceElementsData getNCInterfaceElementsData(int _id);
};
#endif