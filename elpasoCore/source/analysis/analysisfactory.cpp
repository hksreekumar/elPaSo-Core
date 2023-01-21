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

#include "analysisfactory.h"

#include "../misc/parser/femparserinterface.h"
#include "../misc/postprocess/outputavs.h"
#include "../misc/postprocess/outputba.h"
#include "../misc/postprocess/outputfile.h"
#include "../misc/postprocess/outputlev.h"
#include "../misc/postprocess/outputstp.h"
#include "../misc/postprocess/outputstp2.h"
#include "../misc/postprocess/outputtec.h"
#include "../misc/postprocess/outputvtk.h"
#include "./analysis.h"
#include "./frequency/analysisfrequencybasic.h"

cAnalysisFactory::cAnalysisFactory() {
  // empty
}

cAnalysisFactory::~cAnalysisFactory() {
  // empty
}

cAnalysis* cAnalysisFactory::createAnalysis(std::string _analysisIdentifier,
                                            cFemParserInterface* _parser,
                                            cMesh& _mesh) {
  cAnalysis* ptr = 0;
  std::stringstream data_stream;

  if (_analysisIdentifier == "frequency-basic") {
    ptr = new cAnalysisFrequencyBasic;
    sAnalysisEntities analysis_data = _parser->getAnalysisData();
    data_stream << analysis_data.start << " " << analysis_data.steps << " "
                << analysis_data.delta << " " << analysis_data.swebem << " "
                << analysis_data.sbfem << " " << analysis_data.similarv4 << " ";
    data_stream >> *ptr;
  }

  return ptr;
}

void cAnalysisFactory::prepareAnalysisOutputs(cAnalysis* _analysis,
                                              cFemParserInterface* _parser) {
  ParserOutputData readOutputData = _parser->getOutputData();
  // set flag for output of step files
  cOutputSTP* ptrOutputSTP = new cOutputSTP();
  ptrOutputSTP->setWantOutput(readOutputData.writestp);
  _analysis->insertOutput(ptrOutputSTP);

  // set flag for output of step files
  cOutputSTP2* ptrOutputSTP2 = new cOutputSTP2();
  ptrOutputSTP2->setWantOutput(readOutputData.writestp2);
  _analysis->insertOutput(ptrOutputSTP2);

  // set flag for output of vtk files
  cOutputVTK* ptrOutputVTK = new cOutputVTK();
  ptrOutputVTK->setWantOutput(readOutputData.writevtk);
  _analysis->insertOutput(ptrOutputVTK);
}