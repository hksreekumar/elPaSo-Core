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

#include "./femparserhdf5.h"

#include "../hdf5/rulesh5.h"
#include "../hdf5/inputh5.h"
#include "../hdf5/handlerhdf5.h"
#include "../../analysis/analysisfactory.h"
#include "femparserhdf5.h"

cFemParserHDF5::cFemParserHDF5() : cFemParserInterface() {
    m_AnalysisFactory           = new cAnalysisFactory;
    m_MaterialFactory           = new cMaterialFactory;
    m_ElementFactory            = new cElementFactory;
    m_ElementLoadFactory        = new cElementLoadFactory;
    m_BoundaryConditionFactory  = new cBoundaryConditionFactory;
    m_NodalForceFactory         = new cNodalForceFactory;
    m_NodalMomentFactory        = new cNodalMomentFactory;
    m_ElementInterfaceFactory   = new cElementInterfaceFactory;
    m_NCElementInterfaceFactory = new cNCElementInterfaceFactory;
}

cFemParserHDF5::~cFemParserHDF5() {
    // do nothing
}

void cFemParserHDF5::openInputFile(std::string _filename)
{
    // -----------------------------------------------------------------------
    //   open HDF5 File
    // -----------------------------------------------------------------------
    cInputSingletonH5::getInstance()->setGlobalFileName(_filename);
    cInputSingletonH5::getInstance()->openContainer(ELPASO_H5_READWRITE);
    trace("  HDF5 file handler opened");
}

void cFemParserHDF5::closeInputFile()
{
    // -----------------------------------------------------------------------
    //   close HDF5 File
    // -----------------------------------------------------------------------
    cInputSingletonH5::getInstance()->closeContainer();
    trace("  HDF5 file handler closed");
}

std::string cFemParserHDF5::getAnalysisType() 
{
    return cInputSingletonH5::getInstance()->readStringAttributeFromGroup(nRulesHDF5::m_cAnalysis_sType, nRulesHDF5::m_cAnalysis);
}

sAnalysisEntities cFemParserHDF5::getAnalysisData()
{
    sAnalysisEntities analysisData;

    // basic
    analysisData.start = cInputSingletonH5::getInstance()->readDoubleAttributeFromGroup(nRulesHDF5::m_cAnalysis_sStart, nRulesHDF5::m_cAnalysis);
    analysisData.steps = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sSteps, nRulesHDF5::m_cAnalysis);
    analysisData.delta = cInputSingletonH5::getInstance()->readDoubleAttributeFromGroup(nRulesHDF5::m_cAnalysis_sDelta, nRulesHDF5::m_cAnalysis);

    // time
    analysisData.alpha = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sAlpha, nRulesHDF5::m_cAnalysis);
    analysisData.alpham = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sAlpham, nRulesHDF5::m_cAnalysis);
    analysisData.beta = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sBeta, nRulesHDF5::m_cAnalysis);
    analysisData.gamma = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sGamma, nRulesHDF5::m_cAnalysis);
    analysisData.raya = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sRayA, nRulesHDF5::m_cAnalysis);
    analysisData.rayb = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sRayB, nRulesHDF5::m_cAnalysis);

    std::strcpy(analysisData.sbfem, cInputSingletonH5::getInstance()->readStringAttributeFromGroup(nRulesHDF5::m_cAnalysis_sSbFem, nRulesHDF5::m_cAnalysis).c_str());
    std::strcpy(analysisData.similarv4, cInputSingletonH5::getInstance()->readStringAttributeFromGroup(nRulesHDF5::m_cAnalysis_sSimilarv4, nRulesHDF5::m_cAnalysis).c_str());
    std::strcpy(analysisData.tbem, cInputSingletonH5::getInstance()->readStringAttributeFromGroup(nRulesHDF5::m_cAnalysis_sTBEM, nRulesHDF5::m_cAnalysis).c_str());

    //  eigen
    analysisData.factor = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sFactor, nRulesHDF5::m_cAnalysis);
    analysisData.computeEigenvectors = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sVectors, nRulesHDF5::m_cAnalysis);

    // geoopt
    analysisData.geoOptStep = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sStepGeo, nRulesHDF5::m_cAnalysis);
    analysisData.geoOptRR0 = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sRR0, nRulesHDF5::m_cAnalysis);
    analysisData.geoOptER = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sER, nRulesHDF5::m_cAnalysis);
    analysisData.geoOptMaxRR = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sMAXRR, nRulesHDF5::m_cAnalysis);

    std::strcpy(analysisData.geoOptType, cInputSingletonH5::getInstance()->readStringAttributeFromGroup(nRulesHDF5::m_cAnalysis_sOptType, nRulesHDF5::m_cAnalysis).c_str());
    std::strcpy(analysisData.geoOptStressType, cInputSingletonH5::getInstance()->readStringAttributeFromGroup(nRulesHDF5::m_cAnalysis_sStressType, nRulesHDF5::m_cAnalysis).c_str());

    // mor-offline
    analysisData.morMethod = cInputSingletonH5::getInstance()->readStringAttributeFromGroup(nRulesHDF5::m_cAnalysis_sMorMethod, nRulesHDF5::m_cAnalysis);
    analysisData.morSetting = cInputSingletonH5::getInstance()->readStringAttributeFromGroup(nRulesHDF5::m_cAnalysis_sMorSetting, nRulesHDF5::m_cAnalysis);
    analysisData.morSystem = cInputSingletonH5::getInstance()->readStringAttributeFromGroup(nRulesHDF5::m_cAnalysis_sMorSystem, nRulesHDF5::m_cAnalysis);
    analysisData.morEstimator = cInputSingletonH5::getInstance()->readStringAttributeFromGroup(nRulesHDF5::m_cAnalysis_sMorEstimator, nRulesHDF5::m_cAnalysis);

    analysisData.morStart = cInputSingletonH5::getInstance()->readDoubleAttributeFromGroup(nRulesHDF5::m_cAnalysis_sMorFreqStart, nRulesHDF5::m_cAnalysis);
    analysisData.morSteps = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sMorFreqSteps, nRulesHDF5::m_cAnalysis);
    analysisData.morEnd = cInputSingletonH5::getInstance()->readDoubleAttributeFromGroup(nRulesHDF5::m_cAnalysis_sMorFreqEnd, nRulesHDF5::m_cAnalysis);

    analysisData.errorTol = cInputSingletonH5::getInstance()->readDoubleAttributeFromGroup(nRulesHDF5::m_cAnalysis_sMorErrorTol, nRulesHDF5::m_cAnalysis);
    analysisData.sigmaTol = cInputSingletonH5::getInstance()->readDoubleAttributeFromGroup(nRulesHDF5::m_cAnalysis_sMorSigmaTol, nRulesHDF5::m_cAnalysis);
    analysisData.maxOrder = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sMorMaxOrder, nRulesHDF5::m_cAnalysis);

    analysisData.pointLoading = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sMorLoadingPoint, nRulesHDF5::m_cAnalysis);
    analysisData.constantLoading = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sMorLoadingConst, nRulesHDF5::m_cAnalysis);
    analysisData.freqDepLoading = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sMorLoadingFreq, nRulesHDF5::m_cAnalysis);

    analysisData.errorQuantity = cInputSingletonH5::getInstance()->readStringAttributeFromGroup(nRulesHDF5::m_cAnalysis_sMorErrorQuantity, nRulesHDF5::m_cAnalysis);

    cInputSingletonH5::getInstance()->readDoubleVector(analysisData.inputnodes, nRulesHDF5::m_cNodesets, "", nRulesHDF5::m_cNodesets_sMorInput);
    cInputSingletonH5::getInstance()->readDoubleVector(analysisData.outputnodes, nRulesHDF5::m_cNodesets, "", nRulesHDF5::m_cNodesets_sMorOutput);
    cInputSingletonH5::getInstance()->readIntegerVectorFromAttributeGroup(analysisData.activedofs, nRulesHDF5::m_cAnalysis, nRulesHDF5::m_cAnalysis_vMorActiveDofs);

    // pmor-offline
    analysisData.numFoms = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sPmorNumFom, nRulesHDF5::m_cAnalysis);
    analysisData.numRoms = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sPmorNumRom, nRulesHDF5::m_cAnalysis);

    for (size_t ifile = 0; ifile < analysisData.numFoms; ifile++) {
        analysisData.filenameFomsHdf5.push_back(cInputSingletonH5::getInstance()->readStringAttributeFromGroup(nRulesHDF5::m_cAnalysis_sPmorFom + std::to_string(ifile+1), nRulesHDF5::m_cAnalysis));
    }

    for (size_t ifile = 0; ifile < analysisData.numRoms; ifile++) {
        analysisData.filenameRomsHdf5.push_back(cInputSingletonH5::getInstance()->readStringAttributeFromGroup(nRulesHDF5::m_cAnalysis_sPmorRom + std::to_string(ifile+1), nRulesHDF5::m_cAnalysis));
    }

    return analysisData;
}

ParserOutputData cFemParserHDF5::getOutputData()
{
    ParserOutputData readOutputData;
    readOutputData.writestp2 = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sOutSTP, nRulesHDF5::m_cAnalysis);
    readOutputData.writevtk = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sOutVTK, nRulesHDF5::m_cAnalysis);
    
    return readOutputData;
}

void cFemParserHDF5::readOrientationFromFile(void) 
{
#ifdef HAVE_HDF5
    std::vector<double>  mat_freqvariable;
    int m, n;
    cInputSingletonH5::getInstance()->readDenseDoubleMatrix(mat_freqvariable, m, n, "", "", m_SectionsOrientationFilename);

    m_lines_in_orient_file = m;

    m_SectionsOrientationVectors.resize(m_lines_in_orient_file, 9);

    for (PetscInt k = 0; k < m_lines_in_orient_file; k++)
    {
        for (PetscInt m = 0; m < 9; m++)
        {
            m_SectionsOrientationVectors(k, m) = mat_freqvariable[k * n + m];
        }
    }
#endif
}

void cFemParserHDF5::parseDescriptionEntity(cProblem& userData)
{
    trace("  reading description ...");
    std::string desc;

    desc = cInputSingletonH5::getInstance()->readStringAttributeFromGroup(nRulesHDF5::m_cAnalysis_sDescription, nRulesHDF5::m_cAnalysis);

    if (!desc.length() > 1)
        desc = "no description";

    cMesh* ptr = userData.getMesh();
    ptr->setDescription(desc);
    ptr = 0;
}

void cFemParserHDF5::parseRevisionEntity()
{
    trace("  reading revision ...");
    // --- read the data within Revision tag
    m_Revision = cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cAnalysis_sRevision, nRulesHDF5::m_cAnalysis);

    if (m_Revision < 6)
    {
        message("*******************************************************************\n");
        message("* The format of the input file is not supported by elPaSo anymore.\n");
        message("* The input file's revision is %d, you need at least revision 6.\n", m_Revision);
        message("*******************************************************************\n");
        ExitApp();
    }
}

ParserNodeData cFemParserHDF5::getNodalData()
{
    ParserNodeData readData;
    cInputSingletonH5::getInstance()->readNodeCompoundData(readData.nodalIds, readData.nodalCoordinates, nRulesHDF5::m_cNodes, "", nRulesHDF5::m_cNodes_dFemNodes);

    return readData;
}

int cFemParserHDF5::getNumberOfMaterials()
{
    return cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cMaterials_sN, nRulesHDF5::m_cMaterials);
}

ParserMaterialData cFemParserHDF5::getMaterialData(int _materialId)
{
    ParserMaterialData readMatData;

    std::string datasetname = cInputSingletonH5::getInstance()->getMemberName(nRulesHDF5::m_cMaterials, _materialId - 1);
    readMatData.matType     = cInputSingletonH5::getInstance()->readStringAttributeFromDataset(nRulesHDF5::m_cMaterials_dMaterial_sMaterialType, nRulesHDF5::m_cMaterials, datasetname);

    readMatData.matName     = cInputSingletonH5::getInstance()->readStringAttributeFromDataset(nRulesHDF5::m_cMaterials_dMaterial_sName, nRulesHDF5::m_cMaterials, datasetname);

    readMatData.matParameters[(eMaterialTags)Id] = std::to_string(cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cMaterials_dMaterial_sId, nRulesHDF5::m_cMaterials, datasetname));

    std::vector<HDF5::MaterialHandler> handlers = {
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sPointmass,readMatData.matParameters[(eMaterialTags)M],         datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sE,        readMatData.matParameters[(eMaterialTags)E],         datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sEx,       readMatData.matParameters[(eMaterialTags)Ex],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sEy,       readMatData.matParameters[(eMaterialTags)Ey],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sEz,       readMatData.matParameters[(eMaterialTags)Ez],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sExMem,    readMatData.matParameters[(eMaterialTags)ExMem],     datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sEyMem,    readMatData.matParameters[(eMaterialTags)EyMem],     datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sNu,       readMatData.matParameters[(eMaterialTags)nu],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sNuxy,     readMatData.matParameters[(eMaterialTags)nuxy],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sNuxz,     readMatData.matParameters[(eMaterialTags)nuxz],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sNuyz,     readMatData.matParameters[(eMaterialTags)nuyz],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sNuxyMem,  readMatData.matParameters[(eMaterialTags)nuxyMem],   datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sGxy,      readMatData.matParameters[(eMaterialTags)Gxy],       datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sGxz,      readMatData.matParameters[(eMaterialTags)Gxz],       datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sGyz,      readMatData.matParameters[(eMaterialTags)Gyz],       datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sGxyMem,   readMatData.matParameters[(eMaterialTags)GxyMem],    datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sA,        readMatData.matParameters[(eMaterialTags)A],         datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sI,        readMatData.matParameters[(eMaterialTags)Iy],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sIx,       readMatData.matParameters[(eMaterialTags)Ix],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sIy,       readMatData.matParameters[(eMaterialTags)Iy],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sIz,       readMatData.matParameters[(eMaterialTags)Iz],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sRho,      readMatData.matParameters[(eMaterialTags)rho],       datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sRhof,     readMatData.matParameters[(eMaterialTags)rhof],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sC,        readMatData.matParameters[(eMaterialTags)c],         datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sCf,       readMatData.matParameters[(eMaterialTags)cf],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sT,        readMatData.matParameters[(eMaterialTags)t],         datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sFi,       readMatData.matParameters[(eMaterialTags)Fi],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_spreX,     readMatData.matParameters[(eMaterialTags)preX],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_spreY,     readMatData.matParameters[(eMaterialTags)preY],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sTyp,      readMatData.matParameters[(eMaterialTags)mtyp],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sEta,      readMatData.matParameters[(eMaterialTags)eta],       datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sG,        readMatData.matParameters[(eMaterialTags)G],         datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sAlpha1,   readMatData.matParameters[(eMaterialTags)alpha],     datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sAlpha,    readMatData.matParameters[(eMaterialTags)alpha],     datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sAlinf,    readMatData.matParameters[(eMaterialTags)alinf],     datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sKs,       readMatData.matParameters[(eMaterialTags)Ks],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sKf,       readMatData.matParameters[(eMaterialTags)Kf],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sKappa,    readMatData.matParameters[(eMaterialTags)kappa],     datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sPhi,      readMatData.matParameters[(eMaterialTags)phi],       datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sR,        readMatData.matParameters[(eMaterialTags)R],         datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sKk,       readMatData.matParameters[(eMaterialTags)K],         datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sK,        readMatData.matParameters[(eMaterialTags)K],         datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sks,       readMatData.matParameters[(eMaterialTags)ks],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sRL,       readMatData.matParameters[(eMaterialTags)RL],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sSt,       readMatData.matParameters[(eMaterialTags)st],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sSv,       readMatData.matParameters[(eMaterialTags)sv],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sAt,       readMatData.matParameters[(eMaterialTags)at],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sSigma,    readMatData.matParameters[(eMaterialTags)sigma],     datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sLambda,   readMatData.matParameters[(eMaterialTags)lambda],    datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sLambdat,  readMatData.matParameters[(eMaterialTags)lambdat],   datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sFactorE,  readMatData.matParameters[(eMaterialTags)factorE],   datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sPhicm,    readMatData.matParameters[(eMaterialTags)phi_cm],    datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sCcm,      readMatData.matParameters[(eMaterialTags)c_cm],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sDcm,      readMatData.matParameters[(eMaterialTags)D_cm],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sWcm,      readMatData.matParameters[(eMaterialTags)W_cm],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sXcm,      readMatData.matParameters[(eMaterialTags)X_cm],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sRicl,     readMatData.matParameters[(eMaterialTags)R_i_cl],    datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sRocl,     readMatData.matParameters[(eMaterialTags)R_o_cl],    datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sCSxcl,    readMatData.matParameters[(eMaterialTags)CSx_cl],    datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sCSycl,    readMatData.matParameters[(eMaterialTags)CSy_cl],    datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sCSzcl,    readMatData.matParameters[(eMaterialTags)CSz_cl],    datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sCx,       readMatData.matParameters[(eMaterialTags)Cx],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sCy,       readMatData.matParameters[(eMaterialTags)Cy],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sCz,       readMatData.matParameters[(eMaterialTags)Cz],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sCrx,      readMatData.matParameters[(eMaterialTags)Crx],       datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sCry,      readMatData.matParameters[(eMaterialTags)Cry],       datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sCrz,      readMatData.matParameters[(eMaterialTags)Crz],       datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sCreal,    readMatData.matParameters[(eMaterialTags)creal],     datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sCimag,    readMatData.matParameters[(eMaterialTags)cimag],     datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sRhoreal,  readMatData.matParameters[(eMaterialTags)rhoreal],   datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sRhoimag,  readMatData.matParameters[(eMaterialTags)rhoimag],   datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sE0,       readMatData.matParameters[(eMaterialTags)E0],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sE1,       readMatData.matParameters[(eMaterialTags)E1],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sE2,       readMatData.matParameters[(eMaterialTags)E2],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sE3,       readMatData.matParameters[(eMaterialTags)E3],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sEta0,     readMatData.matParameters[(eMaterialTags)Eta0],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sEta1,     readMatData.matParameters[(eMaterialTags)eta1],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sEta2,     readMatData.matParameters[(eMaterialTags)eta2],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sEta3,     readMatData.matParameters[(eMaterialTags)eta3],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sBeta,     readMatData.matParameters[(eMaterialTags)beta],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sB,        readMatData.matParameters[(eMaterialTags)B],         datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sBtyp,     readMatData.matParameters[(eMaterialTags)Btyp],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sNu1,      readMatData.matParameters[(eMaterialTags)nu1],       datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sNu2,      readMatData.matParameters[(eMaterialTags)nu2],       datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sNu3,      readMatData.matParameters[(eMaterialTags)nu3],       datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sG2,       readMatData.matParameters[(eMaterialTags)G2],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sH1,       readMatData.matParameters[(eMaterialTags)H1],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sH2,       readMatData.matParameters[(eMaterialTags)H2],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sH3,       readMatData.matParameters[(eMaterialTags)H3],        datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sRho1,     readMatData.matParameters[(eMaterialTags)rho1],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sRho2,     readMatData.matParameters[(eMaterialTags)rho2],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sRho3,     readMatData.matParameters[(eMaterialTags)rho3],      datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sWERKU,    readMatData.matParameters[(eMaterialTags)wE_RKU],    datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sWElam,    readMatData.matParameters[(eMaterialTags)wE_Laminate],  datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sWNuRKU,   readMatData.matParameters[(eMaterialTags)wNu_RKU],   datasetname),
        HDF5::MaterialHandler(nRulesHDF5::m_cMaterials_dMaterial_sWNuLam,   readMatData.matParameters[(eMaterialTags)wNu_Laminate], datasetname)
    };// If you add some parameters here - check num_mat_tags

    sMaterialAttributeHandler(&handlers[0], handlers.size());

    return readMatData;
}

void cFemParserHDF5::sMaterialAttributeHandler(HDF5::MaterialHandler* _handlers, int _size)
{
    for (size_t i = 0; i < _size; i++) {
        try {
            _handlers[i].setContainer(cInputSingletonH5::getInstance()->readStringAttributeFromDataset(_handlers[i].m_attributename, nRulesHDF5::m_cMaterials, _handlers[i].m_dataset));
        }
        catch (...) {
            // ignore
        }
    }
}

int cFemParserHDF5::getNumberOfElementBlocks()
{
    return cInputSingletonH5::getInstance()->getNumberOfMembersInGroup(nRulesHDF5::m_cElements);
}

ParserElementBlockData cFemParserHDF5::getElementBlockData(int _blockId)
{
    int iElemGroup = _blockId;
    ParserElementBlockData readData;
    std::string datasetname = cInputSingletonH5::getInstance()->getMemberName(nRulesHDF5::m_cElements, iElemGroup);
    readData.blockNumElements = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cElements_dFemElemGroup_sN, nRulesHDF5::m_cElements, datasetname);
    readData.blockName = cInputSingletonH5::getInstance()->readStringAttributeFromDataset(nRulesHDF5::m_cElements_dFemElemGroup_sName, nRulesHDF5::m_cElements, datasetname);
    readData.blockOrientation = cInputSingletonH5::getInstance()->readStringAttributeFromDataset(nRulesHDF5::m_cElements_dFemElemGroup_sOrientation, nRulesHDF5::m_cElements, datasetname);
    readData.blockOrientationFilename = cInputSingletonH5::getInstance()->readStringAttributeFromDataset(nRulesHDF5::m_cElements_dFemElemGroup_sOrientationFile, nRulesHDF5::m_cElements, datasetname);

    readData.blockId = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cElements_dFemElemGroup_sId, nRulesHDF5::m_cElements, datasetname);
    readData.blockMaterial = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cElements_dFemElemGroup_sMaterial, nRulesHDF5::m_cElements, datasetname);

    readData.blockElementType = cInputSingletonH5::getInstance()->readStringAttributeFromDataset(nRulesHDF5::m_cElements_dFemElemGroup_sElementType, nRulesHDF5::m_cElements, datasetname);

    cInputSingletonH5::getInstance()->readDenseDoubleMatrix(readData.blockConnectiviy, readData.blockConnectiviyNumRows, readData.blockConnectiviyNumCols, nRulesHDF5::m_cElements, "", datasetname);

    return readData;
}

int cFemParserHDF5::getNumberOfNodeConstraints()
{
    return cInputSingletonH5::getInstance()->getNumberOfMembersInGroup(nRulesHDF5::m_cNodeConstraints);
}

std::string cFemParserHDF5::getNodeConstraintsType(int _id)
{
    std::string NodeCon_dataset = cInputSingletonH5::getInstance()->getMemberName(nRulesHDF5::m_cNodeConstraints, _id);
    return cInputSingletonH5::getInstance()->readStringAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sConType, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
}

ParserNodeConstraintStructureData cFemParserHDF5::getNodeConstraintsStructureData(int _id)
{
    ParserNodeConstraintStructureData readData;
    std::string NodeCon_dataset = cInputSingletonH5::getInstance()->getMemberName(nRulesHDF5::m_cNodeConstraints, _id);

    readData.nconstraintName = cInputSingletonH5::getInstance()->readStringAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sConType, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintId = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sId, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);

    readData.nconstraintFlag.resize(readData.tagcount);
    readData.nconstraintFlag[0] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sU1, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[1] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sU2, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[2] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sU3, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[3] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sW1, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[4] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sW2, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[5] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sW3, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[6] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sP0, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[7] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sP1, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[8] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sP2, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[9] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sP3, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[10] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sXD3, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[11] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sWD1, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[12] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sWD2, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[13] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sdwdx, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[14] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sdwdy, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[15] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sdwdxy, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[16] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sFluid, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[17] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sZ1, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[18] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sZ3, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[19] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sX12, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintFlag[20] = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sX22, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);

    readData.nconstraintVals.resize(readData.tagcount);
    readData.nconstraintVals[0] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValU1, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[1] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValU2, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[2] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValU3, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[3] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValW1, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[4] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValW2, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[5] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValW3, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[6] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValP0, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[7] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValP1, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[8] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValP2, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[9] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValP3, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[10] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValXD3, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[11] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValWD1, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[12] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValWD2, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[13] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValdwdx, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[14] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValdwdy, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[15] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValdwdxy, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[16] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValFluid, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[17] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValZ1, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[18] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValZ3, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[19] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValX12, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintVals[20] = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValX22, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);

    return readData;
}

ParserNodeConstraintAcousticData cFemParserHDF5::getNodeConstraintsAcousticData(int _id)
{
    ParserNodeConstraintAcousticData readData;
    std::string NodeCon_dataset = cInputSingletonH5::getInstance()->getMemberName(nRulesHDF5::m_cNodeConstraints, _id);
    
    readData.nconstraintName = cInputSingletonH5::getInstance()->readStringAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sConType, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintId = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sId, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintValue = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sValP, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);

    return readData;
}

ParserNodeConstraintIdsData cFemParserHDF5::getNodeConstraintsIdsData(int _id)
{   
    ParserNodeConstraintIdsData readData;
    std::string NodeCon_dataset = cInputSingletonH5::getInstance()->getMemberName(nRulesHDF5::m_cNodeConstraints, _id);
    readData.nconstraintNodeId = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sNodeId, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    readData.nconstraintBcId = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeConstraints_sNodeConstraint_sId, nRulesHDF5::m_cNodeConstraints, NodeCon_dataset);
    
    return readData;

}

int cFemParserHDF5::getNumberOfElementLoads()
{
    return cInputSingletonH5::getInstance()->getNumberOfMembersInGroup(nRulesHDF5::m_cElemLoads);
}

ParserElementLoadsData cFemParserHDF5::getElementLoadsData(int _id)
{
    ParserElementLoadsData readData;

    std::string ElLoad_datasetname = cInputSingletonH5::getInstance()->getMemberName(nRulesHDF5::m_cElemLoads, _id);
    readData.eloadDataset = ElLoad_datasetname;
    readData.eloadType = cInputSingletonH5::getInstance()->readStringAttributeFromDataset(nRulesHDF5::m_cElemLoads_dFemElemLoad_sLoadType, nRulesHDF5::m_cElemLoads, ElLoad_datasetname);
    readData.eloadId = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cElemLoads_dFemElemLoad_sId, nRulesHDF5::m_cElemLoads, ElLoad_datasetname);
    // read velocity from hdf5
    readData.eloadVelocity = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cElemLoads_dFemElemLoad_sVn, nRulesHDF5::m_cElemLoads, ElLoad_datasetname);
    // read face info from hdf5
    readData.eloadFace = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cElemLoads_dFemElemLoad_sFace, nRulesHDF5::m_cElemLoads, ElLoad_datasetname);

    readData.eloadElementId = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cElemLoads_dFemElemLoad_sElementId, nRulesHDF5::m_cElemLoads, ElLoad_datasetname);
    return readData;
}

int cFemParserHDF5::getNumberOfNodeLoads()
{
    return cInputSingletonH5::getInstance()->getNumberOfMembersInGroup(nRulesHDF5::m_cNodeLoads);
}

ParserNodeLoadsData cFemParserHDF5::getNodeLoadsData(int _id)
{
    ParserNodeLoadsData readData;
    std::string nodeloads_dataset = cInputSingletonH5::getInstance()->getMemberName(nRulesHDF5::m_cNodeLoads, _id);
    readData.nloadType = cInputSingletonH5::getInstance()->readStringAttributeFromDataset(nRulesHDF5::m_cNodeLoads_dFemNodeLoad_sLoadType, nRulesHDF5::m_cNodeLoads, nodeloads_dataset);
    readData.nloadId = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeLoads_dFemNodeLoad_sId, nRulesHDF5::m_cNodeLoads, nodeloads_dataset);
    readData.nloadSteps = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeLoads_dFemNodeLoad_sSteps, nRulesHDF5::m_cNodeLoads, nodeloads_dataset);
    readData.nloadDeltaT = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeLoads_dFemNodeLoad_sDeltaT, nRulesHDF5::m_cNodeLoads, nodeloads_dataset);
    readData.nloadFx = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeLoads_dFemNodeLoad_sFx, nRulesHDF5::m_cNodeLoads, nodeloads_dataset);
    readData.nloadFy = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeLoads_dFemNodeLoad_sFy, nRulesHDF5::m_cNodeLoads, nodeloads_dataset);
    readData.nloadFz = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeLoads_dFemNodeLoad_sFz, nRulesHDF5::m_cNodeLoads, nodeloads_dataset);

    readData.nloadPf = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeLoads_dFemNodeLoad_sP, nRulesHDF5::m_cNodeLoads, nodeloads_dataset);

    readData.nloadNodeId = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeLoads_dFemNodeLoad_sNodeId, nRulesHDF5::m_cNodeLoads, nodeloads_dataset);

    return readData;
}

int cFemParserHDF5::getNumberOfNodeMoments()
{
    return cInputSingletonH5::getInstance()->getNumberOfMembersInGroup(nRulesHDF5::m_cNodeMoments);
}

ParserNodeMomentsData cFemParserHDF5::getNodeMomentsData(int _id)
{
    ParserNodeMomentsData readData;

    std::string nodemoments_dataset = cInputSingletonH5::getInstance()->getMemberName(nRulesHDF5::m_cNodeMoments, _id);
    readData.nmomentType = cInputSingletonH5::getInstance()->readStringAttributeFromDataset(nRulesHDF5::m_cNodeMoments_dFemNodeMoment_sLoadType, nRulesHDF5::m_cNodeMoments, nodemoments_dataset);
    readData.nmomentId = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeMoments_dFemNodeMoment_sId, nRulesHDF5::m_cNodeMoments, nodemoments_dataset);
    readData.nmomentSteps = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeMoments_dFemNodeMoment_snSteps, nRulesHDF5::m_cNodeMoments, nodemoments_dataset);
    readData.nmomentDeltaT = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeMoments_dFemNodeMoment_sDeltaT, nRulesHDF5::m_cNodeMoments, nodemoments_dataset);

    readData.nmomentMx = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeMoments_dFemNodeMoment_sMx, nRulesHDF5::m_cNodeMoments, nodemoments_dataset);
    readData.nmomentMy = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeMoments_dFemNodeMoment_sMy, nRulesHDF5::m_cNodeMoments, nodemoments_dataset);
    readData.nmomentMz = cInputSingletonH5::getInstance()->readDoubleAttributeFromDataset(nRulesHDF5::m_cNodeMoments_dFemNodeMoment_sMz, nRulesHDF5::m_cNodeMoments, nodemoments_dataset);
    readData.nmomentNodeId = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNodeMoments_dFemNodeMoment_sNodeId, nRulesHDF5::m_cNodeMoments, nodemoments_dataset);

    return readData;
}

int cFemParserHDF5::getNumberOfInterfaceElements()
{
    return cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cInterElems_sN, nRulesHDF5::m_cInterElems);
}

int cFemParserHDF5::getNumberOfInterfaces()
{
    return cInputSingletonH5::getInstance()->getNumberOfMembersInGroup(nRulesHDF5::m_cInterElems);
}

ParserInterfaceElementsData cFemParserHDF5::getInterfaceElementsData(int _id)
{
    ParserInterfaceElementsData readData;

    std::string ielem_dataset = cInputSingletonH5::getInstance()->getMemberName(nRulesHDF5::m_cInterElems, _id);
    readData.Type = cInputSingletonH5::getInstance()->readStringAttributeFromDataset(nRulesHDF5::m_cInterElems_dInterfaceElements_sType, nRulesHDF5::m_cInterElems, ielem_dataset);
    cInputSingletonH5::cInputSingletonH5::getInstance()->readDenseIntegerMatrix(readData.groupdata, readData.groupdata_m, readData.groupdata_n, nRulesHDF5::m_cInterElems, "", ielem_dataset);

    readData.MatF = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cInterElems_dInterfaceElements_sMatF, nRulesHDF5::m_cInterElems, ielem_dataset);
    readData.MatS = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cInterElems_dInterfaceElements_sMatS, nRulesHDF5::m_cInterElems, ielem_dataset);

    readData.nnod_NF = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cInterElems_dInterfaceElements_sNfluidNodes, nRulesHDF5::m_cInterElems, ielem_dataset);
    readData.nnod_NS = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cInterElems_dInterfaceElements_sNstructNodes, nRulesHDF5::m_cInterElems, ielem_dataset);
        
    return readData;
}

int cFemParserHDF5::getNumberOfNCInterfaceElements()
{
    return cInputSingletonH5::getInstance()->readIntegerAttributeFromGroup(nRulesHDF5::m_cNCInterElems_sN, nRulesHDF5::m_cNCInterElems);
}

int cFemParserHDF5::getNumberOfNCInterfaces()
{
    return cInputSingletonH5::getInstance()->getNumberOfMembersInGroup(nRulesHDF5::m_cNCInterElems);
}

ParserNCInterfaceElementsData cFemParserHDF5::getNCInterfaceElementsData(int _id)
{
    ParserNCInterfaceElementsData readData;

    std::string ielem_dataset = cInputSingletonH5::getInstance()->getMemberName(nRulesHDF5::m_cNCInterElems, _id);
    readData.Type = cInputSingletonH5::getInstance()->readStringAttributeFromDataset(nRulesHDF5::m_cNCInterElems_dInterfaceElements_sType, nRulesHDF5::m_cNCInterElems, ielem_dataset);

    std::string m_cNCNodes_dNCNodes_group = nRulesHDF5::m_cNCNodes_dNCNodes;// + std::to_string(i);
    cInputSingletonH5::getInstance()->readNodeCompoundData(readData.vec_ids, readData.matvec_coord, nRulesHDF5::m_cNCNodes,"", nRulesHDF5::m_cNCNodes_dNCNodes);
        
    cInputSingletonH5::cInputSingletonH5::getInstance()->readDenseIntegerMatrix(readData.groupdata, readData.groupdata_m, readData.groupdata_n, nRulesHDF5::m_cNCInterElems, "", ielem_dataset);

    // material info
    readData.MatF = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNCInterElems_dInterfaceElements_sMatF, nRulesHDF5::m_cNCInterElems, ielem_dataset);
    readData.MatS = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNCInterElems_dInterfaceElements_sMatS, nRulesHDF5::m_cNCInterElems, ielem_dataset);
    
    // Number of nodes
    readData.nnod_NF = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNCInterElems_dInterfaceElements_sNfluidNodes, nRulesHDF5::m_cNCInterElems, ielem_dataset);
    readData.nnod_NS = cInputSingletonH5::getInstance()->readIntegerAttributeFromDataset(nRulesHDF5::m_cNCInterElems_dInterfaceElements_sNstructNodes, nRulesHDF5::m_cNCInterElems, ielem_dataset);
}
