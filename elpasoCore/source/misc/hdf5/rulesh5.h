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

//! @brief Globally accessible HDF5 rules used in elPaSo
//! @author Harikrishnan Sreekumar
//! @date 10.08.2020
//!
//! cRulesHDF5 writes data to a h5 file.

#include <string>

namespace nRulesHDF5 {

// Analysis container
static const std::string m_cAnalysis = "/Analysis";
static const std::string m_cAnalysis_sDescription = "description";
static const std::string m_cAnalysis_sType = "type";
static const std::string m_cAnalysis_sSolver = "solver";
static const std::string m_cAnalysis_sRevision = "revision";
static const std::string m_cAnalysis_sStart = "start";
static const std::string m_cAnalysis_sSteps = "steps";
static const std::string m_cAnalysis_sDelta = "delta";
static const std::string m_cAnalysis_sAlpha = "alpha";
static const std::string m_cAnalysis_sAlpham = "alpham";
static const std::string m_cAnalysis_sBeta = "beta";
static const std::string m_cAnalysis_sGamma = "gamma";
static const std::string m_cAnalysis_sRayA = "raya";
static const std::string m_cAnalysis_sRayB = "rayb";

static const std::string m_cAnalysis_sFactor = "factor";
static const std::string m_cAnalysis_sVectors = "vectors";
static const std::string m_cAnalysis_sSweBem = "swebem";
static const std::string m_cAnalysis_sSbFem = "sbfem";
static const std::string m_cAnalysis_sSimilarv4 = "similarv4";
static const std::string m_cAnalysis_sTBEM = "tbem";
static const std::string m_cAnalysis_sStepGeo = "step";

static const std::string m_cAnalysis_sRR0 = "RR0";
static const std::string m_cAnalysis_sER = "ER";
static const std::string m_cAnalysis_sMAXRR = "MAXRR";
static const std::string m_cAnalysis_sOptType = "optType";
static const std::string m_cAnalysis_sStressType = "stressType";

static const std::string m_cAnalysis_sMorMethod = "method";
static const std::string m_cAnalysis_sMorSystem = "system";
static const std::string m_cAnalysis_sMorSetting = "setting";
static const std::string m_cAnalysis_sMorEstimator = "estimator";
static const std::string m_cAnalysis_sMorFreqStart = "mor_freq_start";
static const std::string m_cAnalysis_sMorFreqEnd = "mor_freq_end";
static const std::string m_cAnalysis_sMorFreqSteps = "mor_freq_steps";
static const std::string m_cAnalysis_sMorMaxOrder = "max_order";
static const std::string m_cAnalysis_sMorSigmaTol = "sigma_tol";
static const std::string m_cAnalysis_sMorErrorTol = "error_tol";
static const std::string m_cAnalysis_sMorLoadingPoint =
    "mor_loading_point_loading";
static const std::string m_cAnalysis_sMorLoadingConst =
    "mor_loading_constant_loading";
static const std::string m_cAnalysis_sMorLoadingFreq =
    "mor_loading_freqdep_loading";
static const std::string m_cAnalysis_sMorErrorQuantity = "quantity";
static const std::string m_cAnalysis_vMorActiveDofs = "mor_active_dofs";

static const std::string m_cAnalysis_sPmorNumFom = "num_foms";
static const std::string m_cAnalysis_sPmorFom = "FOM_File";
static const std::string m_cAnalysis_sPmorNumRom = "num_roms";
static const std::string m_cAnalysis_sPmorRom = "ROM_File";
// Analysis output
static const std::string m_cAnalysis_sOutSTP = "stp";
static const std::string m_cAnalysis_sOutVTK = "vtk";
static const std::string m_cAnalysis_sOutHDF5 = "hdf5";

// Node container
static const std::string m_cNodes = "/Nodes";
static const std::string m_cNodes_dFemNodes = "/mtxFemNodes";

// Elements container
static const std::string m_cElements = "/Elements";
static const std::string m_cElements_dFemElemGroup = "/mtxFemElemGroup";
static const std::string m_cElements_dFemElemGroup_sN = "N";
static const std::string m_cElements_dFemElemGroup_sId = "Id";
static const std::string m_cElements_dFemElemGroup_sName = "Name";
static const std::string m_cElements_dFemElemGroup_sMaterial = "MaterialId";
static const std::string m_cElements_dFemElemGroup_sElementType = "ElementType";
static const std::string m_cElements_dFemElemGroup_sOrientation = "Orientation";
static const std::string m_cElements_dFemElemGroup_sOrientationFile =
    "OrientationFile";

// Materials container
static const std::string m_cMaterials = "/Materials";
static const std::string m_cMaterials_sN = "N";
static const std::string m_cMaterials_dMaterial = "/material";
static const std::string m_cMaterials_dMaterial_sMaterialType = "MaterialType";
static const std::string m_cMaterials_dMaterial_sId = "Id";
static const std::string m_cMaterials_dMaterial_sName = "Name";
static const std::string m_cMaterials_dMaterial_sPointmass = "M";
static const std::string m_cMaterials_dMaterial_sE = "E";
static const std::string m_cMaterials_dMaterial_sEx = "Ex";
static const std::string m_cMaterials_dMaterial_sEy = "Ey";
static const std::string m_cMaterials_dMaterial_sEz = "Ez";
static const std::string m_cMaterials_dMaterial_sExMem = "ExMem";
static const std::string m_cMaterials_dMaterial_sEyMem = "EyMem";
static const std::string m_cMaterials_dMaterial_sNu = "nu";
static const std::string m_cMaterials_dMaterial_sNuxy = "nuxy";
static const std::string m_cMaterials_dMaterial_sNuxz = "nuxz";
static const std::string m_cMaterials_dMaterial_sNuyz = "nuyz";
static const std::string m_cMaterials_dMaterial_sNuxyMem = "nuxyMem";
static const std::string m_cMaterials_dMaterial_sGxy = "Gxy";
static const std::string m_cMaterials_dMaterial_sGxz = "Gxz";
static const std::string m_cMaterials_dMaterial_sGyz = "Gyz";
static const std::string m_cMaterials_dMaterial_sGxyMem = "GxyMem";
static const std::string m_cMaterials_dMaterial_sA = "A";
static const std::string m_cMaterials_dMaterial_sI = "I";
static const std::string m_cMaterials_dMaterial_sIx = "Ix";
static const std::string m_cMaterials_dMaterial_sIy = "Iy";
static const std::string m_cMaterials_dMaterial_sIz = "Iz";
static const std::string m_cMaterials_dMaterial_sRho = "rho";
static const std::string m_cMaterials_dMaterial_sRhof = "rhof";
static const std::string m_cMaterials_dMaterial_sC = "c";
static const std::string m_cMaterials_dMaterial_sCf = "cf";
static const std::string m_cMaterials_dMaterial_sT = "t";
static const std::string m_cMaterials_dMaterial_sFi = "Fi";
static const std::string m_cMaterials_dMaterial_spreX = "preX";
static const std::string m_cMaterials_dMaterial_spreY = "preY";
static const std::string m_cMaterials_dMaterial_sTyp = "type";
static const std::string m_cMaterials_dMaterial_sEta = "eta";
static const std::string m_cMaterials_dMaterial_sG = "G";
static const std::string m_cMaterials_dMaterial_sAlpha1 = "alpha1";
static const std::string m_cMaterials_dMaterial_sAlpha = "alpha";
static const std::string m_cMaterials_dMaterial_sAlinf = "alinf";
static const std::string m_cMaterials_dMaterial_sKs = "Ks";
static const std::string m_cMaterials_dMaterial_sKf = "Kf";
static const std::string m_cMaterials_dMaterial_sKappa = "kappa";
static const std::string m_cMaterials_dMaterial_sPhi = "phi";
static const std::string m_cMaterials_dMaterial_sR = "R";
static const std::string m_cMaterials_dMaterial_sKk = "Kk";
static const std::string m_cMaterials_dMaterial_sK = "K";
static const std::string m_cMaterials_dMaterial_sks = "ks";
static const std::string m_cMaterials_dMaterial_sRL = "RL";
static const std::string m_cMaterials_dMaterial_sSt = "st";
static const std::string m_cMaterials_dMaterial_sSv = "sv";
static const std::string m_cMaterials_dMaterial_sAt = "at";
static const std::string m_cMaterials_dMaterial_sSigma = "sigma";
static const std::string m_cMaterials_dMaterial_sLambda = "lambda";
static const std::string m_cMaterials_dMaterial_sLambdat = "lambdat";
static const std::string m_cMaterials_dMaterial_sFactorE = "factorE";
static const std::string m_cMaterials_dMaterial_sPhicm = "phi_cm";
static const std::string m_cMaterials_dMaterial_sCcm = "c_cm";
static const std::string m_cMaterials_dMaterial_sDcm = "D_cm";
static const std::string m_cMaterials_dMaterial_sWcm = "W_cm";
static const std::string m_cMaterials_dMaterial_sXcm = "X_cm";
static const std::string m_cMaterials_dMaterial_sRicl = "R_i_cl";
static const std::string m_cMaterials_dMaterial_sRocl = "R_o_cl";
static const std::string m_cMaterials_dMaterial_sCSxcl = "CSx_cl";
static const std::string m_cMaterials_dMaterial_sCSycl = "CSy_cl";
static const std::string m_cMaterials_dMaterial_sCSzcl = "CSz_cl";
static const std::string m_cMaterials_dMaterial_sCx = "Cx";
static const std::string m_cMaterials_dMaterial_sCy = "Cy";
static const std::string m_cMaterials_dMaterial_sCz = "Cz";
static const std::string m_cMaterials_dMaterial_sCrx = "Crx";
static const std::string m_cMaterials_dMaterial_sCry = "Cry";
static const std::string m_cMaterials_dMaterial_sCrz = "Crz";
static const std::string m_cMaterials_dMaterial_sCreal = "creal";
static const std::string m_cMaterials_dMaterial_sCimag = "cimag";
static const std::string m_cMaterials_dMaterial_sRhoreal = "rhoreal";
static const std::string m_cMaterials_dMaterial_sRhoimag = "rhoimag";
static const std::string m_cMaterials_dMaterial_sE0 = "E0";
static const std::string m_cMaterials_dMaterial_sE1 = "E1";
static const std::string m_cMaterials_dMaterial_sE2 = "E2";
static const std::string m_cMaterials_dMaterial_sE3 = "E3";
static const std::string m_cMaterials_dMaterial_sEta0 = "Eta0";
static const std::string m_cMaterials_dMaterial_sEta1 = "eta1";
static const std::string m_cMaterials_dMaterial_sEta2 = "eta2";
static const std::string m_cMaterials_dMaterial_sEta3 = "eta3";
static const std::string m_cMaterials_dMaterial_sBeta = "beta";
static const std::string m_cMaterials_dMaterial_sB = "B";
static const std::string m_cMaterials_dMaterial_sBtyp = "Btyp";
static const std::string m_cMaterials_dMaterial_sNu1 = "nu1";
static const std::string m_cMaterials_dMaterial_sNu2 = "nu2";
static const std::string m_cMaterials_dMaterial_sNu3 = "nu3";
static const std::string m_cMaterials_dMaterial_sG2 = "G2";
static const std::string m_cMaterials_dMaterial_sH1 = "H1";
static const std::string m_cMaterials_dMaterial_sH2 = "H2";
static const std::string m_cMaterials_dMaterial_sH3 = "H3";
static const std::string m_cMaterials_dMaterial_sRho1 = "rho1";
static const std::string m_cMaterials_dMaterial_sRho2 = "rho2";
static const std::string m_cMaterials_dMaterial_sRho3 = "rho3";
static const std::string m_cMaterials_dMaterial_sWERKU = "wE_RKU";
static const std::string m_cMaterials_dMaterial_sWElam = "wE_Laminate";
static const std::string m_cMaterials_dMaterial_sWNuRKU = "wNu_RKU";
static const std::string m_cMaterials_dMaterial_sWNuLam = "wNU_Laminate";

// NodeLoads container
static const std::string m_cNodeLoads = "/NodeLoads";
static const std::string m_cNodeLoads_dFemNodeLoad = "/mtxFemNodeLoad";
static const std::string m_cNodeLoads_dFemNodeLoad_sId = "Id";
static const std::string m_cNodeLoads_dFemNodeLoad_sLoadType = "LoadType";
static const std::string m_cNodeLoads_dFemNodeLoad_sMethodType = "MethodType";
static const std::string m_cNodeLoads_dFemNodeLoad_sNodeId = "NodeId";
static const std::string m_cNodeLoads_dFemNodeLoad_sFx = "Fx";
static const std::string m_cNodeLoads_dFemNodeLoad_sFy = "Fy";
static const std::string m_cNodeLoads_dFemNodeLoad_sFz = "Fz";
static const std::string m_cNodeLoads_dFemNodeLoad_sP = "P";
static const std::string m_cNodeLoads_dFemNodeLoad_sSteps = "Steps";
static const std::string m_cNodeLoads_dFemNodeLoad_sDeltaT = "DeltaT";

// NodeMoments container
static const std::string m_cNodeMoments = "/NodeMoments";
static const std::string m_cNodeMoments_dFemNodeMoment = "/mtxFemNodeMoment";
static const std::string m_cNodeMoments_dFemNodeMoment_sId = "Id";
static const std::string m_cNodeMoments_dFemNodeMoment_sLoadType = "LoadType";
static const std::string m_cNodeMoments_dFemNodeMoment_sMethodType =
    "MethodType";
static const std::string m_cNodeMoments_dFemNodeMoment_sNodeId = "NodeId";
static const std::string m_cNodeMoments_dFemNodeMoment_snSteps = "nSteps";
static const std::string m_cNodeMoments_dFemNodeMoment_sDeltaT = "DeltaT";
static const std::string m_cNodeMoments_dFemNodeMoment_sMx = "Mx";
static const std::string m_cNodeMoments_dFemNodeMoment_sMy = "My";
static const std::string m_cNodeMoments_dFemNodeMoment_sMz = "Mz";

// ElemLoads container
static const std::string m_cElemLoads = "/ElemLoads";
static const std::string m_cElemLoads_dFemElemLoad = "/mtxFemElemLoad";
static const std::string m_cElemLoads_dFemElemLoad_sElementId = "ElementId";
static const std::string m_cElemLoads_dFemElemLoad_sFreqCount = "FreqCount";
static const std::string m_cElemLoads_dFemElemLoad_sId = "Id";
static const std::string m_cElemLoads_dFemElemLoad_sLoadType = "LoadType";
static const std::string m_cElemLoads_dFemElemLoad_sMethodType = "MethodType";
static const std::string m_cElemLoads_dFemElemLoad_sVn = "vn";
static const std::string m_cElemLoads_dFemElemLoad_sFace = "Face";

// InterfaceElements container
static const std::string m_cInterElems = "/InterfaceElements";
static const std::string m_cInterElems_sN = "N";
static const std::string m_cInterElems_dInterfaceElements_sNfluidNodes =
    "NfluidNodes";
static const std::string m_cInterElems_dInterfaceElements_sNstructNodes =
    "NstructNodes";
static const std::string m_cInterElems_dInterfaceElements_sMatF = "matF";
static const std::string m_cInterElems_dInterfaceElements_sMatS = "mats";
static const std::string m_cInterElems_dInterfaceElements_sOri = "ori";
static const std::string m_cInterElems_dInterfaceElements_sType = "type";

// Non conforming InterfaceElements container
static const std::string m_cNCInterElems = "/NCInterfaceElements";
static const std::string m_cNCInterElems_sN = "N";
static const std::string m_cNCInterElems_dInterfaceElements_sNfluidNodes =
    "NfluidNodes";
static const std::string m_cNCInterElems_dInterfaceElements_sNstructNodes =
    "NstructNodes";
static const std::string m_cNCInterElems_dInterfaceElements_sMatF = "matF";
static const std::string m_cNCInterElems_dInterfaceElements_sMatS = "mats";
static const std::string m_cNCInterElems_dInterfaceElements_sOri = "ori";
static const std::string m_cNCInterElems_dInterfaceElements_sType = "type";
static const std::string m_cNCNodes = "/InterNodes";
static const std::string m_cNCNodes_dNCNodes = "/mtxFemInterNodes";

// NodeConstraints container
static const std::string m_cNodeConstraints = "/NodeConstraints";
static const std::string m_cNodeConstraints_sNodeConstraint = "nodeConstraint";
static const std::string m_cNodeConstraints_sNodeConstraint_sConType =
    "ConstraintType";
static const std::string m_cNodeConstraints_sNodeConstraint_sName = "Name";
static const std::string m_cNodeConstraints_sNodeConstraint_sId = "Id";
static const std::string m_cNodeConstraints_sNodeConstraint_sMethodType =
    "MethodType";
static const std::string m_cNodeConstraints_sNodeConstraint_sNodeId = "NodeId";
static const std::string m_cNodeConstraints_sNodeConstraint_sU1 = "u1";
static const std::string m_cNodeConstraints_sNodeConstraint_sU2 = "u2";
static const std::string m_cNodeConstraints_sNodeConstraint_sU3 = "u3";
static const std::string m_cNodeConstraints_sNodeConstraint_sW1 = "w1";
static const std::string m_cNodeConstraints_sNodeConstraint_sW2 = "w2";
static const std::string m_cNodeConstraints_sNodeConstraint_sW3 = "w3";

static const std::string m_cNodeConstraints_sNodeConstraint_sP0 = "p0";
static const std::string m_cNodeConstraints_sNodeConstraint_sP1 = "p1";
static const std::string m_cNodeConstraints_sNodeConstraint_sP2 = "p2";
static const std::string m_cNodeConstraints_sNodeConstraint_sP3 = "p3";
static const std::string m_cNodeConstraints_sNodeConstraint_sXD3 = "xd3";
static const std::string m_cNodeConstraints_sNodeConstraint_sWD1 = "wd1";
static const std::string m_cNodeConstraints_sNodeConstraint_sWD2 = "wd2";
static const std::string m_cNodeConstraints_sNodeConstraint_sdwdx = "dwdx";
static const std::string m_cNodeConstraints_sNodeConstraint_sdwdy = "dwdy";
static const std::string m_cNodeConstraints_sNodeConstraint_sdwdxy = "dwdxy";
static const std::string m_cNodeConstraints_sNodeConstraint_sFluid = "fluid";
static const std::string m_cNodeConstraints_sNodeConstraint_sZ1 = "z_1";
static const std::string m_cNodeConstraints_sNodeConstraint_sZ3 = "z_3";
static const std::string m_cNodeConstraints_sNodeConstraint_sX12 = "x1_2";
static const std::string m_cNodeConstraints_sNodeConstraint_sX22 = "x2_2";

static const std::string m_cNodeConstraints_sNodeConstraint_sValU1 = "valu1";
static const std::string m_cNodeConstraints_sNodeConstraint_sValU2 = "valu2";
static const std::string m_cNodeConstraints_sNodeConstraint_sValU3 = "valu3";
static const std::string m_cNodeConstraints_sNodeConstraint_sValW1 = "valw1";
static const std::string m_cNodeConstraints_sNodeConstraint_sValW2 = "valw2";
static const std::string m_cNodeConstraints_sNodeConstraint_sValW3 = "valw3";
static const std::string m_cNodeConstraints_sNodeConstraint_sValP0 = "valp0";
static const std::string m_cNodeConstraints_sNodeConstraint_sValP1 = "valp1";
static const std::string m_cNodeConstraints_sNodeConstraint_sValP2 = "valp2";
static const std::string m_cNodeConstraints_sNodeConstraint_sValP3 = "valp3";
static const std::string m_cNodeConstraints_sNodeConstraint_sValXD3 = "valxd3";
static const std::string m_cNodeConstraints_sNodeConstraint_sValWD1 = "valwd1";
static const std::string m_cNodeConstraints_sNodeConstraint_sValWD2 = "valwd2";
static const std::string m_cNodeConstraints_sNodeConstraint_sValdwdx =
    "valdwdx";
static const std::string m_cNodeConstraints_sNodeConstraint_sValdwdy =
    "valdwdy";
static const std::string m_cNodeConstraints_sNodeConstraint_sValdwdxy =
    "valdwdxy";
static const std::string m_cNodeConstraints_sNodeConstraint_sValFluid =
    "valfluid";
static const std::string m_cNodeConstraints_sNodeConstraint_sValZ1 = "valz_1";
static const std::string m_cNodeConstraints_sNodeConstraint_sValZ3 = "valz_3";
static const std::string m_cNodeConstraints_sNodeConstraint_sValX12 = "valx1_2";
static const std::string m_cNodeConstraints_sNodeConstraint_sValX22 = "valx2_2";

static const std::string m_cNodeConstraints_sNodeConstraint_sValP = "p";
static const std::string m_cNodeConstraints_sNodeConstraint_sValIx = "Ix";
static const std::string m_cNodeConstraints_sNodeConstraint_sValIy = "Iy";
static const std::string m_cNodeConstraints_sNodeConstraint_sValIz = "Iz";
static const std::string m_cNodeConstraints_sNodeConstraint_sValCs = "Cs";

// Nodesets
static const std::string m_cNodesets = "/Nodesets";
static const std::string m_cNodesets_sMorInput = "/vecMorInput";
static const std::string m_cNodesets_sMorOutput = "/vecMorOutput";
}  // namespace nRulesHDF5