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

#include "femparserinterface.h"
#include "../../fedata/mesh.h"

cFemParserInterface::cFemParserInterface() {
    mEstimatedNumberOfEntities.AttributesExistInInputFile = false;
    mEstimatedNumberOfEntities.Interfaces = 0;
    mEstimatedNumberOfEntities.Nodes = 0;
    mEstimatedNumberOfEntities.Elements = 0;
    mEstimatedNumberOfEntities.Materials = 0;
    mEstimatedNumberOfEntities.NodalForces = 0;
    mEstimatedNumberOfEntities.NodalMoments = 0;
    mEstimatedNumberOfEntities.NodalPressures = 0;
    mEstimatedNumberOfEntities.NodalValues = 0;
    mEstimatedNumberOfEntities.LoadedNodes = 0;
    mEstimatedNumberOfEntities.MomentedNodes = 0;
    mEstimatedNumberOfEntities.PressuredNodes = 0;
    mEstimatedNumberOfEntities.NodeBCs = 0;
    mEstimatedNumberOfEntities.FixedNodes = 0;
    mEstimatedNumberOfEntities.ElementLoads = 0;
    mEstimatedNumberOfEntities.LoadedElements = 0;
    mEstimatedNumberOfEntities.InterfaceElements = 0;
    mEstimatedNumberOfEntities.ElementsFF = 0;
    mEstimatedNumberOfEntities.Constraints = 0;
}

cFemParserInterface::~cFemParserInterface() {
    // empty
}

/**
* Newer version have additional attributes at the top level of the different
* sections where the number of expected entries is stored. They are compared
* here to the number of read objects.
*/
void cFemParserInterface::checkCounters(cMesh* ptrMesh)
{
    bool ErrorFlag = false;

    trace("==============================");
    message("  Revision of Input-File : %d\n", m_Revision);
    trace("==============================");

    trace("checking number of read objects .... ");

    if (mEstimatedNumberOfEntities.AttributesExistInInputFile == false)
    {
        trace("  old input file format => not able to check counters!");
        trace("  Therefore, aborting test.");
        return;
    }

    if (ptrMesh->getNumberOfNodes() == mEstimatedNumberOfEntities.Nodes)
        trace("  Number of nodes   ** ok **");
    else
    {
        trace("  *** Number of nodes differs ***");
        message("      expected : %7d\n", mEstimatedNumberOfEntities.Nodes);
        message("      read     : %7d\n", ptrMesh->getNumberOfNodes());
        ErrorFlag = true;
    }

    if (ptrMesh->getNumberOfElements() == mEstimatedNumberOfEntities.Elements)
        trace("  Number of elements   ** ok **");
    else
    {
        trace("  *** Number of elements differs ***");
        message("      expecting : %7d\n", mEstimatedNumberOfEntities.Elements);
        message("      read      : %7d\n", ptrMesh->getNumberOfElements());
        ErrorFlag = true;
    }

    if (ptrMesh->getNumberOfMaterials() == mEstimatedNumberOfEntities.Materials)
        trace("  Number of materials   ** ok **");
    else
    {
        trace("  *** Number of materials differs ***");
        message("      expected : %7d\n", mEstimatedNumberOfEntities.Materials);
        message("      read     : %7d\n", ptrMesh->getNumberOfMaterials());
        ErrorFlag = true;
    }

    const int CountInterface = ptrMesh->getNumberOfInterfaceElements()
        + ptrMesh->getNumberOfCplFemBemElements();

    if (CountInterface == mEstimatedNumberOfEntities.InterfaceElements)
        trace("  Number of interface elements   ** ok **");
    else
    {
        trace("  *** Number of interface elements differs ***");
        message("      expected : %7d\n", mEstimatedNumberOfEntities.InterfaceElements);
        message("      read     : %7d\n", ptrMesh->getNumberOfInterfaceElements());
        ErrorFlag = true;
    }

    if (ptrMesh->getNumberOfBoundaryConditions() == mEstimatedNumberOfEntities.NodeBCs)
        trace("  Number of nodal boundary conditions    ** ok **");
    else
    {
        trace("  *** Number of nodal boundary conditions differs ***");
        message("      expected : %7d\n", mEstimatedNumberOfEntities.NodeBCs);
        message("      read     : %7d\n", ptrMesh->getNumberOfBoundaryConditions());
        ErrorFlag = true;
    }

    if (ptrMesh->getNumberOfNodalForces() == mEstimatedNumberOfEntities.NodalForces)
        trace("  Number of nodal forces    ** ok **");
    else
    {
        trace("  *** Number of nodal forces differs ***");
        message("      expected : %7d\n", mEstimatedNumberOfEntities.NodalForces);
        message("      read     : %7d\n", ptrMesh->getNumberOfNodalForces());
        ErrorFlag = true;
    }

    if (ptrMesh->getNumberOfElementLoads() == mEstimatedNumberOfEntities.ElementLoads)
        trace("  Number of element loads   ** ok **");
    else
    {
        trace("  *** Number of element loads differs ***");
        message("      expected : %7d\n", mEstimatedNumberOfEntities.ElementLoads);
        message("      read     : %7d\n", ptrMesh->getNumberOfElementLoads());
        ErrorFlag = true;
    }

    if (ptrMesh->getNumberOfElementsFF() == mEstimatedNumberOfEntities.ElementsFF)
        trace("  Number of fluid flow elements   ** ok **");
    else
    {
        trace("  *** Number of fluid flow elements differs ***");
        message("      expected : %7d\n", mEstimatedNumberOfEntities.ElementsFF);
        message("      read     : %7d\n", ptrMesh->getNumberOfElementsFF());
        ErrorFlag = true;
    }

    if (ptrMesh->getNumberOfConstraints() == mEstimatedNumberOfEntities.Constraints)
        trace("  Number of constraint    ** ok **");
    else
    {
        trace("  *** Number of constraint differs ***");
        message("      expected : %7d\n", mEstimatedNumberOfEntities.Constraints);
        message("      read     : %7d\n", ptrMesh->getNumberOfConstraints());
        ErrorFlag = true;
    }

    // --- terminate the code if there are no loads applied to the system
    //   if ((mEstimatedNumberOfEntities.ElementLoads == 0) && (mEstimatedNumberOfEntities.NodalForces == 0))
    //   {
    //     logstream << " ************************ ");
    //     logstream << " *** NO LOADS APPLIED *** ");
    //     logstream << " ************************ ");
    //     ErrorFlag = true;
    //   }


    // --- check if symmetry is possible
    //     symmetry does not work for:
    //        FEM-BEM coupling
    //        2d poro
    PetscBool optSymmetric, symError = PETSC_FALSE;
    PetscOptionsHasName(PETSC_NULL, PETSC_NULL, "-symmetric", &optSymmetric);

    if (ErrorFlag == true)
    {
        throw cException("Error in inputfile detected", __FILE__, __LINE__);
    }
    else
        trace("No error found in inputfile !");
}

void cFemParserInterface::parseInputFile(cProblem& myProblem)
{
    try
    {
        message("Inputfile : %s\n", myProblem.getFilename().c_str());
        this->openInputFile(myProblem.getFilename());

        // Entity | Description
        this->parseDescriptionEntity(myProblem);
        // Entity | Revision
        this->parseRevisionEntity();
        // Entity | Analysis
        this->parseAnalysisEntity(myProblem);
        // Entity | Nodes
        this->parseNodeEntity(myProblem);
        // Entity | Materials
        this->parseMaterialsEntity(myProblem);
        // Entity | Elements
        this->parseElementsEntity(myProblem);
        // Entity | ElementLoads
        this->parseElementLoadsEntity(myProblem);
        // Entity | NodeConstraints
        this->parseNodeConstraintsEntity(myProblem);
        // Entity | NodeLoads
        this->parseNodeLoadsEntity(myProblem);
        // Entity | cNodeMoments
        this->parseNodeMomentsEntity(myProblem);
        // Entity | InterfaceElements
        this->parseInterfaceElementsEntity(myProblem);
        // Entity | NCInterfaceElements
        this->parseNCInterfaceElementsEntity(myProblem);

        dumpObjectCounters(logfile); // write tCounter<T>-objects to logfile
        checkCounters(myProblem.getMesh());
    }
    catch (const cException& err)
    {
        PetscPrintf(PETSC_COMM_SELF, "%s\n", err.what().c_str());
        ExitApp();
    }
    catch (...)
    {
        PetscPrintf(PETSC_COMM_SELF, "unknown exceptione (parseInputFile)\n");
        ExitApp();
    }
}

void cFemParserInterface::parseAnalysisEntity(cProblem& userData)
{
    trace("  reading analysis parameters ...");
    std::string     Type;
    cAnalysis*      ptr = 0;

    Type = this->getAnalysisType();

    cMesh* MyMesh = userData.getMesh();
    ptr = m_AnalysisFactory->createAnalysis(Type, this, *MyMesh);
    if(ptr == 0)
    {
        std::cout<<"analysis "<< Type << " not supported"<< std::endl;
        ExitApp();
    }
    userData.setAnalysis(ptr);
    ptr = 0;
    
    trace("  reading output parameters ...");
    cAnalysis* ptrAnalysis = userData.getAnalysis();
    m_AnalysisFactory->prepareAnalysisOutputs(ptrAnalysis, this);

    ptrAnalysis = 0;
}

void cFemParserInterface::parseNodeEntity(cProblem& userData)
{
    trace("  reading nodes ...");
    mEstimatedNumberOfEntities.AttributesExistInInputFile = true;

    ParserNodeData readData = this->getNodalData();
    int _m = readData.nodalIds.size();
    mEstimatedNumberOfEntities.Nodes = _m;
    message("    expecting %d nodes ...\n", mEstimatedNumberOfEntities.Nodes);

    for (size_t i = 0; i < _m; i++) // for all nodes
    {
        cNode* ptr = new cNode;
        cMesh* pMesh = userData.getMesh();
        ptr->setId(readData.nodalIds[i]);
        (*ptr)[0] = readData.nodalCoordinates[i * 3 + 0];
        (*ptr)[1] = readData.nodalCoordinates[i * 3 + 1];
        (*ptr)[2] = readData.nodalCoordinates[i * 3 + 2];
        pMesh->insertNode(ptr);
        pMesh = 0;
        ptr = 0;
        delete pMesh;
        delete ptr;
    }
}

void cFemParserInterface::parseMaterialsEntity(cProblem& userData)
{
    trace("  reading material parameters ...");

    mEstimatedNumberOfEntities.AttributesExistInInputFile = true;
    mEstimatedNumberOfEntities.Materials = this->getNumberOfMaterials();

    message("    expecting %d materials ...\n", mEstimatedNumberOfEntities.Materials);

    for (size_t i = 1; i <= mEstimatedNumberOfEntities.Materials; i++)
    {
        ParserMaterialData readData = this->getMaterialData(i);

        cMaterial* pMaterial = 0;
        pMaterial = m_MaterialFactory->createMaterial(readData);
        if(pMaterial==0)
        {
            std::cout<<"material "<< readData.matType << " not supported"<< std::endl;
            ExitApp();
        }
        
        userData.getMesh()->insertMaterial(pMaterial);
        pMaterial = 0;
    }
}

void cFemParserInterface::parseElementsEntity(cProblem& userData)
{
    int NumElementsInThisGroup = 0;
    int NumElemGroups = this->getNumberOfElementBlocks();
    
    trace("  reading elements ...");
    for (size_t iElemGroup = 0; iElemGroup < NumElemGroups; iElemGroup++)
    {
        ParserElementBlockData readData = this->getElementBlockData(iElemGroup);
        mEstimatedNumberOfEntities.AttributesExistInInputFile = true;
        NumElementsInThisGroup = readData.blockNumElements;
        mEstimatedNumberOfEntities.Elements += NumElementsInThisGroup;

        message("    expecting %d elements for this group\n", NumElementsInThisGroup);

        PetscInt GroupId;
        std::string NameOfGroup = readData.blockName;
        m_SectionsOrientation = readData.blockOrientation;
        m_SectionsOrientationFilename = readData.blockOrientationFilename;

        GroupId = readData.blockId;
        m_SectionsMaterial = readData.blockMaterial;

        message("    identifier of this group .....: %s\n", NameOfGroup.c_str());
        message("    id of this section's material : %d\n", m_SectionsMaterial);
        message("    id of this section's group ...: %d\n", GroupId);
        
        //check, if group is a fluid flow group (elements will not be FEM elements
        if (NameOfGroup.compare(0, 3, "ff_") == 0)
        {
            mEstimatedNumberOfEntities.Elements -= NumElementsInThisGroup;
            mEstimatedNumberOfEntities.ElementsFF += NumElementsInThisGroup;
            m_SectionsMaterial = 0; // Material set to zero as it does not apply to fluid flow elems
        }

        cGroup* ptrGroup = new cGroup;
        ptrGroup->setIdentifier(NameOfGroup);
        ptrGroup->setId(GroupId);
        ptrGroup->setMaterialId(m_SectionsMaterial);
        ptrGroup->setOrientation(m_SectionsOrientation);
        ptrGroup->setOrientationFilename(m_SectionsOrientationFilename);
        if (m_SectionsOrientationFilename != "") // Only read the orientation vectors if a filename is given
        {
            readOrientationFromFile();
        }
        if (m_SectionsOrientation == "user-def" && m_SectionsOrientationFilename == "")
        {
            PetscPrintf(PETSC_COMM_SELF, "Orientation type is user-def but no orientation filename is given!\n");
            throw cException("No orientation filename given", __FILE__, __LINE__);
        }
        userData.getMesh()->insertGroup(ptrGroup);
        userData.getMesh()->makeCurrentGroup(GroupId);
        ptrGroup = NULL;

        // Find the element type and handle the corresponding element handler
        std::string ElementType = readData.blockElementType;

        for (size_t iElement = 0; iElement < NumElementsInThisGroup; iElement++)
        {
            int _m = readData.blockConnectiviyNumRows, _n = readData.blockConnectiviyNumCols;
            cElementFEM* ptrElement = m_ElementFactory->createElement(ElementType);
            if(ptrElement==0)
            {
                std::cout<<"element "<< ElementType << " not supported"<< std::endl;
                ExitApp();
            }

            int ElemId = (int)readData.blockConnectiviy[iElement * _n + 0];
            // -------------------------------------------------------------------------
            //   whole data read - now setting up element
            // -------------------------------------------------------------------------
            ptrElement->setId(ElemId);

            // materials
            int MatId = m_SectionsMaterial;

            if (MatId == 0)
            {
                PetscPrintf(PETSC_COMM_SELF, "id of the element : %d\n", ptrElement->getId());
                throw cException("no material defined for element", __FILE__, __LINE__);
            }
            else if (MatId == -1)
            {
                PetscPrintf(PETSC_COMM_SELF, "id of the element : %d\n", ptrElement->getId());
                PetscPrintf(PETSC_COMM_SELF, "Did you specify a material for this group of elements?");
                throw cException("no material defined for element", __FILE__, __LINE__);
            }
            else
            {
                cMaterial* ptr = NULL;
                ptr = userData.getMesh()->getMaterial(MatId);

                if (ptr == NULL)
                {
                    PetscPrintf(PETSC_COMM_SELF, " Element Id = %d, : assigned material %d not found\n", ElemId, MatId);
                    throw cException("material not found for element", __FILE__, __LINE__);
                }
                else
                {
                    ptrElement->setMaterial(ptr);
                    ptrElement->setOrientation(m_SectionsOrientation);
                }

                // nodes
                for (int i = 0; i < ptrElement->getNumberOfNodes(); i++)
                    ptrElement->setNode(i, userData.getMesh()->getNode(readData.blockConnectiviy[iElement * _n + 1 + i]));

                ptr = NULL;
            }

            this->processElement(ptrElement, userData.getMesh());
            ptrElement = NULL;
        }
    }
}

void cFemParserInterface::parseElementLoadsEntity(cProblem &userData)
{
    trace("  reading elementloads and assigning to elements...");

    mEstimatedNumberOfEntities.AttributesExistInInputFile = true;
    int _m = this->getNumberOfElementLoads();
    mEstimatedNumberOfEntities.ElementLoads = (int)_m;
    mEstimatedNumberOfEntities.LoadedElements = (int)_m;

    message("    expecting %d elementloads ...\n", mEstimatedNumberOfEntities.ElementLoads);

    try {
        for (size_t i = 0; i < mEstimatedNumberOfEntities.ElementLoads; i++)
        {
            // - Identify the elemload type and set the definition of load
            ParserElementLoadsData readData = getElementLoadsData(i);

            cElementLoad* ptrElemLoad = m_ElementLoadFactory->createElementLoad(readData);
            if(ptrElemLoad==0)
            {
                std::cout << "elementload " << readData.eloadType << " not supported" << std::endl;
                ExitApp();
            }

            userData.getMesh()->insertElementLoad(ptrElemLoad);
            ptrElemLoad = 0;

            // assign the cooresponding element with the load
            int ElementId = readData.eloadElementId;
            cElementFEM* ptrElement = userData.getMesh()->getElement(ElementId);
            ptrElement->insertElementLoad(readData.eloadFace, userData.getMesh()->getElementLoad(Id));
            ptrElement = 0;
        }
    }
    catch (cException& err)
    {
        PetscPrintf(PETSC_COMM_SELF, "%s\n", err.what().c_str());
        ExitApp();
    }
}

void cFemParserInterface::processElement(cElementFEM* ptrElement, cMesh* ptrMesh)
{
    // Search for the correct material orientation according to the element's location if a user-defined file is given
    if (m_SectionsOrientation == "user-def")
    {
        // Center of element
        std::vector<PetscReal> element_center(3);
        PetscInt element_nodes = ptrElement->getNumberOfNodes();

        element_center[0] = 0.;
        element_center[1] = 0.;
        element_center[2] = 0.;
        for (int i = 0; i < element_nodes; i++)
        {
            element_center[0] += (ptrElement->getNode(i)->getComponent(0)) / float(element_nodes);
            element_center[1] += (ptrElement->getNode(i)->getComponent(1)) / float(element_nodes);
            element_center[2] += (ptrElement->getNode(i)->getComponent(2)) / float(element_nodes);
        }

        // Check distance to all given vectors and remember the nearest vector to that element
        PetscReal last_distance = 1000000.;
        PetscReal current_distance = 0.;
        PetscInt index_of_nearest_vector = 0;
        for (int i = 0; i < m_lines_in_orient_file; i++)
        {
            current_distance = std::sqrt((m_SectionsOrientationVectors(i, 0) - element_center[0]) * (m_SectionsOrientationVectors(i, 0) - element_center[0]) + (m_SectionsOrientationVectors(i, 1) - element_center[1]) * (m_SectionsOrientationVectors(i, 1) - element_center[1]) + (m_SectionsOrientationVectors(i, 2) - element_center[2]) * (m_SectionsOrientationVectors(i, 2) - element_center[2]));
            if (current_distance < last_distance)
            {
                index_of_nearest_vector = i;
                last_distance = current_distance;
            }
        }

        // Set the material orientation to the nearest vector found
        std::vector<PetscReal> elementOrientation(6);
        elementOrientation[0] = m_SectionsOrientationVectors(index_of_nearest_vector, 3); // Components of x direction
        elementOrientation[1] = m_SectionsOrientationVectors(index_of_nearest_vector, 4);
        elementOrientation[2] = m_SectionsOrientationVectors(index_of_nearest_vector, 5);
        elementOrientation[3] = m_SectionsOrientationVectors(index_of_nearest_vector, 6); // Components of y direction (only used in 3D)
        elementOrientation[4] = m_SectionsOrientationVectors(index_of_nearest_vector, 7);
        elementOrientation[5] = m_SectionsOrientationVectors(index_of_nearest_vector, 8);
        // Norm the vectors
        PetscReal vectorLengthl = std::sqrt(elementOrientation[0] * elementOrientation[0] + elementOrientation[1] * elementOrientation[1] + elementOrientation[2] * elementOrientation[2]);
        PetscReal vectorLengthr = std::sqrt(elementOrientation[3] * elementOrientation[3] + elementOrientation[4] * elementOrientation[4] + elementOrientation[5] * elementOrientation[5]);
        for (int i = 0; i < 3; i++)
        {
            elementOrientation[i] = elementOrientation[i] / vectorLengthl;
            elementOrientation[i + 3] = elementOrientation[i + 3] / vectorLengthr;
        }
        // Save the element orientation
        ptrElement->setOrientationVector(elementOrientation);
    }

    // -------------------------------------------------------------------------
    //   now save the element
    // -------------------------------------------------------------------------
    ptrMesh->insertElement(ptrElement);
    ptrElement = 0;
}

void cFemParserInterface::parseNodeConstraintsEntity(cProblem& userData)
{
    trace("  reading nodal boundary conditions and fixed nodes...");

    mEstimatedNumberOfEntities.AttributesExistInInputFile = true;
    mEstimatedNumberOfEntities.NodeBCs = this->getNumberOfNodeConstraints();
    mEstimatedNumberOfEntities.FixedNodes = mEstimatedNumberOfEntities.NodeBCs;

    message("    expecting %d nodal boundary conditions ...\n", mEstimatedNumberOfEntities.NodeBCs);

    for (size_t i = 0; i < mEstimatedNumberOfEntities.NodeBCs; i++)
    {
        std::string         type;
        cBoundaryCondition* ptrBoundCond = 0;

        type = this->getNodeConstraintsType(i);

        ptrBoundCond = m_BoundaryConditionFactory->createBoundaryCondition(type, i, this);
        if(ptrBoundCond==0)
        {
            std::cout<<"boundary condition "<< type << " not supported"<< std::endl;
            ExitApp();
        }
        userData.getMesh()->insertBoundaryCondition(ptrBoundCond);
        ptrBoundCond = 0;

        ParserNodeConstraintIdsData readIdsData = this->getNodeConstraintsIdsData(i);
        int NodeId, BcId;
        NodeId  = readIdsData.nconstraintNodeId;
        BcId    = readIdsData.nconstraintBcId;

        userData.getMesh()->getNode(NodeId)->insertBoundaryCondition(
            userData.getMesh()->getBoundaryCondition(BcId)
        );
    } 
}

void cFemParserInterface::parseNodeLoadsEntity(cProblem& userData)
{

    trace("  reading nodal forces ...");

    mEstimatedNumberOfEntities.AttributesExistInInputFile = true;

    mEstimatedNumberOfEntities.NodalForces = this->getNumberOfNodeLoads();
    mEstimatedNumberOfEntities.LoadedNodes = mEstimatedNumberOfEntities.NodalForces;
    message("    expecting %d nodal forces ...\n", mEstimatedNumberOfEntities.NodalForces);

    for (size_t i = 0; i < mEstimatedNumberOfEntities.NodalForces; i++)
    {
        ParserNodeLoadsData readData =this->getNodeLoadsData(i);
        cNodalForce* ptr = m_NodalForceFactory->createNodalForce(readData);
        if(ptr==0)
        {
            std::cout << "nodal force "<< readData.nloadType << " not supported" << std::endl;
            ExitApp();
        }
        userData.getMesh()->insertNodalForce(ptr);

        ptr = NULL;

        userData.getMesh()->getNode(readData.nloadNodeId)->setNodalForce(
            userData.getMesh()->getNodalForce(readData.nloadId)
        );
    }
}

void cFemParserInterface::parseNodeMomentsEntity(cProblem& userData)
{
    trace("  reading nodal moments ...");

    mEstimatedNumberOfEntities.AttributesExistInInputFile = true;

    mEstimatedNumberOfEntities.NodalMoments = this->getNumberOfNodeMoments();
    mEstimatedNumberOfEntities.MomentedNodes = mEstimatedNumberOfEntities.NodalMoments;
    message("    expecting %d nodal moments ...\n", mEstimatedNumberOfEntities.NodalMoments);


    for (size_t i = 0; i < mEstimatedNumberOfEntities.NodalMoments; i++)
    {
        ParserNodeMomentsData   readData = this->getNodeMomentsData(i);
        cNodalMoment*           ptr = m_NodalMomentFactory->createNodalMoment(readData);
        if(ptr==0)
        {
            std::cout << "nodal moment "<< readData.nmomentType << " not supported" << std::endl;
            ExitApp();
        }
        userData.getMesh()->insertNodalMoment(ptr);

        ptr = NULL;
        userData.getMesh()->getNode(readData.nmomentNodeId)->setNodalMoment(
            userData.getMesh()->getNodalMoment(readData.nmomentId)
        );
    }
}

void cFemParserInterface::parseInterfaceElementsEntity(cProblem& userData)
{
    trace("  reading interface elements ...");

    mEstimatedNumberOfEntities.AttributesExistInInputFile = true;
    mEstimatedNumberOfEntities.InterfaceElements = this->getNumberOfInterfaceElements();

    message("    expecting %d interface elements ...\n", mEstimatedNumberOfEntities.InterfaceElements);

    int numberOfGroups = this->getNumberOfInterfaces();
    for (size_t i = 0; i < numberOfGroups; i++)
    {
        ParserInterfaceElementsData ielemData = this->getInterfaceElementsData(i);

        int nnod = -1;  // number of nodes in general
        if (ielemData.nnod_NF == ielemData.nnod_NS)
            nnod = ielemData.nnod_NF;

        for (size_t iElemInGroup = 0; iElemInGroup < ielemData.groupdata_m; iElemInGroup++)
        {
            // Define struct
            AppInterfaceElements interElemData;
            interElemData.dataspace = &ielemData.groupdata[iElemInGroup * ielemData.groupdata_n];
            interElemData.dataspace_size = ielemData.groupdata_n;
            interElemData.MatF = ielemData.MatF;
            interElemData.MatS = ielemData.MatS;
            interElemData.nnod = nnod;

            bool found = m_ElementInterfaceFactory->createElementInterface(ielemData.Type, interElemData, &userData);
            if(!found)
            {
                std::cout << "element interface "<< ielemData.Type << " not supported" << std::endl;
                ExitApp();
            }

        }
    }
}

void cFemParserInterface::parseNCInterfaceElementsEntity(cProblem& userData)
{
    trace("  reading non conforming interface elements ...");

    mEstimatedNumberOfEntities.AttributesExistInInputFile = true;
    mEstimatedNumberOfEntities.NCInterfaceElements = this->getNumberOfNCInterfaceElements();

    message("    expecting %d non conforming interface elements ...\n", mEstimatedNumberOfEntities.NCInterfaceElements);

    int numberOfGroups = this->getNumberOfNCInterfaces();

    for (size_t i = 0; i < numberOfGroups; i++)
    {
        ParserNCInterfaceElementsData ielemData = this->getNCInterfaceElementsData(i);
        
        // Reading the virtual nodes created by the intersection elements 
        trace("  reading intersection nodes ...");

        // info on intersection nodes
        std::vector<double> _matvec_coord_elem;

        int nnod = -1;  // number of nodes in general
        if (ielemData.nnod_NF == ielemData.nnod_NS)
            nnod = ielemData.nnod_NF;
  
        for (size_t iElemInGroup = 0; iElemInGroup < ielemData.groupdata_m; iElemInGroup++)
        {
            // Define struct
            AppInterfaceElements interElemData;
            interElemData.dataspace = &ielemData.groupdata[iElemInGroup *  ielemData.groupdata_n];
            interElemData.dataspace_size = ielemData.groupdata_n;
            interElemData.MatF = ielemData.MatF;
            interElemData.MatS = ielemData.MatS;
            interElemData.nnod = nnod;
  
            // The virtual coordinates of the intersection element are also needed
            for (size_t idx = 0; idx < nnod; idx++)
            {
                _matvec_coord_elem.push_back(ielemData.matvec_coord[3 * interElemData.dataspace[2 + idx] - 3]);
                _matvec_coord_elem.push_back(ielemData.matvec_coord[3 * interElemData.dataspace[2 + idx] - 2]);
                _matvec_coord_elem.push_back(ielemData.matvec_coord[3 * interElemData.dataspace[2 + idx] - 1]);
            }
            interElemData._matvec_coord_elem = _matvec_coord_elem;
            _matvec_coord_elem.erase(_matvec_coord_elem.begin(),_matvec_coord_elem.end());
            
            bool found = m_NCElementInterfaceFactory->createNCElementInterface(ielemData.Type, interElemData, &userData);
            if(!found)
            {
                std::cout << "element interface "<< ielemData.Type << " not supported" << std::endl;
                ExitApp();
            }
        }
    }
}