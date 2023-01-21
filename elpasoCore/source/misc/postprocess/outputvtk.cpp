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

#include "outputvtk.h"
#include "../../basics/mpi/mpitools.h"
//#include "../element/fluid/elementfluid3d.h"

cOutputVTK::cOutputVTK()
{  
}

cOutputVTK::cOutputVTK(cOutputVTK &other):
cOutputfile(other)
{
}


cOutputVTK::~cOutputVTK()
{
}


void cOutputVTK::writeResultsToFile(cMesh &MyMesh, const int &NrStep, const double &Step)
{
  if (wantOutput() > 0) {
    exportSolutionToVtk(&MyMesh, NrStep, Step);
  }
}


void cOutputVTK::exportMeshToVtk(cMesh *myMesh)
{
  trace("export mesh to vtk ... ");

  std::ofstream Output;
  std::string   FilenameRoot;
  std::string   Filename;

  // --- join filename of outputfile
  cOutputfile::createFilenameRoot(FilenameRoot, myMesh);
  Filename = FilenameRoot + ".vtk";

  // -- open vtk file
  Output.open(Filename.c_str());
  exportToVtk(*myMesh, Output);
  // --- close file
  Output.close();
}

void cOutputVTK::exportSolutionToVtk(cMesh *myMesh, const PetscInt &NumStep, const PetscReal &Step)
{
  std::string   Filename;
  std::string   FilenameRoot;

  std::stringstream Temp;

  trace("export solution to vtk ...");

  if (PetscGlobalRank == 0)
  {
    // --- join filename of outputfile
    createFilenameRoot(FilenameRoot, myMesh);

    Temp << ".";
    Temp.width(8);
    Temp.fill('0');
    Temp << NumStep; 
    Filename=FilenameRoot;

    // -- open vtk file
    Temp << ".vtk";
    Filename += Temp.str();
    VtkOutput.open(Filename.c_str());
    nodeConOutput.open("nodeConnect.info");
  }

  exportToVtk(*myMesh, VtkOutput, nodeConOutput, Step);

  // --- close file
  if (PetscGlobalRank == 0)
  {
    VtkOutput.close();
    nodeConOutput.close();
  }
}

void cOutputVTK::nodesNeeded( cElementFEM *ptrElement,  int &NodesNeeded, int &VtkElementType )
{
  const eElementShape &Shape = ptrElement->getElementShape();
  switch(Shape)
  {
  case Point:
    VtkElementType = 1; //VTK_VERTEX
    NodesNeeded    = 1;
    break;
  case Beam:
    VtkElementType = 3; //VTK_LINE
    NodesNeeded    = 2;
    break;
  case Tria:
    VtkElementType = 5; //VTK_TRIAGLE
    NodesNeeded    = 3;
    break;
  case Quadrilateral:
    if(ptrElement->getNumberOfNodes()==8 || ptrElement->getNumberOfNodes()==9)
    {
      VtkElementType = 23; //VTK_QUADRATIC_QUAD
      NodesNeeded    = 8;
    }
    else
    {
      VtkElementType = 9; //VTK_QUAD
      NodesNeeded    = 4;
    }
    break;
  case Tetrahedron:
    if(ptrElement->getNumberOfNodes()==10)
    {
      VtkElementType = 24; //VTK_QUADRATIC_TETRA
      NodesNeeded    = 10;
    }
    else
    {
      VtkElementType = 10; //VTK_TETRA
      NodesNeeded    = 4;
    }
    break;
  case Hexahedron:
    if(ptrElement->getNumberOfNodes()==20 || //Brick20
      ptrElement->getNumberOfNodes()==27 )  //Brick27
    {
      VtkElementType = 25; //VTK_QUADRATIC_HEXAHEDRON
      NodesNeeded    = 20;
    }
    else
    {
      VtkElementType = 12; //VTK_HEXAHEDRON
      NodesNeeded    = 8;
    }
    break;
  case QuadSerendipity:
    trace(" warning: serendipity elements will not be");
    trace("          exported to vtk file !!");
    return;
    break;
  default:
    trace(" warning: eElementShape for vtk-output not implemented");
    break;
  }
}

void cOutputVTK::exportToVtk(cMesh &MyMesh, std::ostream &output, std::ostream &nodeConOutput, const PetscReal &Step)
{

  // ---  VTK Output of stresses at NODES only implemented for the first 3 dofs (sigma11, sigma22, sigma33) (mw);
  if(m_InitPostProcess==false) initPostProcess();

  int	VtkElementType;					// Vtk's identifier for an elementtype
  int NodesNeeded = 0;				// number of nodes that build the element
  std::stringstream nodeOutput;					// string, which contents the output of nodes
  std::stringstream elementOutput;				// string, which contents the output of elements
  std::stringstream elementtypeOutput;			// string, which contents the output of element types
  std::vector<std::stringstream*> solutionOutput;	// vector of strings, which contents the output of solutions
#ifdef PETSC_USE_COMPLEX
  std::vector<std::stringstream*> solutionOutput_imag;	// vector of strings, which contents the imaginary parts of the output of solutions 
#endif

  std::stringstream stressOutputNodes;			// string, which contents the output of stresses at the nodes
  std::stringstream shearStressOutputNodes;		// string, which contents the output of stresses at the nodes
  std::stringstream gradPOutputNodes;			// string, which contents the output of gradP at the nodes
  std::stringstream stressOutputCells;			// string, which contents the output of stresses at the cells
  std::stringstream stressOutputCellsImag;		// string, which contents the output of stresses at the cells (imaginary part)
  std::stringstream stressOutputCellsSec2;		// string, which contents the output of stresses at the cells - section 2 (for plates, lower side)
  std::stringstream stressOutputCellsSec2Imag;		// string, which contents the output of stresses at the cells (imaginary part) - section 2 (for plates, lower side)
  std::stringstream shearStressOutputCells;		// string, which contents the output of shear stresses
  std::stringstream shearStressOutputCellsImag;		// string, which contents the output of shear stresses (imaginary part)
  std::stringstream shearStressOutputCellsSec2;		// string, which contents the output of shear stresses - section 2 (for plates, lower side)
  std::stringstream shearStressOutputCellsSec2Imag;	// string, which contents the output of shear stresses (imaginary part) - section 2 (for plates, lower side)
  std::stringstream stressVonMisesOutputCells;		// string, which contents the output of von Mises stresses
  std::stringstream materialOutput;				// material id of element
  std::stringstream elemGroupOutput;				// group of each element
  std::stringstream rankOutput;					// the rank to which the element belongs

  std::set<PetscInt> NodeIds;					// ids of the nodes for elements of specific type
  std::set<PetscInt> ElementIds;				// ids of elements for a specific type

  cElementFEM *ptrElement = NULL;

  // --- loop over elements and look for the specified shape
  //     If this shape is found read as much nodes as needed to draw the element.
  //     The mid-side-nodes of quadratic elements are ignored.
  for (ItMapElements it = MyMesh.getFirstElement(); it != MyMesh.getLastElement(); it++)
  {
    ptrElement = MyMesh.getElement( it->second->getId() );
    //void material filter!!!
    if(ptrElement->getMaterial()->getIdentifier()!="void")
    {
      // first is used for export only the mesh, second one if we want to export
      // the solution at the nodes
      ElementIds.insert(it->second->getId());
      nodesNeeded(ptrElement, NodesNeeded, VtkElementType);

      if(VtkElementType==25 && NodesNeeded==20)
      {
        for (int k=0; k<12; k++) NodeIds.insert( it->second->getNode(k)->getId() );
        for (int k=16; k<20; k++) NodeIds.insert( it->second->getNode(k)->getId() );
        for (int k=12; k<16; k++) NodeIds.insert( it->second->getNode(k)->getId() );
      }
      else
      {
        for (int k=0; k<NodesNeeded; k++) NodeIds.insert( it->second->getNode(k)->getId() );
      }
    }
  }
  ptrElement = NULL;

  // --- leave if no appropriate elements are found
  //     (you need to check on all ranks)
  PetscInt GlobalSum = 0;
  PetscInt LocalSum = ElementIds.size();

  MPI_Allreduce(&LocalSum, &GlobalSum, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD);
  if (GlobalSum == 0) {
    return;
  }

  // --- collect all nodes needed to display the
  //     mesh on rank 0. First, the number of nodes found
  //     on all procs is collected on rank 0
  PetscInt LocalSize = NodeIds.size();
  std::vector<PetscInt> NodeIdsPerProc(PetscGlobalSize, 0);
  MPI_Gather( &LocalSize, 1, MPIU_INT, &(NodeIdsPerProc[0]), 1, MPIU_INT, 0, PETSC_COMM_WORLD);

  // --- now all nodes are sent to rank 0 and added to the set of nodes
  for (int rank=1; rank<PetscGlobalSize; rank++)
  {
    if (PetscGlobalRank == rank)
    {
      std::vector<PetscInt> buffer(NodeIds.size());
      std::copy(NodeIds.begin(), NodeIds.end(), buffer.begin());
      MPI_Send( &(buffer[0]), buffer.size(), MPIU_INT, 0, rank, PETSC_COMM_WORLD );
    }
    else if (PetscGlobalRank == 0)
    {
      MPI_Status status;
      std::vector<PetscInt> buffer( NodeIdsPerProc[rank], 0 );
      MPI_Recv( &(buffer[0]), buffer.size(), MPIU_INT, rank, rank, PETSC_COMM_WORLD, &status );

      for (int k=0; k<NodeIdsPerProc[rank]; k++)
      {
        NodeIds.insert( buffer[k] );
      }
    }
  }

  //Scatter solution
  Vec FullSolution;
  if(m_OutSolution!=NULL)
  {
    VecScatter ctx;
    VecScatterCreateToZero(m_OutSolution[0], &ctx, &FullSolution);   
    VecScatterBegin(ctx,m_OutSolution[0],FullSolution,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx,m_OutSolution[0],FullSolution,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterDestroy(&ctx);
  }

  Vec FullStressesNodes;
  if(m_OutStressesNodes!=NULL)
  {

    VecScatter ctx;
    VecScatterCreateToAll(m_OutStressesNodes[0], &ctx, &FullStressesNodes);
    VecScatterBegin(ctx,m_OutStressesNodes[0],FullStressesNodes,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx,m_OutStressesNodes[0],FullStressesNodes,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterDestroy(&ctx);
  }
  // --- writing header of current zone
  std::set<PetscInt>::const_iterator itSet;

  if (PetscGlobalRank == 0)
  {

    //head nodeOutput
    nodeOutput << "# vtk DataFile Version 1.0" << std::endl;
    nodeOutput << "by elPaSo" << std::endl;
    nodeOutput << "ASCII" << std::endl;
    nodeOutput << "DATASET UNSTRUCTURED_GRID" << std::endl;
    nodeOutput << "POINTS " << NodeIds.size() <<" double" << std::endl;

    //head solutionOutput
    solutionOutput.resize(cstNumberOfKnownDofs);
#ifdef PETSC_USE_COMPLEX
    solutionOutput_imag.resize(cstNumberOfKnownDofs);
#endif
    //run over all  dof
    for(int l=0; l<cstNumberOfKnownDofs; ++l)
      if(m_DispVector[l]==true)
      {
        std::stringstream *ss = new std::stringstream;
        solutionOutput[l]=ss;
        solutionOutput[l]->setf(std::ios::scientific);
#ifdef PETSC_USE_COMPLEX
        std::stringstream *ss_imag = new std::stringstream;
        solutionOutput_imag[l]=ss_imag;
        solutionOutput_imag[l]->setf(std::ios::scientific);
#endif
      }
      //disp_x
      if(m_dispX>0)
      {
        if(!solutionOutput[0])
        {
          std::stringstream *ss = new std::stringstream;
          solutionOutput[0]=ss;
          solutionOutput[0]->setf(std::ios::scientific);
        }
#ifdef PETSC_USE_COMPLEX
        if(!solutionOutput_imag[0])
        {
          std::stringstream *ss_imag = new std::stringstream;
          solutionOutput_imag[0]=ss_imag;
          solutionOutput_imag[0]->setf(std::ios::scientific);
        }
#endif
        for(int i=0;i<3;++i)m_DispVector[i] = true;
        (*solutionOutput[0]) << "VECTORS disp_x float"<< std::endl;
#ifdef PETSC_USE_COMPLEX
        for(int i=0;i<3;++i)m_DispVector[i] = true;
        (*solutionOutput_imag[0]) << "VECTORS disp_x_imag float"<< std::endl;
#endif

      }
      //disp_w
      if(m_dispW>0)
      {
        if(!solutionOutput[3])
        {
          std::stringstream *ss = new std::stringstream;
          solutionOutput[3]=ss;
          solutionOutput[3]->setf(std::ios::scientific);
        }
#ifdef PETSC_USE_COMPLEX
        if(!solutionOutput_imag[3])
        {
          std::stringstream *ss_imag = new std::stringstream;
          solutionOutput_imag[3]=ss_imag;
          solutionOutput_imag[3]->setf(std::ios::scientific);
        }
#endif

        for(int i=3;i<6;++i)m_DispVector[i] = true;
        (*solutionOutput[3])<< "VECTORS disp_w float"<< std::endl;
#ifdef PETSC_USE_COMPLEX
        for(int i=3;i<6;++i)m_DispVector[i] = true;
        (*solutionOutput_imag[3]) << "VECTORS disp_w_imag float"<< std::endl;
#endif

      }
      //all others
      for(int l=6; l<cstNumberOfKnownDofs; ++l){
        if(m_DispVector[l]==true)
        {
          (*solutionOutput[l]) <<"SCALARS "<<sKnownDofsIdentifiers[l]<<" float 1"<<std::endl;
          (*solutionOutput[l]) <<"LOOKUP_TABLE default"<<std::endl;
#ifdef PETSC_USE_COMPLEX
          (*solutionOutput_imag[l]) <<"SCALARS "<<sKnownDofsIdentifiers[l]<<"_imag float 1"<<std::endl;
          (*solutionOutput_imag[l]) <<"LOOKUP_TABLE default"<<std::endl;
#endif

        }
      }
      nodeOutput.setf(std::ios::scientific);
      // --- output of nodal coordinates
      //     we are only printing the nodes we've found
      //     in the sought for elements

      cNode *ptrNode = NULL;
      for (itSet=NodeIds.begin(); itSet!=NodeIds.end(); itSet++)
      {
        ptrNode = MyMesh.getNode( *itSet );
        nodeConOutput << ptrNode->getId() << std::endl;

        // --- collect output of node's coordinates and solutions
        for(int l=0; l<cstNumberOfKnownDofs; ++l)
        {
          PetscInt    row;
          PetscScalar value = 0.;
          if(m_DispVector[l] == true)
          {
            if(ptrNode->checkIfActive( l ) == true)
            {
              row = ptrNode->getGlobalRow( l );
              if(m_OutSolution) VecGetValues(FullSolution, 1, &row, &value);
#ifdef PETSC_USE_COMPLEX
              if(l<3)
              {
                (*solutionOutput[0]) << value.real() << " ";
                (*solutionOutput_imag[0]) << value.imag() << " ";
              }
              else if(l<6)
              {
                (*solutionOutput[3]) << value.real() << " ";
                (*solutionOutput_imag[3]) << value.imag() << " ";
              }
              else
              {
                (*solutionOutput[l]) << value.real() << " ";
                (*solutionOutput_imag[l]) << value.imag() << " ";
              }
#else
              if(l<3)      (*solutionOutput[0]) << value << " ";
              else if(l<6) (*solutionOutput[3]) << value << " ";
              else         (*solutionOutput[l]) << value << " ";
#endif
            }
            else
            {
              if(l<3)      (*solutionOutput[0]) << "0.0 ";
              else if(l<6) (*solutionOutput[3]) << "0.0 ";
              else		     (*solutionOutput[l]) << "0.0 ";
#ifdef PETSC_USE_COMPLEX
              if(l<3)	     (*solutionOutput_imag[0]) << "0.0 ";
              else if(l<6) (*solutionOutput_imag[3]) << "0.0 ";
              else         (*solutionOutput_imag[l]) << "0.0 ";
#endif
            }
            if(l<3)
            {
              PetscReal xi = (*ptrNode)[l];
              // here is the output for the nodes
              // do not add the displacement values!!!
#ifdef PETSC_USE_COMPLEX
              nodeOutput << std::showpoint << xi << " ";
#else
              nodeOutput << std::showpoint << xi  << " ";
#endif
            }
          }
          else if(l<3) //testing
          {
#ifdef PETSC_USE_COMPLEX
            PetscScalar xi = (*ptrNode)[l];
            nodeOutput << std::showpoint << xi.real() << " ";
#endif
          }
        }
        nodeOutput << std::endl;
        if(m_DispVector[0] == true) (*solutionOutput[0]) << std::endl;
#ifdef PETSC_USE_COMPLEX
        if(m_DispVector[0] == true) (*solutionOutput_imag[0]) << std::endl;
#endif
        if(m_DispVector[3] == true) (*solutionOutput[3]) << std::endl;
#ifdef PETSC_USE_COMPLEX
        if(m_DispVector[3] == true) (*solutionOutput_imag[3]) << std::endl;
#endif
        for(int l=6; l<cstNumberOfKnownDofs; ++l) {
          if(m_DispVector[l] == true) (*solutionOutput[l]) << std::endl;
#ifdef PETSC_USE_COMPLEX
          if(m_DispVector[l] == true) (*solutionOutput_imag[l]) << std::endl;
#endif
        }
        //--- stress Nodes
        if(m_OutStressesNodes!=NULL)
        {
          stressOutputNodes.setf(std::ios::scientific);
          shearStressOutputNodes.setf(std::ios::scientific);
          gradPOutputNodes.setf(std::ios::scientific);
          PetscInt pos;
          PetscScalar value;
          for(int i=0; i<3; ++i) 
          {
            pos=(ptrNode->getGlobalSeqId())*9+i;
            VecGetValues(FullStressesNodes, 1, &pos, &value);

#ifdef PETSC_USE_COMPLEX
            stressOutputNodes << value.real() << " ";
#else
            stressOutputNodes << value << " ";
#endif
          }
          for(int i=3; i<6; ++i) 
          {
            pos=(ptrNode->getGlobalSeqId())*9+i;
            VecGetValues(FullStressesNodes, 1, &pos, &value);
#ifdef PETSC_USE_COMPLEX
            shearStressOutputNodes << value.real() << " ";
#else
            shearStressOutputNodes << value << " ";
#endif
          }
          for(int i=6; i<9; ++i) 
          {
            pos=(ptrNode->getGlobalSeqId())*9+i;
            VecGetValues(FullStressesNodes, 1, &pos, &value);
#ifdef PETSC_USE_COMPLEX
            gradPOutputNodes << value.real() << " ";
#else
            gradPOutputNodes << value << " ";
#endif
          }
          stressOutputNodes << std::endl;
          shearStressOutputNodes << std::endl;
          gradPOutputNodes << std::endl;
        }
      }
      ptrNode = NULL;
  }//end rank 0 only

  //--- FullSolution not longer needed
  if(m_OutSolution!=NULL)
  {
    VecDestroy(&FullSolution);
  }
  if(m_OutStressesNodes!=NULL)
  {
    VecDestroy(&FullStressesNodes);
  }
  // --- copy node ids into a vector, first only on
  //     rank 0. In order to get correct connectivity
  //     tables the vector must be sent to all other
  //     ranks.
  PetscInt VecSizeRank0 = NodeIds.size(); // correct value currently only on rank 0 !!
  MPI_Bcast ( &VecSizeRank0, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  std::vector<PetscInt> VecNodeIds( VecSizeRank0 );

  if (PetscGlobalRank == 0)
    std::copy(NodeIds.begin(), NodeIds.end(), VecNodeIds.begin());

  MPI_Bcast ( &(VecNodeIds[0]), VecSizeRank0, MPIU_INT, 0, PETSC_COMM_WORLD);

  //--- output of the elements
  ptrElement = NULL;
  //  PetscInt  puffer[8];

  // --- output of the elements and stresses
  //Scatter stresses
  Vec FullStressesCells;
  Vec FullStressesCellsSec2;
  if(m_OutStressesCells!=NULL)
  {
    VecScatter ctx;
    VecScatterCreateToAll(m_OutStressesCells[0], &ctx, &FullStressesCells);
    VecScatterBegin(ctx,m_OutStressesCells[0],FullStressesCells,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx,m_OutStressesCells[0],FullStressesCells,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterDestroy(&ctx);
    // std::cout << stressOutputCells << std::endl;
    stressOutputCells.setf(std::ios::scientific);
    stressOutputCellsImag.setf(std::ios::scientific);
    shearStressOutputCells.setf(std::ios::scientific);
    shearStressOutputCellsImag.setf(std::ios::scientific);
    // std::cout << stressOutputCells << std::endl;
  }
  if(m_OutStressesCellsSec2!=NULL)
  {
    VecScatter ctx;
    VecScatterCreateToAll(m_OutStressesCellsSec2[0], &ctx, &FullStressesCellsSec2);
    VecScatterBegin(ctx,m_OutStressesCellsSec2[0],FullStressesCellsSec2,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx,m_OutStressesCellsSec2[0],FullStressesCellsSec2,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterDestroy(&ctx);
    // std::cout << stressOutputCells << std::endl;
    stressOutputCellsSec2.setf(std::ios::scientific);
    stressOutputCellsSec2Imag.setf(std::ios::scientific);
    shearStressOutputCellsSec2.setf(std::ios::scientific);
    shearStressOutputCellsSec2Imag.setf(std::ios::scientific);
    // std::cout << stressOutputCells << std::endl;
  }

  if (PetscGlobalRank == 0)
  {
    elementtypeOutput << "CELL_TYPES " << GlobalSum << std::endl;
  }
  int i=0; //counts number of values
  for (itSet=ElementIds.begin(); itSet!=ElementIds.end(); itSet++)
  {
    ptrElement = MyMesh.getElement( *itSet );

    nodesNeeded(ptrElement, NodesNeeded, VtkElementType);

    elementOutput << NodesNeeded << " ";
    elementtypeOutput << VtkElementType << " ";

    if(VtkElementType==25 && NodesNeeded==20)
    {
      for (int k=0; k<12; k++) elementOutput << (int)(std::find(VecNodeIds.begin(), VecNodeIds.end(), ptrElement->getNode(k)->getId()) - VecNodeIds.begin()) << " ";
      for (int k=16; k<20; k++) elementOutput << (int)(std::find(VecNodeIds.begin(), VecNodeIds.end(), ptrElement->getNode(k)->getId()) - VecNodeIds.begin()) << " ";
      for (int k=12; k<16; k++) elementOutput << (int)(std::find(VecNodeIds.begin(), VecNodeIds.end(), ptrElement->getNode(k)->getId()) - VecNodeIds.begin()) << " ";
    }
    else
    {
      for (int k=0; k<NodesNeeded; k++)
      {
        // --- now compute the index of the value within the vector
        //     thx to de.comp.lang.iso-c++ :-)
        elementOutput << (int)(std::find(VecNodeIds.begin(), VecNodeIds.end(), ptrElement->getNode(k)->getId()) - VecNodeIds.begin()) << " ";
      }
    }
    elementOutput << " " << std::endl;

    //--- stress
    if(m_OutStressesCells!=NULL)
    {
      PetscInt pos[6];
      PetscScalar value[6];
      //PetscScalar pre_value[6];
      for(int i=0; i<6; ++i)
      {
        pos[i]=(ptrElement->getId0n())*9+i;
      }
      VecGetValues(FullStressesCells, 6, &pos[0], &value[0]);
      // Export real and imag part part if elpasoC is used; Else the given value as real part. (C.Blech)
#ifdef PETSC_USE_COMPLEX
      for(int i=0; i<6; ++i)
      {
        if(i<3)
        {
          stressOutputCells << PetscRealPart(value[i]) << " ";
          stressOutputCellsImag << PetscImaginaryPart(value[i]) << " ";
        }
        else
        {
          shearStressOutputCells << PetscRealPart(value[i]) << " ";
          shearStressOutputCellsImag << PetscImaginaryPart(value[i]) << " ";
        }
      }
      stressVonMisesOutputCells <<
        std::sqrt(PetscRealPart(value[0])*PetscRealPart(value[0])+PetscRealPart(value[1])*PetscRealPart(value[1])+PetscRealPart(value[2])*PetscRealPart(value[2])
      -PetscRealPart(value[0])*PetscRealPart(value[1])-PetscRealPart(value[0])*PetscRealPart(value[2])-PetscRealPart(value[1])*PetscRealPart(value[2])
      +3.*(PetscRealPart(value[3])*PetscRealPart(value[3])+PetscRealPart(value[4])*PetscRealPart(value[4])+PetscRealPart(value[5])*PetscRealPart(value[5])))<< " ";
      stressOutputCells << std::endl;
      stressOutputCellsImag << std::endl;
      shearStressOutputCells << std::endl;
      shearStressOutputCellsImag << std::endl;
      stressVonMisesOutputCells << std::endl;
#else
      for(int i=0; i<6; ++i)
      {
        if(i<3)
        {
          stressOutputCells << value[i] << " ";
        }
        else
        {
          shearStressOutputCells << value[i] << " ";
        }
      }
      stressVonMisesOutputCells <<
        std::sqrt(value[0]*value[0]+value[1]*value[1]+value[2]*value[2]
      -value[0]*value[1]-value[0]*value[2]-value[1]*value[2]
      +3.*(value[3]*value[3]+value[4]*value[4]+value[5]*value[5]))<< " ";
      stressOutputCells << std::endl;
      shearStressOutputCells << std::endl;
      stressVonMisesOutputCells << std::endl;
#endif

    }
    //--- stress for section 2
    if(m_OutStressesCellsSec2!=NULL)
    {
      PetscInt pos[6];
      PetscScalar value[6];
      //PetscScalar pre_value[6];
      for(int i=0; i<6; ++i)
      {
        pos[i]=(ptrElement->getId0n())*9+i;
      }
      VecGetValues(FullStressesCellsSec2, 6, &pos[0], &value[0]);
      // Export real and imag part part if elpasoC is used; Else the given value as real part. (C.Blech)
#ifdef PETSC_USE_COMPLEX
      for(int i=0; i<6; ++i)
      {
        if(i<3)
        {
          stressOutputCellsSec2 << PetscRealPart(value[i]) << " ";
          stressOutputCellsSec2Imag << PetscImaginaryPart(value[i]) << " ";
        }
        else
        {
          shearStressOutputCellsSec2 << PetscRealPart(value[i]) << " ";
          shearStressOutputCellsSec2Imag << PetscImaginaryPart(value[i]) << " ";
        }
      }
      stressOutputCellsSec2 << std::endl;
      stressOutputCellsSec2Imag << std::endl;
      shearStressOutputCellsSec2 << std::endl;
      shearStressOutputCellsSec2Imag << std::endl;
#else
      for(int i=0; i<6; ++i)
      {
        if(i<3)
        {
          stressOutputCellsSec2 << value[i] << " ";
        }
        else
        {
          shearStressOutputCellsSec2 << value[i] << " ";
        }
      }
      stressOutputCellsSec2 << std::endl;
      shearStressOutputCellsSec2 << std::endl;
#endif

    }

    //--- material ID and element's group ID in ak3
    if(ptrElement->getHistory()!=NULL)
    {
      if(ptrElement->getHistory()->getYielding()) materialOutput<<-1<< " ";
      else materialOutput<<ptrElement->getMaterial()->getId() << " ";
    }
    else materialOutput<<ptrElement->getMaterial()->getId() << " ";
    ItMapGroups itGroups;
    std::set<PetscInt>::const_iterator itElSet;
    for (itGroups = MyMesh.getFirstGroup(); itGroups != MyMesh.getLastGroup(); itGroups++)
    {
      for (itElSet = itGroups->second->getFirstElementId(); itElSet != itGroups->second->getLastElementId(); itElSet++)
      {
        if(MyMesh.getElement( *itElSet )->getId()==ptrElement->getId())
        {
          elemGroupOutput<< itGroups->second->getId() << " ";
        }
      }
    }
    //--- rank
    rankOutput<< PetscGlobalRank << " ";
    i=i+1+NodesNeeded;
  }
  elementtypeOutput << std::endl;
  ptrElement = NULL;

  //--- FullStressesCells not longer needed
  if(m_OutStressesCells!=NULL)
  {
    VecDestroy(&FullStressesCells);
  }
  if(m_OutStressesCellsSec2!=NULL)
  {
    VecDestroy(&FullStressesCellsSec2);
  }

  // add i from all ranks
  PetscInt NumOfValues = 0;
  MPI_Allreduce(&i, &NumOfValues, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD);

  //--- wirte nodes into output file
  if (PetscGlobalRank == 0)
  {
    output << nodeOutput.str();
  }

  //--- wirte cells into output file
  std::stringstream cells;
  cells << "CELLS "<< GlobalSum << " " << NumOfValues;
  collectOutputStrings(output, elementOutput, cells.str(), "");

  //--- write cell type into output file
  collectOutputStrings(output, elementtypeOutput, "", "");

  //--- Solution Point Data
  if( PetscGlobalRank == 0)
  {
    if(m_OutSolution!=NULL)
    {
      output << "POINT_DATA " << NodeIds.size() << std::endl;
      for(int l=0; l<cstNumberOfKnownDofs; ++l) {
        if(solutionOutput[l]) output << (*solutionOutput[l]).str();
#ifdef PETSC_USE_COMPLEX
        if(solutionOutput_imag[l]) output << (*solutionOutput_imag[l]).str();
#endif

      }
    }

    if(m_OutStressesNodes!=NULL)
    {
      output << "VECTORS stressNodes float" <<std::endl;
      output << stressOutputNodes.str();

      output << "VECTORS shearStressNodes float" <<std::endl;
      output << shearStressOutputNodes.str();

      output << "VECTORS gradPNodes float" <<std::endl;
      output << gradPOutputNodes.str();
    }

    //--- delete strings
    for(int l=0; l<cstNumberOfKnownDofs; ++l) {
      if((solutionOutput[l])!=NULL) delete (solutionOutput[l]);
#ifdef PETSC_USE_COMPLEX
      if((solutionOutput_imag[l])!=NULL) delete (solutionOutput_imag[l]);
#endif

    }
  }

  //--- CELL_DATA...
  output << "CELL_DATA "<< GlobalSum <<std::endl;
  if(m_OutStressesCells!=NULL)
  {
    //--- stresses at cells
    collectOutputStrings(output, stressOutputCells, "VECTORS stressCells float", "");
    collectOutputStrings(output, stressOutputCellsImag, "VECTORS stressCells_imag float", "");
    //--- shear stresses at cells
    collectOutputStrings(output, shearStressOutputCells, "VECTORS shearstress float", "");
    collectOutputStrings(output, shearStressOutputCellsImag, "VECTORS shearstress_imag float", "");
    collectOutputStrings(output, stressVonMisesOutputCells, "SCALARS vonMises float", "LOOKUP_TABLE default");
  }
  if(m_OutStressesCellsSec2!=NULL)
  {
    //--- stresses at cells (for section 2)
    collectOutputStrings(output, stressOutputCellsSec2, "VECTORS stressCellsSec2 float", "");
    collectOutputStrings(output, stressOutputCellsSec2Imag, "VECTORS stressCellsSec2_imag float", "");
    //--- shear stresses at cells (for section 2)
    collectOutputStrings(output, shearStressOutputCellsSec2, "VECTORS shearstressSec2 float", "");
    collectOutputStrings(output, shearStressOutputCellsSec2Imag, "VECTORS shearstressSec2_imag float", "");
  }
  //output << std::endl;
  //collect material information
  collectOutputStrings(output, materialOutput, "SCALARS MatID float", "LOOKUP_TABLE default");
  output << std::endl;
  //collect elem group ID information
  collectOutputStrings(output, elemGroupOutput, "SCALARS elemGroupID float", "LOOKUP_TABLE default");
  output << std::endl;
  //collect rank information
  collectOutputStrings(output, rankOutput,     "SCALARS rank float",  "LOOKUP_TABLE default");
  output << std::endl;
}

void cOutputVTK::collectOutputStrings(std::ostream &output, std::stringstream &string, std::string tag1, std::string tag2)
{
  //--- rank 0
  if( PetscGlobalRank == 0)
  {
    if(tag1.size()!=0)output << tag1 << std::endl;
    if(tag2.size()!=0)output << tag2 << std::endl;
    output << string.str();
  }

  //--- collect information
  for (int proc=1; proc<PetscGlobalSize; proc++)
  {
    if (PetscGlobalRank == proc)
    {
      cMpiTools::sendStringValue(string.str(),0,9600+proc,PETSC_COMM_WORLD);
    }
    else if(PetscGlobalRank == 0)
    {
      std::string recvbuf;
      recvbuf = cMpiTools::receiveStringValue(proc,9600+proc, PETSC_COMM_WORLD);
      output<<recvbuf;
    }
  }
}

void cOutputVTK::writeMeshToFile(cMesh &MyMesh)
{
  if (wantOutput() > 0) {

    exportMeshToVtk(&MyMesh);
  }
}

void cOutputVTK::writeMeshToFileGroup(cMesh &MyMesh)
{
  if (wantOutput() > 0) {

    trace("no mesh output for groups possible in VTK-format ...");
  }
}

//! Creates a copy of current Object
cOutputfile* cOutputVTK::copy() {

  cOutputfile* copy_object = new cOutputVTK(*this);

  return copy_object;

}

std::ostream& cOutputVTK::write(std::ostream &os) const
{
  os << " want a vtk file                 : " << wantOutput() << std::endl;

  return os;
}


std::istream& cOutputVTK::read(std::istream &is)
{
  return is;
}


std::ostream& cOutputVTK::writeXml(std::ostream &os) const
{
  os << "  <vtk>" << wantOutput() << "</vtk>" << std::endl;

  return os;
}
