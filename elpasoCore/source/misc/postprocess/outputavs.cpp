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

#include "outputavs.h"

cOutputAVS::cOutputAVS()
{  
}

cOutputAVS::cOutputAVS(cOutputAVS &other):
cOutputfile(other)
{
}

cOutputAVS::~cOutputAVS()
{
}


void cOutputAVS::writeResultsToFile(cMesh &MyMesh, const int &NrStep, const double &Step)
{
  if (wantOutput() > 0) {

    exportSolutionToAVS(&MyMesh, NrStep, Step);
  }
}

void cOutputAVS::writeMeshToFile(cMesh &MyMesh)
{
  if (wantOutput() > 0) {

    exportMeshToAVS(&MyMesh);
  }
}

void cOutputAVS::writeMeshToFileGroup(cMesh &MyMesh)
{
  if (wantOutput() > 0) {

    trace("no mesh output for groups possible in AVS-format ...");
  }
}


void cOutputAVS::exportMeshToAVS(cMesh *myMesh)
{
  trace("export solution to AVS ...");

  std::ofstream Output;
  std::string   FilenameRoot;
  std::string   Filename;

  std::stringstream Temp;

  // --- join filename of outputfile
  createFilenameRoot(FilenameRoot, myMesh);
  Temp << "_r";
  Temp << PetscGlobalRank;
  Temp << ".inp";
  Filename = FilenameRoot + Temp.str();


  // -- open AVS INP file
  Output.open(Filename.c_str());
  exportToAVS(*myMesh, Output);
  // --- close file
  Output.close();
}


void cOutputAVS::exportSolutionToAVS(cMesh *myMesh, const PetscInt &NumStep, const PetscReal &Step)
{
  std::string   Filename;
  std::string   FilenameRoot;

  std::stringstream fileExt;

  trace("export solution to avs ...");

  // --- join filename of outputfile
  createFilenameRoot(FilenameRoot, myMesh);
  fileExt << "_r";
  fileExt << PetscGlobalRank;
  fileExt << ".";
  fileExt.width(8);
  fileExt.fill('0');
  fileExt << NumStep;
  fileExt << ".inp";

  Filename=FilenameRoot;
  Filename += fileExt.str();

  // -- open avs file
  AVSOutput.open(Filename.c_str());//, std::ios_base::out|std::ios_base::binary);
  exportToAVS(*myMesh, AVSOutput, Step);
  // --- close file
  AVSOutput.close();
}

void cOutputAVS::exportToAVS(cMesh &MyMesh, std::ostream &output, const PetscReal &Step /* = 0.0 */, ePhysicsType PhysicsType /* = Undefined */)
{
  if(m_InitPostProcess==false) initPostProcess();

  std::string	  AVSElementType;			// AVS's identifier for an elementtype
  int			      NodesNeeded = 0;    // number of nodes that build the element
  std::stringstream headOutput;		  // string, which contains file information
  std::stringstream nodeOutput;		  // string, which contains the output of nodes
  std::stringstream elementOutput;	// string, which contains the output of elements
  std::stringstream solutionOutput;	// string, which contains the output of solutions

  std::set<PetscInt> NodeIds;			// ids of the nodes for elements of specific type
  std::set<PetscInt> ElementIds;		// ids of elements for a specific type

  cElementFEM *ptrElement = NULL;

  // --- loop over elements and look for the specified shape
  //     If this shape is found read as much nodes as needed to draw the element.
  //     The mid-side-nodes of quadratic elements are ignored.
  for (ItMapElements it = MyMesh.getFirstElement(); it != MyMesh.getLastElement(); it++) {
    ptrElement = MyMesh.getElement( it->second->getId() );
    //void material filter!!!
    if(ptrElement->getMaterial()->getIdentifier()!="void")
    {
      // first is used for export only the mesh, second one if we want to export
      // the solution at the nodes
      if ((PhysicsType == 0) || (it->second->getPhysicsType() == PhysicsType)) {
        ElementIds.insert(it->second->getId());
        NodesNeeded = it->second->getNumberOfNodes();
        for (int k=0; k<NodesNeeded; k++) {
          NodeIds.insert( it->second->getNode(k)->getId() );
        }
      }	
    }
  }
  ptrElement = NULL;

  Vec FullSolution;

  if(m_OutSolution!=NULL)
  {
    VecScatter ctx;
    VecScatterCreateToAll(m_OutSolution[0], &ctx, &FullSolution);
    VecScatterBegin(ctx,m_OutSolution[0],FullSolution,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx,m_OutSolution[0],FullSolution,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterDestroy(&ctx);
  }

  // --- writing header of current zone
  std::set<PetscInt>::const_iterator itSet;

  int i = 1;

  //head nodeOutput
  headOutput << "# UCD-File created by ELPASO" << std::endl;
  headOutput << NodeIds.size() <<" ";
  headOutput << ElementIds.size() <<" ";
  headOutput << m_numberofDofs << " ";  //Number of known degrees of freedom
  headOutput << "0 0" << std::endl;

  //head solutionOutput
  solutionOutput << m_numberofDofs ;
  if(m_dispX>0) solutionOutput<< " " << m_dispX;
  if(m_dispW>0) solutionOutput<< " " << m_dispW;
  for(int l=6; l<cstNumberOfKnownDofs; ++l) if(m_DispVector[l]==true)solutionOutput <<" 1";
  solutionOutput  << std::endl;
  if(m_dispX>0) solutionOutput <<"disp_x, "<< std::endl;
  if(m_dispW>0) solutionOutput <<"disp_w, "<< std::endl;
  for(int l=6; l<cstNumberOfKnownDofs; ++l) if(m_DispVector[l]==true)solutionOutput<<sKnownDofsIdentifiers[l]<<", "<<std::endl;

  nodeOutput.setf(std::ios::scientific);
  solutionOutput.setf(std::ios::scientific);

  // --- output of nodal coordinates
  //     we are only printing the nodes we've found
  //     in the sought for elements
  cNode *ptrNode = NULL;

  for (itSet=NodeIds.begin(); itSet!=NodeIds.end(); itSet++)
  {
    ptrNode = MyMesh.getNode( *itSet );

    nodeOutput << i << " ";
    solutionOutput << i << " ";
    i++;
    // --- collect output of node's coordinates and solutions
    for(int l=0; l<cstNumberOfKnownDofs; ++l)
    {
      PetscInt    row;
      PetscScalar value = 0.;
      if(m_DispVector[l]==true)
      {
        if(ptrNode->checkIfActive( l ) == true)
        {
          row = ptrNode->getGlobalRow( l );
          if(m_OutSolution) VecGetValues(FullSolution, 1, &row, &value);
#ifdef PETSC_USE_COMPLEX
          solutionOutput << value.real() << " ";
#else
          solutionOutput << value << " ";
#endif
        }
        else
        {
          solutionOutput << "0.0 ";
        }
        if(l<3)
        {
          PetscReal xi = (*ptrNode)[l];
#ifdef PETSC_USE_COMPLEX
          nodeOutput << std::showpoint << xi + value.real() << " ";	
#else
          nodeOutput << std::showpoint << xi + value << " ";	
#endif
        }
      }
      else if(l<3) //testing
      {
        PetscScalar xi = (*ptrNode)[l];
#ifdef PETSC_USE_COMPLEX
        nodeOutput << std::showpoint << xi.real() << " ";
#else
        nodeOutput << std::showpoint << xi << " ";
#endif
      }
    }
    nodeOutput << std::endl;
    solutionOutput << std::endl;
  }
  ptrNode = NULL;


  if(m_OutSolution!=NULL)
  {
    VecDestroy(&FullSolution);
  }

  // --- copy node ids into a vector, first only on
  //     rank 0. In order to get correct connectivity
  //     tables the vector must be sent to all other
  //     ranks.
  PetscInt VecSizeRank0 = NodeIds.size(); // correct value currently only on rank 0 !!
  MPI_Bcast ( &VecSizeRank0, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  std::vector<PetscInt> VecNodeIds( VecSizeRank0 );

  std::copy(NodeIds.begin(), NodeIds.end(), VecNodeIds.begin());

  // --- output of the elements	
  i = 1;
  for (itSet=ElementIds.begin(); itSet!=ElementIds.end(); itSet++)
  {
    ptrElement = MyMesh.getElement( *itSet );
    const eElementShape &Shape = ptrElement->getElementShape();
    switch(Shape)
    {
    case Beam:
      AVSElementType = "line";
      NodesNeeded        = 2;
      break;
    case Tria:
      AVSElementType = "tri";
      NodesNeeded        = 3;
      break;
    case Quadrilateral:
      AVSElementType = "quad";
      NodesNeeded        = 4;
      break;
    case Tetrahedron:
      AVSElementType = "tet";
      NodesNeeded        = 4;
      break;
    case Hexahedron:
      AVSElementType = "hex";
      NodesNeeded    = 8;
      break;
    case QuadSerendipity:
      trace(" warning: serendipity elements will not be");
      trace("          exported to AVS inp file !!");
      return;
      break;
    }
    elementOutput << i << " " << ptrElement->getMaterial()->getId() << " " << AVSElementType << " ";

    for (int k=0; k<NodesNeeded; k++)
    {
      // --- now compute the index of the value within the vector
      //     thx to de.comp.lang.iso-c++ :-)
      elementOutput << (int)(std::find(VecNodeIds.begin(), VecNodeIds.end(), ptrElement->getNode(k)->getId()) - VecNodeIds.begin()) + 1 << " ";
    }
    elementOutput << " " << std::endl;
    i++;
  }
  ptrElement = NULL;

  output << headOutput.str();
  output << nodeOutput.str();
  output << elementOutput.str();
  output << solutionOutput.str();
}


//! Creates a copy of current Object
cOutputfile* cOutputAVS::copy() {

  cOutputfile* copy_object = new cOutputAVS(*this);

  return copy_object;

}


std::ostream& cOutputAVS::write(std::ostream &os) const
{
  os << " want a AVS file                 : " << wantOutput() << std::endl;

  return os;
}


std::istream& cOutputAVS::read(std::istream &is)
{
  return is;
}


std::ostream& cOutputAVS::writeXml(std::ostream &os) const
{
  os << "  <avs>" << wantOutput() << "</avs>" << std::endl;

  return os;
}
