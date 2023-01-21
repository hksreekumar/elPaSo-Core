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

#include "elementcoupling.h"

cElementCoupling::cElementCoupling(const short &NumberOfNodes) :
cElement(NumberOfNodes),
m_IdsMatchingNodes(NumberOfNodes, 0),
m_MatchingBemNodes(NumberOfNodes)
{
	setOrientation( 0. );

	for (int k=0; k < NumberOfNodes; k++)
		m_MatchingBemNodes[k] = NULL;
}


cElementCoupling::cElementCoupling(const cElementCoupling &other) :
cElement( other ),
m_IdsMatchingNodes( other.getNumberOfNodes() ),
m_MatchingBemNodes( other.getNumberOfNodes() )
{
	setOrientation( other.getOrientation() );
	for (int k=0; k < other.getNumberOfNodes(); k++)
	{
		setMatchingBemNode(k, other.getMatchingBemNode(k));
		setIdMatchingBemNode( k, other.getIdMatchingBemNode(k) );
	}
}


cElementCoupling::~cElementCoupling()
{
	for (int k=0; k < getNumberOfNodes(); k++)
		m_MatchingBemNodes[k] = NULL;
}


swebem::cColPoint* cElementCoupling::getMatchingBemNode(int Index) const
{
#ifdef PETSC_USE_DEBUG
	return m_MatchingBemNodes.at(Index);
#else
	return m_MatchingBemNodes[Index];
#endif
}


void cElementCoupling::setMatchingBemNode(int Index, swebem::cColPoint *ptrNode)
{
#ifdef PETSC_USE_DEBUG
	m_MatchingBemNodes.at(Index) = ptrNode;
#else
	m_MatchingBemNodes[Index] = ptrNode;
#endif
}

/// @todo implement for triangles
/// @todo implement for arbitrary configurations (consider all displacement directions)
void cElementCoupling::computeCouplingMatrix(Mat &K, Vec &F)
{
#ifdef HAVE_SWEBEM
	const PetscInt  nnod = getNumberOfNodes();  // number of nodes of this element
	const int       ngp  = (nnod == 4 ? 2 : 3); // number of Gauss points used
	cPoint          gp;                         // current Gauss point
	PetscReal       N[9], Nxi[9], Neta[9];      // evaluated testfunctions
	cVector         a1(3);                      // basis of tangential plane at ...
	cVector         a2(3);                      // ... current Gauss point
	PetscReal       detJac;                     // Jacobian
	PetscReal       weight;                     // weights of numerical integration
	PetscReal       wdJOr;                      // = weight * detJac * orientation of normals
	cElementMatrix  C(3*nnod, nnod);		// Coupling matrix for all 3 dimensions

	std::vector<PetscInt> rows(3*nnod, -1); // corresponding to FEM domain
	std::vector<PetscInt> cols(nnod, -1); // corresponding to BEM domain

	for (int k=0; k<nnod; k++)
	{
		rows[k]		= m_Nodes[k]->getGlobalRow(disp_x1);
		rows[k+nnod]	= m_Nodes[k]->getGlobalRow(disp_x2);
		rows[k+2*nnod]	= m_Nodes[k]->getGlobalRow(disp_x3);
		cols[k] = m_MatchingBemNodes[k]->getGlobalRow();
	}

	// -------------------------------------------------------------------------
	//   perform numerical integration
	//   numerical integration is performed for x1 direction and 
	//   then copied to other directions x2 and x3
	// -------------------------------------------------------------------------
	for (int n=0; n < ngp*ngp; n++)
	{
		// --- evaluate testfunctions
		gp     = m_GaussPoints.getGaussPoint2D(ngp, n);
		weight = m_GaussPoints.getGaussWeight2D(ngp, n);

		for (int k=0; k<nnod; k++)
		{
			N[k]    = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nnod, k, N_fun, gp);
			Nxi[k]  = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nnod, k, N_xi , gp);
			Neta[k] = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nnod, k, N_eta, gp);
		}

		// --- compute Jacobian
		a1.setValue(0.0);
		a2.setValue(0.0);

		for (int k=0; k<nnod; k++)
		{
			for (int d=0; d<3; d++)
			{
				a1[d] += Nxi[k]  * (*m_Nodes[k])[d];
				a2[d] += Neta[k] * (*m_Nodes[k])[d];
			}
		}

		detJac = sqrt( a1.abs2() * a2.abs2() - (a1.dot(a2))*(a1.dot(a2)));
		wdJOr  = weight * detJac * getOrientation();

		for (int z=0; z<nnod; z++)
			for (int s=0; s<nnod; s++)
				C(z,s) += N[z] * N[s] * wdJOr;
	} // end numerical integration

	// --- copy coupling terms to x2 and x3 part of the coupling matrix
	for (int z=0; z<nnod; z++)
		for (int s=0; s<nnod; s++)
		{
			C(z+nnod  ,s) = C(z,s);
			C(z+2*nnod,s) = C(z,s);
		}

		/// @todo apply inhomogenous nodal boundary conditions here (may need to modify loadvector)
		// --- check boundary conditions
		cBoundaryConditionStructure *ptrBC = NULL;
		for (int k=0; k<nnod; k++)
		{
			// add factor for all three dimensions 
			// Factor was "1" when only x3 direction (and accordingly pressure in that direction)
			// was considered.
			cCoord normalAtNode=m_MatchingBemNodes[k]->getNormal();

			/*   // Output for testing purposes
			std::cout<<"GetNormal:  " << normalAtNode.x[0] << "   " << normalAtNode.x[1] << "   "<< normalAtNode.x[2] << " ... ";
			*/ 

			// since the "normalisation" is done with help of the BEM normal which does not 
			// necessarily point in the same direction as the FEM normal, the length is 
			// multiplied by the orientation factor (derived from coupling elements);
			PetscReal len=normalAtNode.len() * getOrientation();

			for (int s=0; s<nnod; s++)
			{
				C(k       ,s) *= (normalAtNode.x[0]/len);
				C(k+nnod  ,s) *= (normalAtNode.x[1]/len);
				C(k+2*nnod,s) *= (normalAtNode.x[2]/len);
			}

			ItNodeBCs it;
			for (it = m_Nodes[k]->getFirstBC(); it != m_Nodes[k]->getLastBC(); it++)
			{
				ptrBC = dynamic_cast<cBoundaryConditionStructure *>(it->second);
				if (ptrBC != NULL)
				{
					// BC fixed in x1
					if (ptrBC->checkIfFixed( disp_x1 )) 
					{
						for (int s=0; s<nnod; s++)
							C(k,       s) = 0.0;
					}   
					// BC fixed in x2
					if (ptrBC->checkIfFixed( disp_x2 )) 
					{
						for (int s=0; s<nnod; s++)
							C(k+nnod,  s) = 0.0;
					}
					// BC fixed in x3
					if (ptrBC->checkIfFixed( disp_x3 )) 
					{
						for (int s=0; s<nnod; s++)
							C(k+2*nnod,s) = 0.0;
					}
				}
			}

		}  // end loop over nodes for scaling and boundary conditions
		ptrBC = NULL;

		// --- add C to global system matrix
		MatSetValues(K, 3*nnod, &(rows[0]), nnod, &(cols[0]), &(C(0,0)), ADD_VALUES);

		/*  // output for testing purposes
		for (int z=0; z<3*nnod; z++)
		{
		for (int s=0; s<nnod; s++)
		std::cout << C(z,s) << "   ";
		std::cout << "\n";
		}
		std::cout<<"\n\n";
		*/
#else
	throw cException("UNABLE TO HANDLE THIS TYPE OF COUPLING ELEMENTS - RECOMPILE THE CODE INCLUDING SWEBEM", __FILE__, __LINE__);
#endif
}

PetscReal cElementCoupling::getAreaOfCouplingElement(void)
{
#ifdef HAVE_SWEBEM
	PetscInt nnod = getNumberOfNodes();
	PetscReal AreaElem;
        // get area of element 
	if (nnod == 3 || nnod == 6)
	{
		cVector ab(3), ac(3), cross(3);
		for (int k=0; k<3; k++)
		{
			ab[k] = (*m_Nodes[1])[k] - (*m_Nodes[0])[k];
			ac[k] = (*m_Nodes[2])[k] - (*m_Nodes[0])[k];
		}

		cross[0] =  ab[1]*ac[2]-ab[2]*ac[1];
		cross[1] = -ab[0]*ac[2]+ab[2]*ac[0];
		cross[2] =  ab[0]*ac[1]-ab[1]*ac[0];

		AreaElem = 0.5 * cross.abs();
	}
	else 
		if (nnod == 4 || nnod == 9)
		{
			cVector ab(3), ad(3), cd(3), cb(3), cross1(3), cross2(3);
			for (int k=0; k<3; k++)
			{
				ab[k] = (*m_Nodes[1])[k] - (*m_Nodes[0])[k];
				ad[k] = (*m_Nodes[3])[k] - (*m_Nodes[0])[k];
				cd[k] = (*m_Nodes[3])[k] - (*m_Nodes[2])[k];
				cb[k] = (*m_Nodes[1])[k] - (*m_Nodes[2])[k];
			}

			cross1[0] =  ab[1]*ad[2]-ab[2]*ad[1];
			cross1[1] = -ab[0]*ad[2]+ab[2]*ad[0];
			cross1[2] =  ab[0]*ad[1]-ab[1]*ad[0];

			cross2[0] =  cd[1]*cb[2]-cd[2]*cb[1];
			cross2[1] = -cd[0]*cb[2]+cd[2]*cb[0];
			cross2[2] =  cd[0]*cb[1]-cd[1]*cb[0]; 

			AreaElem = 0.5*(cross1.abs() + cross2.abs());
		}
		else
		{
			trace("  Coupling Element messy\n. ");
			PetscPrintf(PETSC_COMM_SELF, "Coupling element messy!\n");
			ExitApp();
		}
	return AreaElem;

#else
	throw cException("UNABLE TO GET AREA OF COUPLING ELEMENTS - RECOMPILE THE CODE INCLUDING SWEBEM", __FILE__, __LINE__);
#endif
}



std::ostream& cElementCoupling::write(std::ostream &os) const
{
	os << "interfacelement FEM<->BEM" << std::endl;
	os << "   id = " << getId() << std::endl;
	os << "   ori = " << getOrientation() << std::endl;
	for (int k=0; k<getNumberOfNodes(); k++)
		os << *m_Nodes[k] << std::endl;

	for (int k=0; k<getNumberOfNodes(); k++)
		os << *m_MatchingBemNodes[k] << std::endl;

	return os;
}


std::ostream& cElementCoupling::writeXml(std::ostream &os) const
{
	os << "<InterfaceSwebem n=\"" << getNumberOfNodes() << "\">";
	os << std::endl;
	os << "<Id>" << getId() << "</Id>";
	os << "<ori>" << std::showpoint << getOrientation() << "</ori>" << std::endl;
	for (int k=0; k<getNumberOfNodes(); k++)
		os << "<NF>" << m_Nodes[k]->getId() << "</NF>";
	os << std::endl;
	for (int k=0; k<getNumberOfNodes(); k++)
		os << "<NS>" << getIdMatchingBemNode(k) << "</NS>";
	os << std::endl;
	os << "</InterfaceSwebem>";
	return os;
}
