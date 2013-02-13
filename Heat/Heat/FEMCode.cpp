#include "stdafx.h"
#include "FEMObjects.h"
#include "FEMCode.h"

// #define _SCL_SECURE_NO_WARNINGS ! is defined in properties -> preprocessor definitions

void FEMCode(std::vector<elementMesh> elements)
{	
	// the following sections should probably be put into separate functions

	// FEM mesh and boundary generation
	elementFEM_ptr_vector	HexElement(elements.size());
	facetFEM_ptr_vector		Boundary(6*HexElement.size());
	compose_FEM_mesh(elements,HexElement,Boundary);
	elements.~vector();
	
	// ODE system obtaining		
	GlobalMatrices gm(1000); // the actual number of nodes put here instead of 1000, finding ways to avoid preallocating it is preferable	
	BOOST_FOREACH(elementFEM element, HexElement)				
		calc_element_matrices(&element, gm);

	// the imposition of Neumann boundary conditions
	impose_Neumann(Boundary, gm);

	// ODE system solution
}


void compose_FEM_mesh(std::vector<elementMesh> &elements, elementFEM_ptr_vector &HexElements, facetFEM_ptr_vector &Boundaries)
{	
	elementFEM::NaturalCoordinates_Set = false;
	
	index cnt = 0; // may be calculated **
	BOOST_FOREACH(elementMesh elem, elements){	
		elementFEM one = elementFEM(elem,cnt);
		HexElements.push_back(&one);
		BOOST_FOREACH(facetFEM f, one.Facet)
		{
			index nind = 0; size_t NofFacetNds = f.Node.size(); size_t NofElts = 1;			
			while((NofElts == 1) && (nind < NofFacetNds)){
				NofElts = f.Node.at(nind).element.size(); // the number of elements that reference a node of a facet
				nind++;
			}

			if(NofElts == 1) // if all nodes of the facet are referenced by only one element then so is the facet, and hence it's boundary
				Boundaries.push_back(&f);		
		}
		cnt++; // **, for now it's just increment
	}

}

void calc_element_matrices(elementFEM_ptr element, GlobalMatrices &Global)
{
	size_t size = element->Node.size();
	std::vector<index> sctr = std::vector<index>(size);	// for assembling the current matrices into the global ones	
	for(index n = 0; n < size; n++)
		sctr[n] = element->Node.at(n).iGlob;

	element->Matrix.C = element->Matrix.calculate_C(); // save the matrices inside the element for no reason
	element->Matrix.K = element->Matrix.calculate_K();
	element->Matrix.Q = element->Matrix.calculate_Q();

	for (index row = 0; row < size; row++) // the assembling
	{
		Global.Q(sctr[row]) += element->Matrix.Q(row);
		for (index col = 0; col < size; col++)
		{
			Global.C(sctr[row],sctr[col]) += element->Matrix.C(row,col);
			Global.K(sctr[row],sctr[col]) += element->Matrix.K(row,col);			
		}
	}
}

void impose_Neumann(facetFEM_ptr_vector &NeuBdry, GlobalMatrices &Global){
	// the matrices after the imposition are based upon the free ones **
	Global.K_Neu = Global.K; Global.Q_Neu = Global.Q;
	
	BOOST_FOREACH(facetFEM f, NeuBdry){
		size_t size = f.Node.size();
		std::vector<index> sctr = std::vector<index>(size);	// for assembling the current matrices into the global ones		
		for(index n = 0; n < size; n++)
			sctr[n] = f.Node.at(n).iGlob;

		f.K_Neu = f.calc_K_Neu();
		f.Q_Neu = f.calc_Q_Neu();

		for (index row = 0; row < size; row++) // the assembling
		{
			Global.Q_Neu(sctr[row]) += f.Q_Neu(row);
			for (index col = 0; col < size; col++)							
				Global.K_Neu(sctr[row],sctr[col]) += f.K_Neu(row,col);		// with additional terms **			
		}

	}

}

//void impose_Dirichlet(std::vector<facetFEM_ptr> &DirBdry, GlobalMatrices &Global){
//	// if the Natural boundary is also imposed then it goes first and the Essential boundary is applied to K_Neu, Q_Neu; 
//	// for now suppose Neumann and Dirichlet are separate 
//	Global.K_Dir = Global.K; Global.Q_Dir = Global.Q;
//	
//	// manage indexes
//	size_t AIs = Global.K_Dir.size2(); // All Indexes size
//	std::vector<index> AI = std::vector<index>(AIs); 	
//	for (index i = 0; i < AIs; i++) //fill the all indexes vector
//		AI[i] = i;
//	std::vector<index> BI = std::vector<index>(AIs); // Boundary Indexes
//	std::vector<temperature> TI = std::vector<temperature>(AIs); // Temperature at boundary Indexes
//	BOOST_FOREACH(elementFEM::facetFEM* f, DirBdry){
//		BOOST_FOREACH(node* n, f->Node){
//			BI.push_back(n->iGlob);
//			TI.push_back(f->T_dir);			
//		}
//	}
//	BI.shrink_to_fit(); TI.shrink_to_fit();	size_t BIs = BI.size();
//
//	std::vector<index> BIcopy = std::vector<index>(BI); // a copy of BI to not mess up the correlation between BI and TI
//	std::sort(BIcopy.begin(),BIcopy.end()); // the sorted copy is used to construct RI
//
//	std::vector<index> RI = std::vector<index>(AIs); // the Remaining Indexes
//	RI = std::set_difference(AI.begin(),AI.end(),BIcopy.begin(),BIcopy.end(),RI); // RI := AI\BI
//	BIcopy.~vector(); AI.~vector(); // are no longer necessary
//	RI.shrink_to_fit(); size_t RIs = RI.size();
//
//	// now to the elemination procedure necessary for the imposition	
//	index ri, bi;
//	for (index s = 0; s < RIs; s++)
//	{
//		ri = RI[s];
//		for (index i = 0; i < BIs; i++)
//		{
//			bi = BI[i];
//			Global.Q_Dir(bi) = Global.K_Dir(bi,bi)*TI[bi];
//			for (index j = 0; j < AIs; j++)		
//			{
//				if (bi != j)
//					Global.K_Dir(bi,j) = 0;
//			}
//			Global.Q_Dir(ri) -= Global.K_Dir(ri,bi)*TI[bi];
//			Global.K_Dir(ri,bi) = 0;
//		}
//	}
//}

GlobalMatrices::GlobalMatrices(size_t NumOfNds)
{
	this->NofNds = NumOfNds;
	this->C(NumOfNds,NumOfNds); this->K(NumOfNds,NumOfNds); this->Q(NumOfNds);	
	this->K_Dir(NumOfNds,NumOfNds); this->Q_Dir(NumOfNds);
	this->K_Neu(NumOfNds,NumOfNds); this->Q_Neu(NumOfNds);
}

GlobalMatrices::~GlobalMatrices()
{
}
