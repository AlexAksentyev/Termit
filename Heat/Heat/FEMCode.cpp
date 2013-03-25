#include "stdafx.h"
#include "FEMCode.h"
#include <boost\bind.hpp>
#include <boost\ref.hpp>
#include <crtdbg.h>


// #define _SCL_SECURE_NO_WARNINGS ! is defined in properties -> preprocessor definitions

void FEMCode(std::vector<elementMesh_ptr> elements) // instead of the vector of elements, indices of the elements will be transfered into the code; the (full) mesh model is transfered only first time
{	
	// the following sections should probably be put into separate functions

	// the FEM mesh and boundary generation
	FEM_Model Model(elements);
	elements.~vector(); // destroys the elements of the vector but keeps the referenced mesh model elements untouched for further use (i've no right to destroy them)
	
	// the ODE system obtaining		
	ODE_System System(1000); // the actual number of nodes put here instead of 1000, FINDING WAYS TO AVOID PREALLOCATING IT IS PREFERABLE !!!
	Thread_Pool Pool(boost::thread::hardware_concurrency()-1);
	BOOST_FOREACH(elementFEM element, Model.HexElement)				
		Pool.enqueue( boost::bind(calc_element_matrices,boost::ref(element), boost::ref(System)) ); 
	
	// the imposition of Neumann boundary conditions
	BOOST_FOREACH(facetFEM f, Model.Boundary)
		Pool.enqueue( boost::bind(impose_Neumann,boost::ref(f), boost::ref(System)) ); 
	Pool.~Thread_Pool(); // close the pool

	// the ODE system solution
	double step = 1e-3; double span [2] = {0, 10}; std::string filename = std::string("C:\output");
	System.parameters.Set(step, span, filename);

	typedef boost::numeric::ublas::vector<temperature>								Tvector;
	Tvector Tini = Tvector(System.GM.K.size2()); // initializes it with zeros
	
	System.integrate(Tini);	

	_CrtDumpMemoryLeaks();
}

//void compose_FEM_mesh(std::vector<elementMesh_ptr> &elements, elementFEM_ptr_vector &HexElements, facetFEM_ptr_vector &Boundaries)
//{	
//	elementFEM::NaturalCoordinates_Set = false;
//	
//	index cnt = 0; // may (and should) be calculated **
//	BOOST_FOREACH(elementMesh_ptr elem_ptr, elements){			
//		elementFEM_ptr one = new elementFEM(elem_ptr,cnt);
//		HexElements.push_back(one);
//		BOOST_FOREACH(facetFEM f, one->Facet)
//		{
//			index nind = 0; size_t NofFacetNds = f.Node.size(); size_t NofElts = 1;			
//			while((NofElts == 1) && (nind < NofFacetNds)){
//				NofElts = f.Node.at(nind)->element.size(); // the number of elements that reference a node of a facet
//				nind++;
//			}
//
//			if(NofElts == 1) // if all nodes of the facet are referenced by only one element then so is the facet, and hence it's boundary
//				Boundaries.push_back(&f);		
//		}
//		cnt++; // **, for now it's just increment
//	}
//
//}

void calc_element_matrices(elementFEM &element, ODE_System &Sys) 
{
	size_t size = element.Node.size();
	std::vector<index> sctr = std::vector<index>(size);	// for assembling the current matrices into the global ones	
	for(index n = 0; n < size; n++)
		sctr[n] = element.Node.at(n)->iGlob;

	element.Matrix.C = element.Matrix.calculate_C(); // save the matrices inside the element for no reason
	element.Matrix.K = element.Matrix.calculate_K();
	element.Matrix.Q = element.Matrix.calculate_Q();

	for (index row = 0; row < size; row++) // the assembling
	{
		Sys.GM.Q(sctr[row]) += element.Matrix.Q(row);
		for (index col = 0; col < size; col++)
		{
			Sys.GM.C(sctr[row],sctr[col]) += element.Matrix.C(row,col);
			Sys.GM.K(sctr[row],sctr[col]) += element.Matrix.K(row,col);			
		}
	}
}

void impose_Neumann(facetFEM &face, ODE_System &Sys){
	// the matrices after the imposition are based upon the free ones **
	Sys.GM.K_Neu = Sys.GM.K; Sys.GM.Q_Neu = Sys.GM.Q;
	
	size_t size = face.Node.size();
	std::vector<index> sctr = std::vector<index>(size);	// for assembling the current matrices into the global ones		
	for(index n = 0; n < size; n++)
		sctr[n] = face.Node.at(n)->iGlob;

	face.K_Neu = face.calc_K_Neu();
	face.Q_Neu = face.calc_Q_Neu();

	for (index row = 0; row < size; row++) // the assembling
	{
		Sys.GM.Q_Neu(sctr[row]) += face.Q_Neu(row);
		for (index col = 0; col < size; col++)							
			Sys.GM.K_Neu(sctr[row],sctr[col]) += face.K_Neu(row,col);		// with additional terms **			
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

FEM_Model::FEM_Model(std::vector<elementMesh_ptr>& elements)
{
	elementFEM::NaturalCoordinates_Set = false;

	HexElement = elementFEM_ptr_vector(elements.size());
	// preallocation of Boundary is desirable but i've yet no idea of how many boundary elements there are
	
	index cnt = 0; // may (and should) be calculated **
	BOOST_FOREACH(elementMesh_ptr elem_ptr, elements) // SHOULD BE PARALLELED !!!!
	{			
		elementFEM_ptr one = new elementFEM(elem_ptr,cnt);
		HexElement.push_back(one);
		BOOST_FOREACH(facetFEM f, one->Facet)
		{
			index nind = 0; size_t NofFacetNds = f.Node.size(); size_t NofElts = 1;			
			while((NofElts == 1) && (nind < NofFacetNds)){
				NofElts = f.Node.at(nind)->element.size(); // the number of elements that reference a node of a facet
				nind++;
			}

			if(NofElts == 1) // if all nodes of the facet are referenced by only one element then so is the facet, and hence it's boundary
				Boundary.push_back(&f);		
		}
		cnt++; // **, for now it's just increment
	}

}

void Thread_Pool::Worker::operator()() { pool.service.run(); }

Thread_Pool::Thread_Pool(size_t NumOfThrds) : work(service)
{
	for(index t = 0; t != NumOfThrds; t++)
		threads.push_back(std::unique_ptr<boost::thread>(new boost::thread(Worker(*this))));
}

Thread_Pool::~Thread_Pool()
{
	service.stop();
	for(size_t t = 0; t != threads.size(); t++)
		threads[t]->join();
}

template<class F>
void Thread_Pool::enqueue(F task) { service.post(task); }