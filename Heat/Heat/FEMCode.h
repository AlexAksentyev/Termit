#include <vector>
#include "FEMObjects.h"
#include <boost\numeric\ublas\matrix_sparse.hpp>


void FEMCode(std::vector<elementMesh_ptr>); // the base code computing the temperature distribution

void compose_FEM_mesh(std::vector<elementMesh_ptr>&, elementFEM_ptr_vector&, facetFEM_ptr_vector&); // given a collection of mesh elements, produce a FEM mesh

struct GlobalMatrices;

void calc_element_matrices(elementFEM&, GlobalMatrices&); // calculates an element's C,K,Q

void impose_Neumann(facetFEM_ptr_vector&, GlobalMatrices&); // calculates all the additional matrices and incorporates them into the global matrices;

void impose_Dirichlet(facetFEM_ptr_vector&, GlobalMatrices&); // does entries elimination in the global K,Q; WORK ON AI,BI,RI IS STILL REQUIRED !!!

struct GlobalMatrices
{
	boost::numeric::ublas::compressed_matrix<mx_elem> C, K, K_Dir, K_Neu;
	vector Q, Q_Dir, Q_Neu;
	size_t NofNds; //may not be required

	GlobalMatrices(size_t);
	~GlobalMatrices();

};