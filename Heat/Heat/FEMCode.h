#include <vector>
#include "FEMObjects.h"



void FEMCode(std::vector<elementMesh>); // the base code computing the temperature distribution

std::vector<elementFEM> compose_FEM_mesh(std::vector<elementMesh>); // given a collection of mesh elements, produce a FEM mesh

struct GlobalMatrices;

void calc_element_matrices(elementFEM&/*, GlobalMatrices&*/); // calculates an element's C,K,Q

struct GlobalMatrices
{
	sym_matrix C, K, K_Dir, K_Neu;
	vector Q, Q_Dir, Q_Neu;

	GlobalMatrices(size_t);
	~GlobalMatrices();

};