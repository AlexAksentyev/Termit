#include "stdafx.h"
#include "FEMObjects.h"
#include "FEMCode.h"

// #define _SCL_SECURE_NO_WARNINGS ! is defined in properties -> preprocessor definitions

void FEMCode(std::vector<elementMesh> elements)
{	
	// FEM mesh generation
	std::vector<elementFEM> HexElement = compose_FEM_mesh(elements);
	
	// ODE system obtaining		
	//GlobalMatrices gm(1000); // the actual number of nodes put here instead of 1000	
	BOOST_FOREACH(elementFEM element, HexElement)			
		calc_element_matrices(element/*, gm*/);


	// ODE system solution
}


std::vector<elementFEM> compose_FEM_mesh(std::vector<elementMesh> elements)
{	
	elementFEM::NaturalCoordinates_Set = false;

	std::vector<elementFEM> HexElements;// = std::vector<elementFEM>(elements.size());
	
	int cnt = 0; // may be calculated **
	BOOST_FOREACH(elementMesh elem, elements){		
		HexElements.push_back(elementFEM(elem,cnt));
		cnt++; // **
	}

	return HexElements;
}

void calc_element_matrices(elementFEM &element/*, GlobalMatrices &Global*/)
{
	std::vector<size_t> sctr = std::vector<size_t>();
	/*for(size_t n = 0; n < element.Node.size(); n++){

	}*/
	element.Matrix.C = element.Matrix.calculate_C(); // save the matrices inside the element for no reason
	//Global.C = element.Matrix.C; // this is actually where they're needed
	element.Matrix.K = element.Matrix.calculate_K();
	//Global.K( = element.Matrix.K;
	element.Matrix.Q = element.Matrix.calculate_Q();
	//Global.Q = element.Matrix.Q;
}

GlobalMatrices::GlobalMatrices(size_t NumOfNds)
{
	this->C(NumOfNds,NumOfNds); this->K(NumOfNds,NumOfNds); this->Q(NumOfNds);
	this->K_Dir(NumOfNds,NumOfNds); this->Q_Dir(NumOfNds);
	this->K_Neu(NumOfNds,NumOfNds); this->Q_Neu(NumOfNds);
}

GlobalMatrices::~GlobalMatrices()
{
}
