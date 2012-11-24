#ifndef FEMOBJECTS_H
#define FEMOBJECTS_H

#include <vector>
#include <boost\numeric\ublas\matrix.hpp>
#include <boost\numeric\ublas\symmetric.hpp>
#include <boost\foreach.hpp>

typedef	double			coordinate;
typedef	std::size_t		index;
typedef double			mx_elem;

typedef boost::numeric::ublas::matrix<mx_elem>				matrix;
typedef boost::numeric::ublas::symmetric_matrix<mx_elem>	sym_matrix;
typedef boost::numeric::ublas::vector<mx_elem>				vector;


class elementFEM;

class node
{	
public:
	coordinate x[3]; // for the present assume they're given explicitly, and not via the i,j,k indexes on the grid
	index iLoc, iGlob; // an element's local and global indexes; may be set at elementMesh mesh generation
	std::vector<elementFEM*> element;	// the elements a node belongs to
};

class material;

// what is passed into my code to be attributed with FEM parameters like indexes
class elementMesh
{
public:
	material *Stuff;	// a pointer to the object containing the material properties of an element
	std::vector<node*>	Node;	// an array of the nodes of an element

	elementMesh();
	~elementMesh();

private:

};

class elementFEM : elementMesh // adds FEM - required properties to elementMesh
{
public:
	
	index iGlob; // this one may already be in elementMesh; Or rather may not be required. 
	// instead, might want to add some sort of boundary reference for the boundary conditions imposition
	static bool NaturalCoordinates_Set; // to avoid multiple initialization of NC
	
	// these properties are ODE System - related, the local matrices to be scattered into the system. For now place them here
	struct Matrix_Container
	{
		elementFEM *host;
		sym_matrix C, K;		
		vector Q;

		// functions to compute C, K, Q
		sym_matrix calculate_C();
		sym_matrix calculate_K();			
		vector calculate_Q();
	};

	Matrix_Container Matrix; 					
			
	struct Shape
	{	
		// might want to change this one to boost::vector to remove boost\assign\std\vector in the cpp, gives the same functionality (by the looks of it)
		// but doesn't require extra dependencies
		static std::vector<std::vector<coordinate>> NC;	// natural coordinates of the element's nodes (3 coordinates, 8 nodes);	
		elementFEM *host;

		void set_form_functions(); // for now only sets NC; for other sorts of functions might add something else	
		mx_elem	function_i(coordinate,coordinate,coordinate,index);	
		vector gradn_function_i(coordinate,coordinate,coordinate,index); // changes with the concrete type of form functions
		coordinate transform(std::vector<coordinate>,coordinate,coordinate,coordinate);
		matrix Jacobian(coordinate,coordinate,coordinate); 
	};
	
	Shape form;

	elementFEM(elementMesh, int);
	~elementFEM();	
	
};


class material
{
public:

	matrix Conduct; // thermal conductivity of material

	material();
	~material();

private:

};


struct facet
{
	std::vector<node> Node; // already flattened
	sym_matrix K_Neu;
	vector Q_Neu;

	facet();
	~facet();
	
	sym_matrix calc_K_Neu();
	vector calc_Q_Neu();	

private:

	void flatten(); // need it for boundary condition imposition

	mx_elem function_i(coordinate,coordinate,index);
	vector gradn_function_i(coordinate,coordinate,index); // gradient in natural coordinates. grad_xyz = Jacobian^-1 * gardn
	matrix Jacobian(coordinate,coordinate); // 2D jacobian on the surface of an element

};

#endif


