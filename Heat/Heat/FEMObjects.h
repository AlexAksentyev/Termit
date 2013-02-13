#ifndef FEMOBJECTS_H
#define FEMOBJECTS_H

#include <vector>
#include <boost\numeric\ublas\matrix.hpp>
#include <boost\numeric\ublas\symmetric.hpp>
#include <boost\foreach.hpp>
#include <boost\ptr_container\ptr_vector.hpp>

typedef	double			coordinate;
typedef double			temperature;
typedef	std::size_t		index;
typedef double			mx_elem;

typedef boost::numeric::ublas::matrix<mx_elem>				matrix;
typedef boost::numeric::ublas::symmetric_matrix<mx_elem>	sym_matrix;
typedef boost::numeric::ublas::vector<mx_elem>				vector;


class elementFEM;

typedef elementFEM*									elementFEM_ptr;
typedef boost::ptr_vector<elementFEM>				elementFEM_ptr_vector;

class node
{	
public:
	coordinate x[3]; // for the present assume they're given explicitly, and not via the i,j,k indices on the grid
	index iGlob; 
	elementFEM_ptr_vector element;	// the elements a node belongs to

	~node();
};

typedef node*										node_ptr;
typedef boost::ptr_vector<node>						node_ptr_vector;

class material;

typedef material*									material_ptr;

struct facet
{
};

typedef boost::ptr_vector<facet>					facet_ptr_vector;

// what is passed into my code to be attributed with FEM parameters like indices
class elementMesh
{
public:
	material_ptr Stuff;	// a pointer to the object containing the material properties of an element
	node_ptr_vector	Node;	// an array of the nodes of an element, the position of a node in the array is the node's local index in the element
	facet_ptr_vector Facet; // the facets of an element, not indexed yet

	elementMesh();
	~elementMesh();

};

typedef elementMesh*								elementMesh_ptr;

struct facetFEM;

typedef facetFEM*									facetFEM_ptr;
typedef boost::ptr_vector<facetFEM>					facetFEM_ptr_vector;

class elementFEM // adds FEM - required properties to elementMesh
{
public:

	material_ptr Stuff;	// a pointer to the object containing the material properties of an element
	node_ptr_vector Node;	// an array of the nodes of an element, the position of a node in the array is the node's local index in the element
	
	index iGlob; 
	static bool NaturalCoordinates_Set; // to avoid multiple initialization of NC	

	// these properties are ODE System - related, the local matrices to be scattered into the system. For now place them here
	struct Matrix_Container
	{
		elementFEM_ptr host;
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
		elementFEM_ptr host;

		void set_form_functions(); // for now only sets NC; for other sorts of functions might add something else	
		mx_elem	function_i(coordinate,coordinate,coordinate,index);	
		vector gradn_function_i(coordinate,coordinate,coordinate,index); // changes with the concrete type of form functions
		coordinate transform(std::vector<coordinate>,coordinate,coordinate,coordinate);
		matrix Jacobian(coordinate,coordinate,coordinate); 
	};
	
	Shape form;

	facetFEM_ptr_vector Facet; // initialized in the elementFEM constructor (should be, not yet)

	elementFEM(elementMesh&, index);
	~elementFEM();	
	
};

struct facetFEM
{
	node_ptr_vector Node; // already flat	
	sym_matrix K_Neu;
	vector Q_Neu;
	
	temperature T_dir, T_amb;
	double h_conv;
	
	sym_matrix calc_K_Neu();
	vector calc_Q_Neu();	

	facetFEM(facet&);
	~facetFEM();

private:

	void flatten(); // need it(?) for boundary condition imposition
	
	mx_elem function_i(coordinate,coordinate,index);
	vector gradn_function_i(coordinate,coordinate,index); // gradient in natural coordinates. grad_xyz = Jacobian^-1 * gardn
	matrix Jacobian(coordinate,coordinate); // 2D jacobian on the surface of an element

};

class material
{
public:

	matrix Conduct; // thermal conductivity of material

	material();
	material(material_ptr);
	~material();

private:

};


#endif


