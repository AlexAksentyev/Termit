#include "stdafx.h"
#include "ElementMeshProvider.h"
#include <boost\assign\std\vector.hpp>



std::vector<elementMesh_ptr> provide_mesh()
{
	using namespace boost::assign;
	std::vector<coordinate> c(3);
	//these can be modified if necessary
	c += 0,0,0;	node n0 = node(c,0); c += 5,0,0; node n1 = node(c,1); c += 5,3,0; node n2 = node(c,2); c += 0,3,0; node n3 = node(c,3); 
	c += 10,0,0; node n5 = node(c,5); c += 10,3,0; node n6 = node(c,6); c += 5,8,0; node n7 = node(c,7); c += 0,8,0; node n8 = node(c,8); //lower plane 
	
	c += 0,0,7;	node n9 = node(c,9); c += 5,0,7; node n10 = node(c,10); c += 5,3,7; node n11 = node(c,11); c += 0,3,7; node n12 = node(c,12); 
	c += 10,0,7; node n13 = node(c,13); c += 10,3,7; node n14 = node(c,14); c += 5,8,7; node n15 = node(c,15); c += 0,8,7; node n16 = node(c,16); //upper plane

	material_ptr S = new material();

	//but don't touch these if it works
	std::vector<elementMesh_ptr> mesh = std::vector<elementMesh_ptr>(3);  node_ptr_vector fnodes = node_ptr_vector(8); // because (NOW) node_ptr_vector is an std::vector, can use boost::assign
	fnodes += new node(n0), new node(n1), new node(n2), new node(n3), new node(n9), new node(n10), new node(n11), new node(n12);
	mesh += new elementMesh(fnodes,S);// 1st cube
	fnodes += new node(n1), new node(n5), new node(n6), new node(n2), new node(n10), new node(n13), new node(n14), new node(n11);
	mesh += new elementMesh(fnodes,S);// 2nd cube
	fnodes += new node(n3), new node(n2), new node(n7), new node(n8), new node(n12), new node(n11), new node(n15), new node(n16);
	mesh += new elementMesh(fnodes,S);// 3rd cube

	return mesh;
}