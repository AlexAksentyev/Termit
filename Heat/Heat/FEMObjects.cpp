#include "stdafx.h"
#include "FEMObjects.h"
#include <boost\assign\std\vector.hpp>
#include "MatrixIntegration.h"

#define ff_num 8


///////////////// NODE ////////////////////////
node::~node()
{	
	this->element.clear();
}


///////////////// ELEMENT_FEM //////////////////

bool elementFEM::NaturalCoordinates_Set;
std::vector<std::vector<coordinate>> elementFEM::Shape::NC;

elementFEM::elementFEM(elementMesh_ptr one, index number)
{
	this->iGlob = number;
	this->Matrix.host = this; this->form.host = this; // might want to change these to boost::shared_pointer
	this->Stuff = one->Stuff; this->Node = one->Node;

	// there's also some facet-related info in a mesh element which translates into facetFem info in an elementFEM element
	BOOST_FOREACH(facet_ptr f, one->Facet)
		this->Facet.push_back(new facetFEM(f));

	for(index n = 0; n != this->Node.size(); n++)
		this->Node[n]->element.push_back(this);	

	if (this->NaturalCoordinates_Set != true)
		this->form.set_form_functions(); // doesn't really set the functions themselves, but rather initializes the necessities	

	this->Matrix.C(ff_num,ff_num); this->Matrix.K(ff_num,ff_num); this->Matrix.Q(ff_num);
}

elementFEM::~elementFEM()
{	
	this->Facet.clear(); // facetFEMs are created on constuction, the rest of the arrays used in a FEM element are just taken from the mesh	
}

// to be changed according to the preferred form functions
void elementFEM::Shape::set_form_functions()
{
	using namespace boost::assign;
	this->NC = std::vector<std::vector<coordinate>>(3);
	NC[0] += -1, -1,  1, -1, -1, -1, 1, -1;
	NC[1] += -1,  1,  1,  1, -1,  1, 1,  1;
	NC[2] += -1, -1, -1, -1,  1,  1, 1,  1;

	elementFEM::NaturalCoordinates_Set = true;
}

// 8 - node type of element assumed
mx_elem elementFEM::Shape::function_i(coordinate u,coordinate v,coordinate g,index inode)
{
	mx_elem f = 1/8 * (1 + u*this->NC[0][inode])*(1 + v*this->NC[1][inode])*(1 + g*this->NC[2][inode]);
	
	return f;
}

vector elementFEM::Shape::gradn_function_i(coordinate u,coordinate v,coordinate g,index inode)
{
	coordinate ui = this->NC[0][inode]; coordinate vi = this->NC[1][inode]; coordinate gi = this->NC[2][inode];
	vector gradn(3); 
	gradn[0] = ui*(1 + v*vi)*(1 + g*gi); gradn[1] = vi*(1 + u*ui)*(1 + g*gi); gradn[2] = gi*(1 + u*ui)*(1 + v*vi);
	return 1/8 * gradn;
}

matrix elementFEM::Shape::Jacobian(coordinate u,coordinate v,coordinate g)
{	
	std::vector<vector> gradn_fs(ff_num);	
	for (index col = 0; col < ff_num; col++)
		gradn_fs[col] = this->gradn_function_i(u,v,g,col);			
	
	matrix gradn_f_mx (3,ff_num), XYZ(ff_num,3);
	for (index row = 0; row < 3; row++)
		for (index col = 0; col < ff_num; col++)
		{
			gradn_f_mx(row,col) = gradn_fs[col][row];
			XYZ(col,row) += this->host->Node[col]->x[row];
		}	

	matrix J(3,3);	
	J = prod(gradn_f_mx,XYZ);

	return J;
}

coordinate elementFEM::Shape::transform(std::vector<coordinate> X,coordinate u,coordinate v,coordinate g)
{
	size_t nn = X.size();	
	coordinate result = 0;

	for(index i = 0; i < nn; i++)
		result += X[i]*this->function_i(u,v,g,i);	

	return result;
}

sym_matrix elementFEM::Matrix_Container::calculate_C()
{			
	// at present, the integrand doesn't use the material properties of its region
	struct local : integrable<elementFEM>
	{	
		local(elementFEM_ptr element,index row,index column) : integrable(element,row,column){}

		mx_elem operator()(std::vector<coordinate> in)
		{
			coordinate u = in[0]; coordinate v = in[1]; coordinate g = in[2];

			mx_elem detJ = abs(det_inv(host->form.Jacobian(u,v,g))); // !! the absolute value

			return host->form.function_i(u,v,g,this->row)*host->form.function_i(u,v,g,this->col)*detJ;
		}
	};	

	// use an element's transform structure to transform the natural coordinates into x,y,z
	// triple sum the integrand for each of row, column of the triangular matrix
	int_domain domain = int_domain(2);
	using namespace boost::assign;
	domain[0] += -1,-1,-1; domain[1] += 1,1,1;
	
	integrand C_integr;	
	sym_matrix C(ff_num,ff_num);
	for (index row = 0; row < ff_num; row++) 
		for (index col = 0; col <= row; col++) //no need to calculate all entries, the matrix is symmetric
		{			
			C_integr = local(this->host,row,col);
			C(row,col) = gaussian_tplquad(C_integr,domain);
		}

	return C;
}

sym_matrix elementFEM::Matrix_Container::calculate_K()
{		
	// at present, the integrand doesn't use the material properties of its region
	struct local : integrable<elementFEM>
	{	
		local(elementFEM_ptr element,index row,index column) : integrable(element,row,column){}

		mx_elem operator()(std::vector<coordinate> in)
		{
			coordinate u = in[0]; coordinate v = in[1]; coordinate g = in[2];

			matrix J = host->form.Jacobian(u,v,g);
			matrix invJ;
			mx_elem detJ = abs(det_inv(J,invJ,true)); // !! the absolute value			

			matrix D = this->host->Stuff->Conduct;
			
			mx_elem Ni = host->form.function_i(u,v,g,this->row); mx_elem Nj = host->form.function_i(u,v,g,this->col);			
			vector gradn_Ni = host->form.gradn_function_i(u,v,g,this->row);
			vector gradn_Nj = host->form.gradn_function_i(u,v,g,this->col);			

			vector v0 = prod(trans(invJ),gradn_Nj);
			vector v1 = prod(D,v0);
			vector v2 = prod(invJ,v1);
			
			return ( inner_prod(gradn_Ni,v2) + Ni*Nj )*detJ; // coefficient before Ni*Nj not there yet, the material props that is
		}
	};		

	int_domain domain = int_domain(2);
	using namespace boost::assign;
	domain[0] += -1,-1,-1; domain[1] += 1,1,1;

	integrand K_integr;	
	sym_matrix K(ff_num,ff_num);
	for (index row = 0; row < ff_num; row++) 
		for (index col = 0; col <= row; col++) //no need to calculate all entries, the matrix is symmetric
		{			
			K_integr = local(this->host,row,col);			
			K(row,col) = gaussian_tplquad(K_integr,domain);			
		}

		return K;
}

vector elementFEM::Matrix_Container::calculate_Q()
{
	// at present, the integrand doesn't use the material properties of its region
	struct local : integrable<elementFEM>
	{	
		local(elementFEM_ptr element,index row, int column = 1) : integrable(element,row,column){}

		mx_elem operator()(std::vector<coordinate> in)
		{
			coordinate u = in[0]; coordinate v = in[1]; coordinate g = in[2];
			
			mx_elem detJ = abs(det_inv(host->form.Jacobian(u,v,g))); // !! the absolute value

			return host->form.function_i(u,v,g,this->row)*detJ;
		}
	};	

	int_domain domain = int_domain(2);
	using namespace boost::assign;
	domain[0] += -1,-1,-1; domain[1] += 1,1,1;

	integrand Q_integr;
	vector Q(ff_num);
	for (index row = 0; row < ff_num; row++)
	{
		Q_integr = local(this->host,row,1);
		Q(row) = gaussian_tplquad(Q_integr,domain);
	}

	return Q;
}

//////////////////////// FACET_FEM //////////////////////////////

mx_elem facetFEM::function_i(coordinate u,coordinate v,index inode)
{
	mx_elem f = 1/4 * (1 + u*elementFEM::Shape::NC[0][inode])*(1 + v*elementFEM::Shape::NC[1][inode]);
	
	return f;
}

vector facetFEM::gradn_function_i(coordinate u,coordinate v,index inode)
{
	coordinate ui = elementFEM::Shape::NC[0][inode]; coordinate vi = elementFEM::Shape::NC[1][inode];
	vector gradn(2); 
	gradn[0] = ui*(1 + v*vi); gradn[1] = vi*(1 + u*ui);
	return 1/4 * gradn;
}

matrix facetFEM::Jacobian(coordinate u,coordinate v)
{
	size_t nds_on_fct = this->Node.size();
	std::vector<vector> gradn_fs(nds_on_fct);	
	for (index col = 0; col < nds_on_fct; col++)
		gradn_fs[col] = this->gradn_function_i(u,v,col);			
	
	matrix gradn_f_mx (2,nds_on_fct), XY(nds_on_fct,2);
	for (index row = 0; row < 2; row++)
		for (index col = 0; col < nds_on_fct; col++)
		{
			gradn_f_mx(row,col) = gradn_fs[col][row];
			XY(col,row) += this->Node[col]->x[row];
		}	

	matrix J(2,2);	
	J = prod(gradn_f_mx,XY);

	return J;
}

sym_matrix facetFEM::calc_K_Neu()
{
	struct local : integrable<facetFEM>
	{			
		local(facetFEM_ptr host,index row,index column) : integrable(host,row,column){}		

		mx_elem operator()(std::vector<coordinate> in)
		{
			coordinate u = in[0]; coordinate v = in[1]; 

			matrix J = this->host->Jacobian(u,v);			
			mx_elem detJ = abs(det_inv(J)); // !! the absolute value			
			mx_elem Ni = this->host->function_i(u,v,this->row);
			mx_elem Nj = this->host->function_i(u,v,this->col);
			
			return ( Ni*Nj )*detJ; // coefficient before Ni*Nj not there yet, the convection_constant that is
		}
	};		

	int_domain domain = int_domain(2);
	using namespace boost::assign;
	domain[0] += -1,-1; domain[1] += 1,1;

	integrand K_integr;	
	size_t nds_on_fct = this->Node.size();
	sym_matrix K(nds_on_fct,nds_on_fct);
	for (index row = 0; row < nds_on_fct; row++) 
		for (index col = 0; col <= row; col++) //no need to calculate all entries, the matrix is symmetric
		{			
			K_integr = local(this,row,col);			
			K(row,col) = gaussian_dblquad(K_integr,domain);			
		}

	return K;
}

vector facetFEM::calc_Q_Neu()
{
	struct local : integrable<facetFEM>
	{	
		local(facetFEM_ptr host,index row, int column = 1) : integrable(host,row,column){}

		mx_elem operator()(std::vector<coordinate> in)
		{
			coordinate u = in[0]; coordinate v = in[1];
			
			mx_elem detJ = abs(det_inv(host->Jacobian(u,v))); // !! the absolute value

			return host->function_i(u,v,this->row)*detJ; // coefficient before Ni not there yet, the convection_constant*Tambnt + Q_boundary that is
		}
	};	

	int_domain domain = int_domain(2);
	using namespace boost::assign;
	domain[0] += -1,-1; domain[1] += 1,1;

	integrand Q_integr;
	size_t nds_on_fct = this->Node.size();
	vector Q(nds_on_fct);
	for (index row = 0; row < nds_on_fct; row++)
	{
		Q_integr = local(this,row,1);
		Q(row) = gaussian_tplquad(Q_integr,domain);
	}

	return Q;
}

//void facetFEM::impose_Neumann(ODE_System &Sys)
//{
//	size_t size = this->Node.size();
//	std::vector<index> sctr = std::vector<index>(size);	// for assembling the current matrices into the global ones		
//	for(index n = 0; n < size; n++)
//		sctr[n] = this->Node.at(n)->iGlob;
//
//	this->K_Neu = this->calc_K_Neu();
//	this->Q_Neu = this->calc_Q_Neu();
//
//	for (index row = 0; row < size; row++) // the assembling
//	{
//		Sys.GM.Q_Neu(sctr[row]) += this->Q_Neu(row);
//		for (index col = 0; col < size; col++)							
//			Sys.GM.K_Neu(sctr[row],sctr[col]) += this->K_Neu(row,col);		// with additional terms **			
//	}
//}

facetFEM::facetFEM(facet &face) // makes a FEM facet out of a simple mesh facet
{
	// add Node, h_conv, Tamb extraction from face
	size_t size = this->Node.size();
	this->K_Neu(size, size); this->Q_Neu(size);

	this->Node = face.Node;
}

facetFEM::facetFEM(facet_ptr face) // makes a FEM facet out of a simple mesh facet
{
	// add Node, h_conv, Tamb extraction from face
	size_t size = this->Node.size();
	this->K_Neu(size, size); this->Q_Neu(size);

	this->Node = face->Node;
}

facetFEM::~facetFEM()
{
	// might want to erase the local matrices
}
//////////// MATERIAL ////////////

material::material()
{
}

material::material(material_ptr)
{
}

material::~material()
{
}
