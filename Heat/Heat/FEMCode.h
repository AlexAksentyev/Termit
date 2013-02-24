#include <vector>
#include "FEMObjects.h"
#include <boost\numeric\ublas\matrix_sparse.hpp>
#include <fstream>
#include <boost\function.hpp>
#include <boost\numeric\odeint.hpp>


void FEMCode(std::vector<elementMesh_ptr>); // the base code computing the temperature distribution

void compose_FEM_mesh(std::vector<elementMesh_ptr>&, elementFEM_ptr_vector&, facetFEM_ptr_vector&); // given a collection of mesh elements, produce a FEM mesh

struct GlobalMatrices;

void calc_element_matrices(elementFEM&, GlobalMatrices&); // calculates an element's C,K,Q

void impose_Neumann(facetFEM_ptr_vector&, GlobalMatrices&); // calculates all the additional matrices and incorporates them into the global matrices;

void impose_Dirichlet(facetFEM_ptr_vector&, GlobalMatrices&); // does entries elimination in the global K,Q; WORK ON AI,BI,RI IS STILL REQUIRED !!!

void ODE_Solver(GlobalMatrices&, boost::numeric::ublas::vector<temperature>, double); // solves the final system of ordinary differential equations

bool InvertMatrix(const matrix& input, matrix& inverse)
{
	typedef boost::numeric::ublas::permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input
	matrix A(input);

	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;

	// create identity matrix of "inverse"
	inverse.assign(boost::numeric::ublas::identity_matrix<mx_elem> (A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}

struct GlobalMatrices
{
	boost::numeric::ublas::compressed_matrix<mx_elem> C, K, K_Dir, K_Neu;
	vector Q, Q_Dir, Q_Neu;
	size_t NofNds; //may not be required

	GlobalMatrices();
	GlobalMatrices(size_t);
	~GlobalMatrices();

};

class RHS
{
	GlobalMatrices GM;
	matrix Cinv;
	
public:

	RHS(GlobalMatrices&);

	void operator() (vector&, vector&, double);
};

struct observer
{
	std::ofstream &output;

	observer(std::ofstream out) : output(out) {}

	void operator() (const vector &x, double t) const
	{
		output << t;
          for( size_t i = 0 ; i < x.size() ; i++ )
              output << "\t" << x[i];
          output << "\n";
	}
};