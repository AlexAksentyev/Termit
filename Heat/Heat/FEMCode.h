#include <vector>
#include "FEMObjects.h"
#include <boost\numeric\ublas\matrix_sparse.hpp>
#include <fstream>
#include <boost\function.hpp>
#include <boost\numeric\odeint.hpp>
#include <string.h>


void FEMCode(std::vector<elementMesh_ptr>); // the base code computing the temperature distribution

void compose_FEM_mesh(std::vector<elementMesh_ptr>&, elementFEM_ptr_vector&, facetFEM_ptr_vector&); // given a collection of mesh elements, produce a FEM mesh

struct GlobalMatrices; struct ODEParameters;

void calc_element_matrices(elementFEM&, ODE_Solver::GlobalMatrices&); // calculates an element's C,K,Q

void impose_Neumann(facetFEM_ptr_vector&, ODE_Solver::GlobalMatrices&); // calculates all the additional matrices and incorporates them into the global matrices;

void impose_Dirichlet(facetFEM_ptr_vector&, ODE_Solver::GlobalMatrices&); // does entries elimination in the global K,Q; WORK ON AI,BI,RI IS STILL REQUIRED !!!

//void ODE_Solver(boost::numeric::ublas::vector<temperature>&, ODEParameters&); // solves the final system of ordinary differential equations

//bool InvertMatrix(const matrix& input, matrix& inverse)
//{
//	typedef boost::numeric::ublas::permutation_matrix<std::size_t> pmatrix;
//
//	// create a working copy of the input
//	matrix A(input);
//
//	// create a permutation matrix for the LU-factorization
//	pmatrix pm(A.size1());
//
//	// perform LU-factorization
//	int res = lu_factorize(A, pm);
//	if (res != 0)
//		return false;
//
//	// create identity matrix of "inverse"
//	inverse.assign(boost::numeric::ublas::identity_matrix<mx_elem> (A.size1()));
//
//	// backsubstitute to get the inverse
//	lu_substitute(A, pm, inverse);
//
//	return true;
//}

//struct GlobalMatrices
//{
//	boost::numeric::ublas::compressed_matrix<mx_elem> C, K, K_Dir, K_Neu;
//	vector Q, Q_Dir, Q_Neu;
//	size_t NofNds; //may not be required
//
//	GlobalMatrices();
//	GlobalMatrices(size_t);
//	~GlobalMatrices();
//
//};

//struct ODEParameters
//{
//	double time_step;
//	GlobalMatrices* matrices;
//	double time_span[2];
//	std::string output_file;
//
//	ODEParameters(double dt, GlobalMatrices* GM, double TS[], std::string name)
//	{
//		time_step = dt; matrices = GM; output_file = name;
//		for (index i = 0; i != 2; i++)
//			time_span[i] = TS[i];
//	}
//};

//class RHS
//{
//	GlobalMatrices* GM;
//	matrix Cinv;
//	
//public:
//
//	RHS(GlobalMatrices*);
//
//	void operator() (vector&, vector&, double);
//};

//struct Observer
//{
//	std::ofstream &output;
//
//	Observer(std::ofstream out) : output(out) {}
//
//	void operator() (const vector &x, double t) const
//	{
//		output << t;
//          for( size_t i = 0 ; i < x.size() ; i++ )
//              output << "\t" << x[i];
//          output << "\n";
//	}
//};

class ODE_Solver
{
public:
	
	struct GlobalMatrices
	{
		boost::numeric::ublas::compressed_matrix<mx_elem> C, Cinv, K, K_Dir, K_Neu;
		vector Q, Q_Dir, Q_Neu;
		size_t NofNds; //may not be required
		
		GlobalMatrices();
		GlobalMatrices(size_t);
		~GlobalMatrices();
	
	} GM;
	struct ODEParameters
	{
		GlobalMatrices* gm;
		double time_step;
		double time_span[2];
		std::string output_file;
		
		ODEParameters();
		ODEParameters(double dt, double TS[], std::string name)
		{			
			time_step = dt; output_file = name;
			for (index i = 0; i != 2; i++)
				time_span[i] = TS[i];
		}
	} parameters;

	ODE_Solver(size_t);

	bool solve(boost::numeric::ublas::vector<temperature>&);
	
private:
	
	struct RHS
	{	
		GlobalMatrices* gm;

		RHS();		

		void operator() (vector&, vector&, double);
	} rhs;
	struct Observer
	{
		std::ofstream &output;
		
		Observer(std::ofstream out) : output(out) {}
		
		void operator() (const vector &x, double t) const
		{
			output << t;
			for( size_t i = 0 ; i < x.size() ; i++ )
				output << "\t" << x[i];
			output << "\n";
		}
	};
	
	bool InvertMatrix(const matrix& input, matrix& inverse);
};