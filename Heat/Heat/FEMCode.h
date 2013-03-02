

#include <vector>
#include "FEMObjects.h"
#include <boost\numeric\ublas\matrix_sparse.hpp>
#include <fstream>
#include <boost\numeric\odeint.hpp>
#include <string.h>
#include <boost\asio.hpp>
#include <boost\thread.hpp>


void FEMCode(std::vector<elementMesh_ptr>); // the base code computing the temperature distribution

//void compose_FEM_mesh(std::vector<elementMesh_ptr>&, elementFEM_ptr_vector&, facetFEM_ptr_vector&); // given a collection of mesh elements, produce a FEM mesh

class ODE_System; // the ODE system manager

void calc_element_matrices(elementFEM&, ODE_System&); // calculates an element's C,K,Q

void impose_Neumann(facetFEM_ptr_vector&, ODE_System&); // calculates all the additional matrices and incorporates them into the global matrices;

void impose_Dirichlet(facetFEM_ptr_vector&, ODE_System&); // does entries elimination in the global K,Q; WORK ON AI,BI,RI IS STILL REQUIRED !!!

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

class ODE_System
{ // contains all necessary info of and integrates the system [C]{dT/dt} = -[K]{T} + {Q}
public:
	
	struct GlobalMatrices
	{
		boost::numeric::ublas::compressed_matrix<mx_elem> C, Cinv, K, K_Dir, K_Neu;
		vector Q, Q_Dir, Q_Neu;
		size_t NofNds; //may not be required
		
		GlobalMatrices() {}
		GlobalMatrices(size_t NumOfNds)
		{
			this->NofNds = NumOfNds;
			this->C(NumOfNds,NumOfNds); this->Cinv(NumOfNds,NumOfNds);
			this->K(NumOfNds,NumOfNds); this->Q(NumOfNds);	
			this->K_Dir(NumOfNds,NumOfNds); this->Q_Dir(NumOfNds);
			this->K_Neu(NumOfNds,NumOfNds); this->Q_Neu(NumOfNds);
		}
		~GlobalMatrices() {}
	
	} GM; // the ODE system matrices
	struct ODEParameters
	{
		GlobalMatrices* gm;
		double time_step;
		double time_span[2];
		std::string output_file;
		
		ODEParameters() {}
		ODEParameters(GlobalMatrices* in)
		{
			gm = in;
		}
		
		void Set(double dt, double TS[], std::string name)
		{			
			time_step = dt; output_file = name;
			for (index i = 0; i != 2; i++)
				time_span[i] = TS[i];
		}
	} parameters; // the parameters of integration

	ODE_System(size_t numOfnds)
	{
		this->GM = GlobalMatrices(numOfnds);
		rhs = RHS(this);
		parameters = ODEParameters(&(this->GM));	
	}

	bool integrate(boost::numeric::ublas::vector<temperature>& initials);
	
private:
	
	struct RHS
	{	
		ODE_System* host;

		RHS() {}
		RHS(ODE_System* in)
		{
			this->host = in;
			bool inverted = host->InvertMatrix(host->GM.C, host->GM.Cinv);
		}

		void operator() (vector &T, vector &dTdt, double t)
		{
			dTdt = prod(host->GM.Cinv, (-prod(host->GM.K, T) + host->GM.Q));
		}
	} rhs; // the right-hand side of the system of ODEs
	struct Observer
	{
		std::ofstream &output;
		
		Observer(std::ofstream out) : output(out) {}
		~Observer() { output.close(); delete output;}
		
		void operator() (const vector &x, double t) const
		{
			output << t;
			for( size_t i = 0 ; i < x.size() ; i++ )
				output << "\t" << x[i];
			output << "\n";
		}
	}; 
		// the object outputting the results of the integration of the ODE system
	
	bool InvertMatrix(const matrix& input, boost::numeric::ublas::compressed_matrix<mx_elem>& inverse);
};

class FEM_Model // will implement boundary conditions imposition
{
public:
	elementFEM_ptr_vector			 HexElement;
	facetFEM_ptr_vector				 Boundary;
	
	FEM_Model() {}
	FEM_Model(std::vector<elementMesh_ptr>&);	
	~FEM_Model() { HexElement.clear(); Boundary.clear(); }
};

class Thread_Pool
{
public:

	Thread_Pool(size_t);
	template<class F>
	void enqueue(F task);
	~Thread_Pool();

private:
	boost::asio::io_service							service;
	boost::asio::io_service::work					work;	
	std::vector<std::unique_ptr<boost::thread>>		threads;

	class Worker
	{
	public:
		Worker(Thread_Pool &one) : pool(one) {}
		void operator()();
	private:
		Thread_Pool &pool;
	};
};

