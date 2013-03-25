#ifndef ODESYSTEM_H
#define ODESYSTEM_H

#include <boost\numeric\ublas\matrix_sparse.hpp>
#include "FEMObjects.h"
#include <fstream>
#include <string>


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

#endif