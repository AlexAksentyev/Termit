#include "FEMObjects.h"
#include <boost\function.hpp>
#include <boost\assign\std\vector.hpp>
#include <boost\numeric\ublas\lu.hpp>

typedef boost::function<mx_elem (std::vector<coordinate>)> integrand;
typedef std::vector<std::vector<coordinate>> int_domain; // since the domain, in our rectangular case, is defined as [a0,b0]x[a1,b1]x...x[an,bn], int_domain is an std::vector<coordinate>[2]


template<class host_type>
struct integrable
	{		
		index row, col;	
		host_type *host;		

		integrable(host_type *host,index i, index j)
		{
			this->row = i; this->col = j; this->host = host;
		}

		virtual mx_elem operator() (std::vector<coordinate>) = 0;
	};

mx_elem gaussian_tplquad(integrand, int_domain); 
mx_elem gaussian_dblquad(integrand, int_domain);
mx_elem det_inv(const matrix&, matrix& inv = matrix(), bool calc_inv = false);
int determinant_sign(const boost::numeric::ublas::permutation_matrix<std ::size_t>&);
