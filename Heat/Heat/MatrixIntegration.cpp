#include "stdafx.h"
#include "MatrixIntegration.h"


mx_elem gaussian_tplquad(integrand f, int_domain X)
{		
	std::size_t nabsc = 2; // the order of the method (the number of knots where the function is evaluated)
	std::vector<mx_elem> weight(nabsc);	std::vector<coordinate> absc(nabsc); // order driven
	std::vector<coordinate> wdth(3), mid(3); // coordinate dimension driven
	using namespace boost::assign;
	weight += 1,1; absc += -0.5773502691896257,0.5773502691896257;
	coordinate vol = 1;
	for (index i = 0; i < 3; i++){
		wdth[i] = 0.5*(X[nabsc-1][i] - X[0][i]); mid[i] = 0.5*(X[nabsc-1][i] + X[0][i]);
		vol *= wdth[i];
	}

	mx_elem sum = 0; std::vector<coordinate> arg(3);
	for (index i = 0; i < nabsc; i++)
		for (index j = 0; j < nabsc; j++)
			for (index k = 0; k < nabsc; k++)
			{
				arg += mid[i] + wdth[i]*absc[i], mid[j] + wdth[j]*absc[j], mid[k] + wdth[k]*absc[k];
				sum += weight[i]*weight[j]*weight[k]*f(arg);
				arg.clear(); 
			}

	return vol*sum;
}

mx_elem gaussian_dblquad(integrand f, int_domain X)
{
	std::size_t nabsc = 2; // the order of the method (the number of knots where the function is evaluated)
	std::vector<mx_elem> weight(nabsc);	std::vector<coordinate> absc(nabsc); // order driven
	std::vector<coordinate> wdth(2), mid(2); // coordinate dimension driven
	using namespace boost::assign;
	weight += 1,1; absc += -0.5773502691896257,0.5773502691896257;
	coordinate area = 1;
	for (index i = 0; i < 2; i++){
		wdth[i] = 0.5*(X[nabsc-1][i] - X[0][i]); mid[i] = 0.5*(X[nabsc-1][i] + X[0][i]);
		area *= wdth[i];
	}

	mx_elem sum = 0; std::vector<coordinate> arg(2);	
	for (index j = 0; j < nabsc; j++)
		for (index k = 0; k < nabsc; k++)
		{
			arg += mid[j] + wdth[j]*absc[j], mid[k] + wdth[k]*absc[k];
			sum += weight[j]*weight[k]*f(arg);
			arg.clear(); 
		}

	return area*sum;
}

mx_elem det_inv(const matrix & m, matrix & inverse, bool calc_inv)
{	
	matrix A(m); // create a working copy of the input

	using namespace boost::numeric::ublas;	
	permutation_matrix<std ::size_t> pm(A.size1());
	mx_elem det = 1.0;
    if( lu_factorize(A,pm) )
	{
        det = 0.0; calc_inv = false;
    } else 
	{
        for(index i = 0; i < A.size1(); i++)
            det *= A(i,i); // multiply by elements on diagonal

        det = det * determinant_sign( pm );
    }

	if (calc_inv){
		// create identity matrix of "inverse"
		inverse.assign(identity_matrix<mx_elem> (A.size1()));
		
		// backsubstitute to get the inverse
		lu_substitute(A, pm, inverse);
	}

    return det;
}

int determinant_sign(const boost::numeric::ublas::permutation_matrix<std ::size_t>& pm)
{
    int pm_sign=1;
    std::size_t size = pm.size();
    for (std::size_t i = 0; i < size; ++i)
        if (i != pm(i))
            pm_sign *= -1.0; // swap_rows would swap a pair of rows here, so we change sign
    return pm_sign;
}
 
//double determinant( boost::numeric::ublas::matrix<double>& m ) {
//    boost::numeric::ublas::permutation_matrix<std ::size_t> pm(m.size1());
//    double det = 1.0;
//    if( boost::numeric::ublas::lu_factorize(m,pm) ) {
//        det = 0.0;
//    } else {
//        for(int i = 0; i < m.size1(); i++)
//            det *= m(i,i); // multiply by elements on diagonal
//        det = det * determinant_sign( pm );
//    }
//    return det;
//}

//template<class T>
//bool InvertMatrix(const matrix<T>& input, matrix<T>& inverse)
//{
//	typedef permutation_matrix<std::size_t> pmatrix;
//
//	// create a working copy of the input
//	matrix<T> A(input);
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
//	inverse.assign(identity_matrix<T> (A.size1()));
//
//	// backsubstitute to get the inverse
//	lu_substitute(A, pm, inverse);
//
//	return true;
//}
 