/*
 * math.cpp
 *
 *  Created on: Dec 6, 2019
 *      Author: dmarce1
 */


#include "math.hpp"


bool matrix_inverse(const std::array<vect, NDIM> &A, std::array<vect, NDIM> &Ainv) {
#if(NDIM==1)
	if( A[0][0] == 0.0 ) {
		return false;
	}
	Ainv[0][0] = 1.0 / A[0][0];
#else
#if(NDIM==2)
	real det = A[0][0] * A[1][1] - A[1][0] * A[0][1];
	if (det == 0.0) {
		return false;
	}
	const real detinv = 1.0 / det;
	Ainv[0][0] = +A[1][1] * detinv;
	Ainv[0][1] = -A[0][1] * detinv;
	Ainv[1][0] = -A[1][0] * detinv;
	Ainv[1][1] = +A[0][0] * detinv;
#else
	real det = A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]);
	if( det == 0.0 ) {
		return false;
	}
	det /**/-= A[0][1] * (A[1][0] * A[2][2] - A[2][0] * A[1][2]);
	det /**/+= A[0][2] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);
	const real detinv = 1.0 / det;
	Ainv[0][0] = (-A[1][2] * A[2][1] + A[1][1] * A[2][2]) * detinv;
	Ainv[0][1] = (+A[0][2] * A[2][1] - A[0][1] * A[2][2]) * detinv;
	Ainv[0][2] = (-A[0][2] * A[1][1] + A[0][1] * A[1][2]) * detinv;
	Ainv[1][0] = (+A[1][2] * A[2][0] - A[1][0] * A[2][2]) * detinv;
	Ainv[1][1] = (-A[0][2] * A[2][0] + A[0][0] * A[2][2]) * detinv;
	Ainv[1][2] = (+A[0][2] * A[1][0] - A[0][0] * A[1][2]) * detinv;
	Ainv[2][0] = (-A[1][1] * A[2][0] + A[1][0] * A[2][1]) * detinv;
	Ainv[2][1] = (+A[0][1] * A[2][0] - A[0][0] * A[2][1]) * detinv;
	Ainv[2][2] = (-A[0][1] * A[1][0] + A[0][0] * A[1][1]) * detinv;
#endif
#endif
	return true;
}
