/*
 * math.hpp
 *
 *  Created on: Dec 6, 2019
 *      Author: dmarce1
 */

#ifndef SRC_MATH_HPP_
#define SRC_MATH_HPP_

#include "vect.hpp"

real bspline(real);
real W(real, real);
real W_norm(real h);
real dW_norm_dh(real);
bool matrix_inverse(const std::array<vect, NDIM>&, std::array<vect, NDIM> &A);
vect rotate_to(const vect &u, const vect &n);
vect rotate_from(const vect &u, vect n);
#if(NDIM==3)
vect cross( const vect& a, const vect& b);
real triple_product(const vect &a, const vect &b, const vect &c);
#endif
real condition_number(const std::array<vect, NDIM> &A, std::array<vect, NDIM> &Ainv);
real triple_product(const vect &a, const vect &b, const vect &c);

#if(NDIM==3)
inline vect cross( const vect& a, const vect& b) {
	vect c;
	c[0] = +a[1] * b[2] - a[2] * b[1];
	c[1] = -a[0] * b[2] + a[2] * b[0];
	c[2] = +a[2] * b[1] - a[1] * b[2];
	return c;
}

inline real triple_product(const vect &a, const vect &b, const vect &c) {
	return a.dot(cross(b, c));
}
#endif


inline real W(real r, real h) {
	auto norm = W_norm(h);
	return bspline(2 * r / h) / norm;
}

inline real bspline(real r) {
	if (r < 1) {
		const auto r2 = r * r;
		return 1.0 - 1.5 * r2 + 0.75 * r2 * r;
	} else if (r <= 2) {
		return 0.25 * std::pow(2 - r, 3);
	} else {
		return 0.0;
	}
}

inline real W_norm(real h) {
	if constexpr (NDIM == 1) {
		return 3.0 * h / 4.0;
	} else if constexpr (NDIM == 2) {
		return 2 * M_PI * 7 * h * h / 80.0;
	} else {
		return 2 * 2.0 * M_PI * h * h * h / 32.0;
	}
}

inline real dW_norm_dh(real h) {
	if constexpr (NDIM == 1) {
		return 3.0 / 4.0;
	} else if constexpr (NDIM == 2) {
		return 2 * M_PI * 15 * h / 80.0;
	} else {
		return 2 * 6.0 * M_PI * h * h / 32.0;
	}
}

#endif /* SRC_MATH_HPP_ */
