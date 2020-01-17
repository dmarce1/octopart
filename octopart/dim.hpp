/*
 * dim.hpp
 *
 *  Created on: Nov 29, 2019
 *      Author: dmarce1
 */

#ifndef DIM_HPP_
#define DIM_HPP_

#include <cmath>

#define NDIM 1
static constexpr int NCHILD = 1 << NDIM;
static constexpr int NSIBLING = 2 * NDIM;
static constexpr int NNEIGHBOR = std::pow(3,NDIM);

#endif /* DIM_HPP_ */
