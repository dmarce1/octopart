/*
 * math.hpp
 *
 *  Created on: Nov 29, 2019
 *      Author: dmarce1
 */

#ifndef RANGE_HPP_
#define RANGE_HPP_

#include "vect.hpp"
#include <limits>

using range = std::pair<vect,vect>;

bool in_range(const vect&, const range&);
bool in_range(const range&, const range&);
bool ranges_intersect(const range&, const range&);
real range_volume(const range&);
range null_range();

#endif /* MATH_HPP_ */

