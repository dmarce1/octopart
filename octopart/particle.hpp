/*
 * particle.hpp
 *
 *  Created on: Dec 5, 2019
 *      Author: dmarce1
 */

#ifndef SRC_PARTICLE_HPP_
#define SRC_PARTICLE_HPP_


#include "state.hpp"



struct particle {
	conserved_state U;
	vect x;
	real V;
	real h;
	void write(FILE*) const;
	int read(FILE*);
	template<class Arc>
	void serialize(Arc&& a, unsigned) {
		a & U;
		a & h;
		a & x;
		a & V;
	}
};

#endif /* SRC_PARTICLE_HPP_ */
