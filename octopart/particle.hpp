/*
 * particle.hpp
 *
 *  Created on: Dec 5, 2019
 *      Author: dmarce1
 */

#ifndef SRC_PARTICLE_HPP_
#define SRC_PARTICLE_HPP_

#include <octopart/fixed_real.hpp>
#include <octopart/state.hpp>
#include <octopart/vect.hpp>
#include <vector>

struct primitive_particle;
struct hydro_particle;
struct timestep_particle;
struct nesting_particle;

struct particle {
	primitive_state W;
	gradient dW;
	qcon_state Q;
	real m0;
	vect x;
	vect vf;
	vect g;
	vect g0;
	real V;
	real h;
	real Nc;
	fixed_real t;
	fixed_real dt;
	fixed_real tmp;
	std::array<vect, NDIM> B;
	void write(FILE*) const;
	int read(FILE*);
	void con_to_prim();
	template<class Arc>
	void serialize(Arc &&a, unsigned) {
		a & t;
		a & dt;
		a & W;
		a & dW;
		a & Q;
		a & m0;
		a & g0;
		a & vf;
		a & g;
		a & x;
		a & V;
		a & h;
		a & B;
		a & Nc;
	}
	particle() = default;
	particle& operator=(const primitive_particle&);
	particle(const primitive_particle &p);
	particle& operator=(const timestep_particle& p);
	particle(const timestep_particle &p);
	particle& operator=(const hydro_particle&);
	particle(const hydro_particle &p);
	particle& operator=(const nesting_particle&);
	particle(const nesting_particle &p);
};

struct nesting_particle {
	fixed_real dt;
	vect x;
	real h;
	template<class Arc>
	void serialize(Arc &&a, unsigned) {
		a & dt;
		a & x;
		a & h;
	}
	nesting_particle() = default;
	nesting_particle& operator=(const particle &p);
	nesting_particle(const particle &p);
};

struct timestep_particle {
	primitive_state W;
	vect vf;
	vect x;
	real h;
	template<class Arc>
	void serialize(Arc &&a, unsigned) {
		a & W;
		a & vf;
		a & h;
		a & x;
	}
	timestep_particle() = default;
	timestep_particle& operator=(const particle &p);
	timestep_particle(const particle &p);
};

struct hydro_particle {
	primitive_state W;
	gradient dW;
	vect g;
	vect x;
	vect vf;
	real V;
	real h;
	fixed_real dt;
	fixed_real t;
	std::array<vect, NDIM> B;
	template<class Arc>
	void serialize(Arc &&a, unsigned) {
		a & B;
		a & W;
		a & dW;
		a & g;
		a & x;
		a & vf;
		a & V;
		a & h;
		a & dt;
		a & t;
	}
	hydro_particle() = default;
	hydro_particle& operator=(const particle &p);
	hydro_particle(const particle &p);
};

struct primitive_particle {
	primitive_state W;
	std::array<vect, NDIM> B;
	vect x;
	real V;
	real h;
	real Nc;
	fixed_real t;
	template<class Arc>
	void serialize(Arc &&a, unsigned) {
		a & W;
		a & x;
		a & V;
		a & h;
		a & B;
		a & t;
		a & Nc;
	}
	primitive_particle() = default;
	primitive_particle& operator=(const particle &p);
	primitive_particle(const particle &p);
};

struct gravity_part {
	real m;
	real h;
	vect x;
	template<class Arc>
	void serialize(Arc &&arc, unsigned) {
		arc & h;
		arc & m;
		arc & x;
	}
};

std::vector<particle> cartesian_particle_set(int);
std::vector<particle> random_particle_set(int);
std::vector<particle> disc_particle_set(int N);

#endif /* SRC_PARTICLE_HPP_ */
