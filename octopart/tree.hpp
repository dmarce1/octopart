/*
 * tree.hpp
 *
 *  Created on: Dec 5, 2019
 *      Author: dmarce1
 */

#ifndef TREE_SERVER_CPP_
#define TREE_SERVER_CPP_

#include "particle.hpp"
#include "range.hpp"

#include <hpx/include/components.hpp>
#include <octopart/tree_stats.hpp>

struct tree_attr {
	bool leaf;
	range box;
	int nparts;
	template<class Arc>
	void serialize(Arc &&arc, unsigned) {
		arc & leaf;
		arc & box;
		arc & nparts;
	}
};

struct neighbor_attr {
	hpx::id_type id;
	range box;
};

class tree: public hpx::components::managed_component_base<tree> {

	std::vector<particle> new_parts;
	std::vector<particle> parts;
	std::vector<gradient> grad_lim;
	std::vector<gradient> grad;
	std::vector<conserved_state> dudt;
	std::vector<real> Ncond;
	std::array<hpx::id_type, NCHILD> children;
	std::vector<neighbor_attr> neighbors;
	hpx::id_type parent;
	hpx::id_type self;
	range box;
	int nparts0;
	bool leaf;
	hpx::lcos::local::mutex mtx;

public:
	tree(std::vector<particle>&&, const range&);

	void compute_drift(real);
	void compute_gradients();
	void compute_next_state(real, real);
	void compute_time_derivatives(real);
	real compute_timestep() const;
	void compute_interactions();
	void create_children();
	std::vector<particle> destroy();
	void form_tree(const hpx::id_type&, const hpx::id_type&, std::vector<hpx::id_type>);
	tree_attr finish_drift();
	void finish_tree(std::vector<hpx::id_type>);
	tree_attr get_attributes() const;
	std::vector<gradient> get_gradients(const range&, const range&) const;
	std::vector<vect> get_particle_positions(const range&) const;
	std::vector<particle> get_particles(const range&, const range&) const;
	std::array<hpx::id_type, NCHILD> get_children() const;
	void send_particles(const std::vector<particle>&);
	tree_stats tree_statistics() const;
	void write_checkpoint(const std::string&);

	HPX_DEFINE_COMPONENT_ACTION(tree,compute_drift);
	HPX_DEFINE_COMPONENT_ACTION(tree,compute_gradients);
	HPX_DEFINE_COMPONENT_ACTION(tree,compute_next_state);
	HPX_DEFINE_COMPONENT_ACTION(tree,compute_time_derivatives);
	HPX_DEFINE_COMPONENT_ACTION(tree,compute_timestep);
	HPX_DEFINE_COMPONENT_ACTION(tree,compute_interactions);
	HPX_DEFINE_COMPONENT_ACTION(tree,destroy);
	HPX_DEFINE_COMPONENT_ACTION(tree,form_tree);
	HPX_DEFINE_COMPONENT_ACTION(tree,get_attributes);
	HPX_DEFINE_COMPONENT_ACTION(tree,get_gradients);
	HPX_DEFINE_COMPONENT_ACTION(tree,get_particle_positions);
	HPX_DEFINE_COMPONENT_ACTION(tree,get_particles);
	HPX_DEFINE_COMPONENT_ACTION(tree,get_children);
	HPX_DEFINE_COMPONENT_ACTION(tree,finish_drift);
	HPX_DEFINE_COMPONENT_ACTION(tree,send_particles);
	HPX_DEFINE_COMPONENT_ACTION(tree,tree_statistics);
	HPX_DEFINE_COMPONENT_ACTION(tree,write_checkpoint);
};

#endif /* TREE_SERVER_CPP_ */
