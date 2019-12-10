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
	template<class Arc>
	void serialize(Arc &&arc, unsigned) {
		arc & leaf;
		arc & box;
	}
};

using gradient = general_vect<primitive_state,NDIM>;

class tree: public hpx::components::managed_component_base<tree> {

	std::vector<particle> parts;
	std::vector<gradient> grads;
	std::vector<real> Ncond;
	std::array<hpx::id_type, NCHILD> children;
	std::vector<hpx::id_type> neighbors;
	hpx::id_type parent;
	hpx::id_type self;
	range box;
	int nparts0;
	bool leaf;

public:
	tree(std::vector<particle>&&, const range&);

	void compute_next_state();
	void compute_gradients();
	void compute_interactions();
	void form_tree(const hpx::id_type&, const hpx::id_type&, std::vector<hpx::id_type>);
	void finish_tree(std::vector<hpx::id_type>);
	tree_attr get_attributes() const;
	std::vector<gradient> get_gradients(const range&, const range&) const;
	std::vector<vect> get_particle_positions(const range&) const;
	std::vector<particle> get_particles(const range&, const range&) const;
	std::array<hpx::id_type,NCHILD> get_children() const;
	tree_stats tree_statistics() const;
	void write_checkpoint(const std::string&);

	HPX_DEFINE_COMPONENT_ACTION(tree,compute_next_state);
	HPX_DEFINE_COMPONENT_ACTION(tree,compute_gradients);
	HPX_DEFINE_COMPONENT_ACTION(tree,compute_interactions);
	HPX_DEFINE_COMPONENT_ACTION(tree,form_tree);
	HPX_DEFINE_COMPONENT_ACTION(tree,get_attributes);
	HPX_DEFINE_COMPONENT_ACTION(tree,get_gradients);
	HPX_DEFINE_COMPONENT_ACTION(tree,get_particle_positions);
	HPX_DEFINE_COMPONENT_ACTION(tree,get_particles);
	HPX_DEFINE_COMPONENT_ACTION(tree,get_children);
	HPX_DEFINE_COMPONENT_ACTION(tree,tree_statistics);
	HPX_DEFINE_COMPONENT_ACTION(tree,write_checkpoint);
};

#endif /* TREE_SERVER_CPP_ */
