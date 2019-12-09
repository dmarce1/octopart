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
	range box;
	range reach;
	bool leaf;
	std::vector<hpx::id_type> children;
	template<class Arc>
	void serialize(Arc &&arc, unsigned) {
		arc & box;
		arc & reach;
		arc & leaf;
		arc & children;
	}
};

class tree: public hpx::components::managed_component_base<tree> {

	std::vector<particle> parts;
	std::array<hpx::id_type, NCHILD> children;
	std::vector<hpx::id_type> neighbors;
	hpx::id_type parent;
	hpx::id_type self;
	range box;
	range reach;
	bool leaf;

public:
	tree(std::vector<particle>&&, const range&);

	void find_neighbors(std::vector<hpx::id_type>, bool = true);
	void form_tree(const hpx::id_type&, const hpx::id_type&);
	tree_attr get_attributes() const;
	std::vector<vect> get_particle_positions(const range&) const;
	void set_self(const hpx::id_type&);
	void smoothing_length_init();
	bool smoothing_length_iter();
	tree_stats tree_statistics() const;
	void write_checkpoint(const std::string&);

	HPX_DEFINE_COMPONENT_ACTION(tree,find_neighbors);
	HPX_DEFINE_COMPONENT_ACTION(tree,form_tree);
	HPX_DEFINE_COMPONENT_ACTION(tree,get_attributes);
	HPX_DEFINE_COMPONENT_ACTION(tree,get_particle_positions);
	HPX_DEFINE_COMPONENT_ACTION(tree,set_self);
	HPX_DEFINE_COMPONENT_ACTION(tree,smoothing_length_init);
	HPX_DEFINE_COMPONENT_ACTION(tree,smoothing_length_iter);
	HPX_DEFINE_COMPONENT_ACTION(tree,tree_statistics);
	HPX_DEFINE_COMPONENT_ACTION(tree,write_checkpoint);

};

#endif /* TREE_SERVER_CPP_ */
