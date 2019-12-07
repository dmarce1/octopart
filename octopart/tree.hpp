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

static constexpr int NPART_MAX = 100;
#if(NDIM==1)
static constexpr int NNGB = 4;
static constexpr real CV = 1.0;
#elif(NDIM==2)
static constexpr int NNGB = 16;
static constexpr real CV = M_PI;
#else
static constexpr int NNGB = 32;
static constexpr real CV = 4.0 * M_PI / 3.0;
#endif

class tree: public hpx::components::managed_component_base<tree> {

	std::vector<particle> parts;
	std::vector<std::shared_ptr<particle>> nparts;
	std::array<hpx::id_type, NCHILD> children;
	hpx::id_type parent;
	hpx::id_type self;
	range box;
	bool leaf;

public:
	tree(std::vector<particle>&&, const range&);
	real first_sweep();
	std::vector<std::shared_ptr<particle>> particle_search(const range&, const hpx::id_type&, const hpx::id_type&);
	std::vector<std::shared_ptr<particle>> particle_gather(const range&, const hpx::id_type&);
	void write_checkpoint(const std::string&);
	void form_tree(const hpx::id_type&, const hpx::id_type&);
	void set_self(const hpx::id_type&);

	HPX_DEFINE_COMPONENT_ACTION(tree,write_checkpoint,write_checkpoint_action);
	HPX_DEFINE_COMPONENT_ACTION(tree,first_sweep,first_sweep_action);
	HPX_DEFINE_COMPONENT_ACTION(tree,particle_search,particle_search_action);
	HPX_DEFINE_COMPONENT_ACTION(tree,particle_gather,particle_gather_action);
	HPX_DEFINE_COMPONENT_ACTION(tree,form_tree,form_tree_action);
	HPX_DEFINE_COMPONENT_ACTION(tree,set_self,set_self_action);

};

HPX_REGISTER_ACTION_DECLARATION(tree::write_checkpoint_action);
HPX_REGISTER_ACTION_DECLARATION(tree::first_sweep_action);
HPX_REGISTER_ACTION_DECLARATION(tree::particle_search_action);
HPX_REGISTER_ACTION_DECLARATION(tree::particle_gather_action);
HPX_REGISTER_ACTION_DECLARATION(tree::form_tree_action);
HPX_REGISTER_ACTION_DECLARATION(tree::set_self_action);

#endif /* TREE_SERVER_CPP_ */
