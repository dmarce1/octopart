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
#include <hpx/runtime/components/server/migrate_component.hpp>
#include <octopart/tree_stats.hpp>

struct tree_attr {
	bool leaf;
	bool dead;
	int nparts;
	range box;
	template<class Arc>
	void serialize(Arc &&arc, unsigned) {
		arc & leaf;
		arc & dead;
		arc & nparts;
		arc & box;
	}
};

struct sibling_attr {
	hpx::id_type id;
	range box;
	vect pshift;
	template<class Arc>
	void serialize(Arc &&arc, unsigned) {
		arc & id;
		arc & box;
		arc & pshift;
	}
};

struct mass_attr {
	vect com;
	real mtot;
	real rmaxs;
	real rmaxb;
	bool leaf;
	template<class Arc>
	void serialize(Arc &&arc, unsigned) {
		arc & com;
		arc & mtot;
		arc & rmaxs;
		arc & rmaxb;
		arc & leaf;
	}
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

class tree: public hpx::components::migration_support<hpx::components::component_base<tree>> {
	std::vector<particle> new_parts;
	std::vector<particle> parts;
	std::vector<gradient> grad_lim;
	std::vector<gradient> grad;
	std::vector<conserved_state> dudt;
	std::vector<vect> mass_flux;
	std::vector<real> Ncond;
	std::array<hpx::id_type, NCHILD> children;
	std::array<int, NCHILD> child_loads;
	std::vector<sibling_attr> siblings;
	hpx::id_type parent;
	hpx::id_type self;
	range root_box;
	range box;
	int nparts0;
	bool leaf;
	bool dead;
	mass_attr mass;
	std::shared_ptr<hpx::lcos::local::mutex> mtx;

public:
	tree();
	tree(std::vector<particle>&&, const range&, const range&);

	mass_attr compute_mass_attributes();
	void compute_drift(real);
	void compute_gradients();
	void compute_gravity(std::vector<hpx::id_type>, std::vector<mass_attr>);
	void compute_next_state(real);
	void compute_time_derivatives(real);
	real compute_timestep() const;
	void compute_interactions();
	int compute_workload();
	void create_children();
	std::vector<particle> destroy();
	void form_tree(std::vector<hpx::id_type>, bool = true);
	tree_attr finish_drift();
	tree_attr get_attributes() const;
	std::array<hpx::id_type, NCHILD> get_children() const;
	std::vector<gradient> get_gradients(range, range, const vect&) const;
	std::vector<gravity_part> get_gravity_particles() const;
	mass_attr get_mass_attributes() const;
	hpx::id_type get_parent() const;
	std::vector<vect> get_particle_positions(range, const vect&) const;
	std::vector<particle> get_particles(range, range, const vect&) const;
	void initialize(const std::string&);
	void redistribute_workload(int, int);
	void send_particles(const std::vector<particle>&, const vect&);
	void set_central_force();
	void set_self_and_parent(const hpx::id_type, const hpx::id_type);
	tree_stats tree_statistics() const;
	void write_checkpoint(const std::string&) const;
	void write_silo(int) const;

	template<class Arc>
	void serialize(Arc &&arc, unsigned) {
		arc & parts;
		arc & children;
		arc & child_loads;
		arc & siblings;
		arc & parent;
		arc & self;
		arc & root_box;
		arc & box;
		arc & leaf;
		arc & mass;
	}

	HPX_DEFINE_COMPONENT_ACTION(tree,compute_mass_attributes);
	HPX_DEFINE_COMPONENT_ACTION(tree,compute_drift);
	HPX_DEFINE_COMPONENT_ACTION(tree,compute_gradients);
	HPX_DEFINE_COMPONENT_ACTION(tree,compute_next_state);
	HPX_DEFINE_COMPONENT_ACTION(tree,compute_time_derivatives);
	HPX_DEFINE_COMPONENT_ACTION(tree,compute_timestep);
	HPX_DEFINE_COMPONENT_ACTION(tree,compute_interactions);
	HPX_DEFINE_COMPONENT_ACTION(tree,compute_workload);
	HPX_DEFINE_COMPONENT_ACTION(tree,destroy);
	HPX_DEFINE_COMPONENT_ACTION(tree,form_tree);
	HPX_DEFINE_COMPONENT_ACTION(tree,initialize);
	HPX_DEFINE_COMPONENT_ACTION(tree,finish_drift);
	HPX_DEFINE_COMPONENT_ACTION(tree,compute_gravity);
	HPX_DEFINE_COMPONENT_ACTION(tree,redistribute_workload);
	HPX_DEFINE_COMPONENT_ACTION(tree,set_central_force);
	HPX_DEFINE_COMPONENT_ACTION(tree,tree_statistics);
	HPX_DEFINE_COMPONENT_ACTION(tree,write_checkpoint);
	HPX_DEFINE_COMPONENT_ACTION(tree,write_silo);
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_attributes);
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_gradients);
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_gravity_particles);
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_mass_attributes);
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_parent);
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_particle_positions);
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_particles);
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_children);
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,send_particles);
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,set_self_and_parent);

};

#endif /* TREE_SERVER_CPP_ */
