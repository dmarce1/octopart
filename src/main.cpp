/*
 * main.cpp
 *
 *  Created on: Nov 28, 2019
 *      Author: dmarce1
 */

#include <hpx/hpx_init.hpp>
#include <octopart/options.hpp>
#include <octopart/tree.hpp>

int N = 256;
int hpx_main(int argc, char *argv[]) {
	real t = 0.0;
	options opts;
	opts.process_options(argc, argv);
	std::vector<particle> parts = random_particle_set(N * N);
//	std::vector<particle> parts = cartesian_particle_set(N);
	auto root = hpx::new_<tree>(hpx::find_here(), std::move(parts), null_range()).get();
	tree::set_self_and_parent_action()(root, root, hpx::invalid_id);
	tree::form_tree_action()(root, std::vector<hpx::id_type>());
	tree::compute_interactions_action()(root);
	tree::initialize_action()(root, "drift");
	tree::write_silo_action()(root, 0);
	auto dt = tree::compute_timestep_action()(root);
	dt *= 2.0 * 0.8;
	tree::compute_drift_action()(root, dt / 2.0);
	tree::finish_drift_action()(root);
	tree::set_self_and_parent_action()(root, root, hpx::invalid_id);
	tree::find_new_neighbors_action()(root);
	tree::write_silo_action()(root, 1);
	tree_stats s = hpx::async<tree::tree_statistics_action>(root).get();
	s.print();
	return hpx::finalize();

}

int main(int argc, char *argv[]) {
	hpx::init(argc, argv);
}
