/*
 * main.cpp
 *
 *  Created on: Nov 28, 2019
 *      Author: dmarce1
 */

#include <hpx/hpx_init.hpp>
#include <octopart/options.hpp>
#include <octopart/tree.hpp>

int hpx_main(int argc, char *argv[]) {
	options opts;
	opts.process_options(argc,argv);
	std::vector<particle> parts = cartesian_particle_set(64);
	auto root = hpx::new_<tree>(hpx::find_here(), std::move(parts), null_range()).get();
	tree::form_tree_action()(root, root, hpx::invalid_id, std::vector<hpx::id_type>());
	tree::compute_interactions_action()(root);
	tree::write_checkpoint_action()(root, "test.chk");
	tree_stats s = hpx::async<tree::tree_statistics_action>(root).get();
	s.print();
	return hpx::finalize();

}

int main(int argc, char *argv[]) {
	hpx::init(argc, argv);
}
