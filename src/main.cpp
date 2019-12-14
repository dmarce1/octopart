#include <hpx/hpx_init.hpp>
#include <octopart/options.hpp>
#include <octopart/tree.hpp>

int N = 100;
int hpx_main(int argc, char *argv[]) {
	real t = 0.0;
	options opts;
	opts.process_options(argc, argv);
//	std::vector<particle> parts = random_particle_set(N * N);
	std::vector<particle> parts = cartesian_particle_set(N);
	auto root = hpx::new_<tree>(hpx::find_here(), std::move(parts), null_range(), null_range()).get();
	tree::set_self_and_parent_action()(root, root, hpx::invalid_id);
	tree::form_tree_action()(root, std::vector<hpx::id_type>());
	tree::compute_interactions_action()(root);
	tree::initialize_action()(root, opts.problem);
	tree::write_silo_action()(root, 0);
	for (int i = 0; i < 1000; i++) {
		auto dt = tree::compute_timestep_action()(root);
		dt *= 0.4;
		tree_stats s = tree::tree_statistics_action()(root);
		printf("Step = %i t = %e  dt = %e Nparts = %i Nleaves = %i Max Level = %i\n", i, t.get(), dt.get(), s.nparts, s.nleaves, s.max_level);
		tree::compute_drift_action()(root, dt / 2.0);
		tree::finish_drift_action()(root);
		tree::set_self_and_parent_action()(root, root, hpx::invalid_id);
		tree::form_tree_action()(root, std::vector<hpx::id_type>());
		tree::compute_interactions_action()(root);
		if (!opts.dust_only) {
			tree::compute_gradients_action()(root);
			tree::compute_time_derivatives_action()(root, dt);
			tree::compute_next_state_action()(root, dt);
		}
		tree::compute_drift_action()(root, dt);
		tree::finish_drift_action()(root);
		tree::set_self_and_parent_action()(root, root, hpx::invalid_id);
		tree::form_tree_action()(root, std::vector<hpx::id_type>());
		tree::compute_interactions_action()(root);
		tree::write_silo_action()(root, i + 1);
		t += dt;
	}
	return hpx::finalize();

}

int main(int argc, char *argv[]) {
	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };

	hpx::init(argc, argv, cfg);
}
