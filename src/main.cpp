#include <hpx/hpx_init.hpp>
#include <octopart/options.hpp>
#include <octopart/profiler.hpp>
#include <octopart/tree.hpp>

int N = 32;
int hpx_main(int argc, char *argv[]) {
	real t = 0.0;
	options opts;
	opts.process_options(argc, argv);
//	std::vector<particle> parts = random_particle_set(N * N);
	std::vector<particle> parts;
	if (opts.problem == "kepler") {
		parts = disc_particle_set(opts.problem_size);
	} else {
		parts = cartesian_particle_set(opts.problem_size);
	}
	range box;
	for (int dim = 0; dim < NDIM; dim++) {
		box.min[dim] = -0.5;
		box.max[dim] = +0.5;
	}
	auto root = hpx::new_<tree>(hpx::find_here(), std::move(parts), box, null_range()).get();
	tree::set_self_and_parent_action()(root, root, hpx::invalid_id);
	tree::form_tree_action()(root, std::vector<hpx::id_type>(1, root), true);
	tree::compute_interactions_action()(root);
	tree::initialize_action()(root, opts.problem);
	if (opts.gravity) {
		tree::compute_mass_attributes_action()(root);
		tree::compute_gravity_action()(root, std::vector<hpx::id_type>(1, root), std::vector<mass_attr>());
	}
	tree::write_silo_action()(root, 0);
	for (int i = 0; t < opts.tmax; i++) {
		auto dt = tree::compute_timestep_action()(root);
		dt *= 0.4;
		tree_stats s = tree::tree_statistics_action()(root);
		printf("Step = %i t = %e  dt = %e Nparts = %i Nleaves = %i Max Level = %i Mass = %e Momentum = ", i, t.get(), dt.get(), s.nparts, s.nleaves,
				s.max_level, s.mass.get());
		for (int dim = 0; dim < NDIM; dim++) {
			printf("%e ", s.momentum[dim].get());
		}
		printf("Energy = %e\n", s.energy.get());
		if (!opts.dust_only) {
			tree::compute_gradients_action()(root);
			tree::compute_time_derivatives_action()(root, dt);
			tree::compute_next_state_action()(root, dt);
		}
		tree::compute_drift_action()(root, dt);
		tree::finish_drift_action()(root);
		tree::set_self_and_parent_action()(root, root, hpx::invalid_id);
		tree::form_tree_action()(root, std::vector<hpx::id_type>(1, root), true);
		tree::compute_interactions_action()(root);
		if (opts.gravity) {
			tree::compute_mass_attributes_action()(root);
			tree::compute_gravity_action()(root, std::vector<hpx::id_type>(1, root), std::vector<mass_attr>());
		}
		tree::write_silo_action()(root, i + 1);
		t += dt;
	}
	FILE *fp = fopen("profile.txt", "wt");
	profiler_output(fp);
	fclose(fp);
	return hpx::finalize();

}

int main(int argc, char *argv[]) {
	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };

	hpx::init(argc, argv, cfg);
}
