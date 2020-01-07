#include <hpx/hpx_init.hpp>
#include <octopart/options.hpp>
#include <octopart/profiler.hpp>
#include <octopart/tree.hpp>

hpx::id_type root;

void solve_gravity(fixed_real t, fixed_real dt, real beta) {
	static const auto opts = options::get();
	if (opts.problem == "kepler" || opts.problem == "rt") {
		tree::set_problem_force_action()(root);
	} else if (opts.gravity) {
		tree::compute_mass_attributes_action()(root);
		tree::compute_gravity_action()(root, std::vector<hpx::id_type>(1, root), std::vector<mass_attr>());
	}
	if (opts.problem == "kepler" || opts.problem == "rt" || opts.gravity) {
		tree::apply_gravity_action()(root, t, dt, beta);
	}
}

void drift(fixed_real t, fixed_real dt) {
	tree::compute_drift_action()(root, t, dt);
	tree::finish_drift_action()(root);
	tree::set_self_and_parent_action()(root, root, hpx::invalid_id);
	tree::form_tree_action()(root, std::vector<hpx::id_type>(1, root), true);
	tree::compute_interactions_action()(root, t);
}

void hydro(fixed_real t, fixed_real dt, real beta) {
	static const auto opts = options::get();
	if (!opts.dust_only) {
		tree::get_neighbor_particles_action()(root);
		tree::compute_gradients_action()(root);
		tree::compute_time_derivatives_action()(root, t, dt, beta);
	}
}

void init(bool t0) {
	static const auto opts = options::get();
	tree::set_self_and_parent_action()(root, root, hpx::invalid_id);
	tree::form_tree_action()(root, std::vector<hpx::id_type>(1, root), true);
	tree::compute_interactions_action()(root, fixed_real(0.0));
	if (t0) {
		tree::initialize_action()(root, opts.problem);
	}
}

void write_checkpoint(int i) {
	tree::write_silo_action()(root, i + 1);
}

fixed_real timestep(fixed_real t) {
	static const auto opts = options::get();
	tree::get_neighbor_particles_action()(root);
	fixed_real dt = tree::compute_timestep_action()(root, t);
	return dt;
}

auto statistics() {
	return tree::tree_statistics_action()(root);
}

int N = 32;
int hpx_main(int argc, char *argv[]) {
	fixed_real t = 0.0;
	options opts;
	opts.process_options(argc, argv);
	std::vector<particle> parts;
	bool t0;
	if (opts.checkpoint != "") {
		printf("Reading %s\n", opts.checkpoint.c_str());
		FILE *fp = fopen(opts.checkpoint.c_str(), "rb");
		if (fp == NULL) {
			printf("Could not find %s\n", opts.checkpoint.c_str());
			abort();
		} else {
			particle p;
			real dummy;
			int cnt = fread(&dummy, sizeof(real), 1, fp);
			if (cnt == 0) {
				printf("Empty checkpoint\n");
				return hpx::finalize();
			}
			while (p.read(fp)) {
				parts.push_back(p);
			}
			fclose(fp);
		}
		t0 = false;
	} else {
//	std::vector<particle> parts = random_particle_set(N * N);
		if (opts.problem == "kepler") {
			parts = disc_particle_set(opts.problem_size);
		} else {
			parts = cartesian_particle_set(opts.problem_size);
		}
		t0 = true;
		for (auto &p : parts) {
			p.x = p.x * opts.grid_size;
		}
	}
	range box;
	for (int dim = 0; dim < NDIM; dim++) {
		box.min[dim] = -opts.grid_size / 2.0;
		box.max[dim] = +opts.grid_size / 2.0;
	}
	root = hpx::new_<tree>(hpx::find_here(), std::move(parts), box, null_range()).get();
	init(t0);
	tree::set_drift_velocity_action()(root);
	solve_gravity(0.0, 0.0, 0.0);
	write_checkpoint(0);
	int oi = 0;
	int i = 0;
	fixed_real last_output = 0.0;
	while (t < fixed_real(opts.tmax)) {
//		printf( "%e\n", double(t.next_bin()));
		fixed_real dt = timestep(t);
		auto s = statistics();
		printf("Step = %i t = %e  dt = %e Nparts = %i Nleaves = %i Max Level = %i Mass = %e Momentum = ", i, double(t), double(dt), s.nparts, s.nleaves,
				s.max_level, s.mass.get());
		for (int dim = 0; dim < NDIM; dim++) {
			printf("%e ", s.momentum[dim].get());
		}
		printf("Energy = %e\n", s.energy.get());
		solve_gravity(t, dt, 0.5);
		tree::set_drift_velocity_action()(root);
		hydro(t, dt, 0.5);
		drift(t, dt);
		hydro(t, dt, 0.5);
		solve_gravity(t, dt, 0.5);
		tree::compute_next_state_action()(root, t, dt);
		t += dt;
		i++;
		if (int((last_output / fixed_real(opts.output_freq))) != int(((t / fixed_real(opts.output_freq))))) {
			last_output = t;
			write_checkpoint(++oi);
			printf("output %i\n", oi);
		}
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
