#include <octopart/math.hpp>
#include <octopart/options.hpp>
#include <fstream>
#include <iostream>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/async.hpp>
#include <boost/program_options.hpp>
#include <hpx/runtime/threads/run_as_os_thread.hpp>
#include <thread>

options options::global;

HPX_PLAIN_ACTION(options::set, set_options_action);

options& options::get() {
	return global;
}

void options::set(options o) {
	global = o;
	if (global.fpe) {
		enable_floating_point_exceptions();
	}
}

bool options::process_options(int argc, char *argv[]) {
	std::thread([&]() {
		namespace po = boost::program_options;

		po::options_description command_opts("options");

		command_opts.add_options() //
		("help", "produce help message") //
		("cfl", po::value<double>(&cfl)->default_value(0.4), "CFL factor") //
		("checkpoint", po::value<std::string>(&checkpoint)->default_value(""), "checkpoint file") //
		("config_file", po::value<std::string>(&config_file)->default_value(""), "configuration file") //
		("dust_only", po::value<bool>(&dust_only)->default_value(false), "treat particles as dust") //
		("eulerian", po::value<bool>(&eulerian)->default_value(false), "treat particles as dust") //
		("kep_eps", po::value<double>(&kep_eps)->default_value(0.05), "softening length for central force") //
		("fgamma", po::value<double>(&fgamma)->default_value(7.0 / 5.0), "gamma for fluid gamma law") //
		("first_order_space", po::value<bool>(&first_order_space)->default_value(false), "use 1st order spatial scheme") //
		("first_order_time", po::value<bool>(&first_order_time)->default_value(false), "use 1st order time integration") //
		("fpe", po::value<bool>(&fpe)->default_value(true), "enable floating point exceptions") //
		("global_time", po::value<bool>(&global_time)->default_value(false), "enable global time-stepping") //
		("gravity", po::value<bool>(&gravity)->default_value(false), "enable gravity") //
		("grid_size", po::value<double>(&grid_size)->default_value(1.0), "size of grid") //
		("output_freq", po::value<double>(&output_freq)->default_value(-1), "output frequency") //
		("parts_per_node", po::value<int>(&parts_per_node)->default_value(1000), "maximum number of particles on a node") //
		("x_periodic", po::value<bool>(&x_periodic)->default_value(false), "enable periodic boundary conditions for x") //
		("y_periodic", po::value<bool>(&y_periodic)->default_value(false), "enable periodic boundary conditions for y") //
		("z_periodic", po::value<bool>(&z_periodic)->default_value(false), "enable periodic boundary conditions for z") //
		("problem_size", po::value<int>(&problem_size)->default_value(100), "problem size") //
		("problem", po::value<std::string>(&problem)->default_value("sod"), "problem name") //
		("x_reflecting", po::value<bool>(&x_reflecting)->default_value(false), "enable reflecting boundary conditions for x") //
		("y_reflecting", po::value<bool>(&y_reflecting)->default_value(false), "enable reflecting boundary conditions for y") //
		("z_reflecting", po::value<bool>(&z_reflecting)->default_value(false), "enable reflecting boundary conditions for z") //
		("theta", po::value<double>(&theta)->default_value(0.35), "theta for Barnes-Hut") //
		("tmax", po::value<double>(&tmax)->default_value(1.0), "time to end simulation") //
				;

		boost::program_options::variables_map vm;
		po::store(po::parse_command_line(argc, argv, command_opts), vm);
		po::notify(vm);
		if (vm.count("help")) {
			std::cout << command_opts << "\n";
			return false;
		}
		if (!config_file.empty()) {
			std::ifstream cfg_fs { vm["config_file"].as<std::string>() };
			if (cfg_fs) {
				po::store(po::parse_config_file(cfg_fs, command_opts), vm);
			} else {
				printf("Configuration file %s not found!\n", config_file.c_str());
				return false;
			}
		}
		po::notify(vm);

		if (output_freq <= 0.0) {
			output_freq = tmax / 100.0;
		}

		if (first_order_space && !first_order_time) {
			printf("ERROR: first_order_time must be enabled for first_order_space\n");
			abort();
		}
		return true;
	}
	).join();
	const auto loc = hpx::find_all_localities();
	const auto sz = loc.size();
	std::vector<hpx::future<void>> futs;
	set(*this);
	for (int i = 1; i < sz; i++) {
		futs.push_back(hpx::async<set_options_action>(loc[i], *this));
	}
	hpx::wait_all(futs);
#define SHOW( opt ) std::cout << std::string( #opt ) << " = " << std::to_string(opt) << '\n';
	SHOW(cfl);
	SHOW(dust_only);
	SHOW(fgamma);
	SHOW(first_order_space);
	SHOW(first_order_time);
	SHOW(fpe);
	SHOW(global_time);
	SHOW(gravity);
	SHOW(eulerian);
	SHOW(parts_per_node);
	SHOW(x_periodic);
	SHOW(y_periodic);
	SHOW(z_periodic);
	SHOW(x_reflecting);
	SHOW(y_reflecting);
	SHOW(z_reflecting);
	SHOW(problem_size);
	SHOW(theta);
	SHOW(tmax);
	return true;
}
