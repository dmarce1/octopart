#include <octopart/math.hpp>
#include <octopart/options.hpp>
#include <fstream>
#include <iostream>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/async.hpp>
#include <hpx/program_options.hpp>


options options::global;

HPX_PLAIN_ACTION(options::set, set_options_action);

options& options::get() {
	return global;
}

void options::set(options o) {
	global = o;
	if( global.fpe) {
		enable_floating_point_exceptions();
	}
}

bool options::process_options(int argc, char *argv[]) {
	namespace po = hpx::program_options;

	po::options_description command_opts("options");

	command_opts.add_options() //
	("help", "produce help message")          //
	("config_file", po::value<std::string>(&config_file)->default_value(""), "configuration file") //
	("dust_only", po::value<bool>(&dust_only)->default_value(false), "treat particles as dust")           //
	("fgamma", po::value<double>(&fgamma)->default_value(7.0/5.0), "gamma for fluid gamma law")           //
	("first_order_space", po::value<bool>(&first_order_space)->default_value(false), "use 1st order spatial scheme")           //
	("first_order_time", po::value<bool>(&first_order_time)->default_value(false), "use 1st order time integration")           //
	("fpe", po::value<bool>(&fpe)->default_value(true), "enable floating point exceptions")           //
	("parts_per_node", po::value<int>(&parts_per_node)->default_value(1000), "maximum number of particles on a node")           //
	("periodic", po::value<bool>(&periodic)->default_value(false), "enable periodic boundary conditions")           //
	("problem_size", po::value<int>(&problem_size)->default_value(100), "problem size")           //
	("problem", po::value<std::string>(&problem)->default_value("sod"), "problem name")           //
	("reflecting", po::value<bool>(&reflecting)->default_value(false), "enable reflecting boundary conditions")           //
	("tmax", po::value<double>(&tmax)->default_value(1.0), "time to end simulation")           //
			;

	hpx::program_options::variables_map vm;
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

	if( first_order_space && !first_order_time) {
		printf( "ERROR: first_order_time must be enabled for first_order_space\n");
		abort();
	}
	if( reflecting && periodic ) {
		printf( "Cannot enable both reflecting and periodic boundary conditions\n");
		abort();
	}

	const auto loc = hpx::find_all_localities();
	const auto sz = loc.size();
	std::vector<hpx::future<void>> futs(sz);
	for( int i = 0; i < sz; i++) {
		futs[i] = hpx::async<set_options_action>(loc[i], *this);
	}
	hpx::wait_all(futs);
#define SHOW( opt ) std::cout << std::string( #opt ) << " = " << std::to_string(opt) << '\n';
	SHOW(dust_only);
	SHOW(fgamma);
	SHOW(first_order_space);
	SHOW(first_order_time);
	SHOW(fpe);
	SHOW(parts_per_node);
	SHOW(periodic);
	SHOW(problem_size);
	SHOW(tmax);
}
