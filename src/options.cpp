#include <octopart/options.hpp>
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>


options options::global;

bool options::process_options(int argc, char *argv[]) {
	namespace po = boost::program_options;

	po::options_description command_opts("options");

	command_opts.add_options() //
	("help", "produce help message")          //
	("config_file", po::value<std::string>(&config_file)->default_value(""), "configuration file") //
	("dust_only", po::value<bool>(&dust_only)->default_value(false), "treat particles as dust")           //
	("first_order_space", po::value<bool>(&first_order_space)->default_value(false), "use 1st order spatial scheme")           //
	("first_order_time", po::value<bool>(&first_order_time)->default_value(false), "use 1st order time integration")           //
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
#define SHOW( opt ) std::cout << std::string( #opt ) << " = " << std::to_string(opt) << '\n';
	SHOW(dust_only);
	SHOW(first_order_space);
	SHOW(first_order_time);
}
