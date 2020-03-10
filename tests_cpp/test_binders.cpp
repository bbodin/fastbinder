/*
 * test_binder.cpp
 *
 *  Created on: Jun 14, 2019
 */

#include <chrono>
#include <iostream>
#include <boost/program_options.hpp>
#include <verbose.h>
#include <binders.h>
#include "binder_data.h"

namespace po = boost::program_options;

long count_occurences (std::string str, std::string pattern) {
	size_t index = 0;
	for (std::size_t found  = str.find(pattern, 0) ; found != std::string::npos ;  found  = str.find(pattern, found + 1) ) {
        std::cout << index << " " << found << std::endl;
        index++;
	}
	return index;
}


arma::uword time_binder (std::function<arma::uword(binder_mat_t,arma::uword)> binder, binder_mat_t data, arma::uword expected ,  double weight) {

	VERBOSE_INFO("The data is " << data.n_rows << " rows " << data.n_cols << " cols");

	auto start = std::chrono::steady_clock::now();

	arma::uword res = binder(data, weight);

	auto end = std::chrono::steady_clock::now();

	COUT_STREAM << "Elapsed time in milliseconds : "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
			<< " ms" << std::endl;

	VERBOSE_ASSERT(res == expected, "Unexpected results " << res);

	return res;
}

int main (int ac, char** av) {

	po::options_description desc("Allowed options");
	desc.add_options()
	    ("help", "Help message")
	    ("verbose", po::value<int>(&VERBOSE_LEVEL)->default_value(LOG_LEVEL), "Verbosity level")
	    ("weight", po::value<double>()->default_value(0.5), "Binder const weight")
	    ("input-file", po::value< std::vector<std::string> >(), "input file")
	;



	po::variables_map vm;
	po::store(po::parse_command_line(ac, av, desc), vm);
	po::notify(vm);
	VERBOSE_LOG("VERBOSE_LEVEL=" << VERBOSE_LEVEL);
	VERBOSE_INFO("INFO_LEVEL ACTIVATED");
	VERBOSE_DEBUG("DEBUG_LEVEL ACTIVATED");

	if (vm.count("help")) {
		COUT_STREAM << desc << "\n";
	    return 0;
	}

	double weight =  vm["weight"].as<double>();

	// produce a big matrix from small one
	//auto binder_data_large = binder_data;
	//for (auto i = 0 ; i < 10 ; i ++) binder_data_large = arma::join_rows (binder_data_large,binder_data);

	//auto binder_data_andrea = binder_data;
	//for (auto i = 0 ; i < 12 ; i ++) binder_data_andrea = arma::join_cols (binder_data_andrea,binder_data);

	//time_binder ( parallel_binder , binder_data_andrea, 117, weight);
	//time_binder ( tiled_binder , binder_data_andrea, 117, weight);
	//time_binder ( opt_binder , binder_data_andrea, 117, weight);
	//time_binder ( naive_binder , binder_data_andrea, 117, weight);

	//time_binder ( opt2_binder , binder_data_large, 117, weight);
	//time_binder ( opt1_binder , binder_data_large, 117, weight);



	// produce a big matrix from small one
	//VERBOSE_DEBUG("BUILD THE DATA");
	//auto binder_data_very_large = binder_data;

	//for (auto i = 0 ; i < (500000 / binder_data.n_cols) ; i ++) binder_data_very_large = arma::join_rows (binder_data_very_large,binder_data);

	//VERBOSE_DEBUG("BUILD THE DATA ROW");
	//for (auto i = 0 ; i <  (10000 / binder_data.n_rows) ; i ++) binder_data_very_large = arma::join_cols (binder_data_very_large,binder_data);

	//VERBOSE_DEBUG("RUN ALGORITHM");

	binder_mat_t  binder_data_very_large (500000, 10000, arma::fill::ones);

	time_binder ( parallel_binder , binder_data_very_large, 117, weight);

}
