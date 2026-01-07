#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include <algorithm>
#include "utils.hpp"
#include "qoi_utils.hpp"
#include "Reconstructor/Reconstructor.hpp"
#include "Synthesizer4GE.hpp"

using namespace MDR;

const std::string s3d_data_file_prefix = "/Users/xuanwu/github/datasets/S3D/data/";
const std::string s3d_rdata_file_prefix = "/Users/xuanwu/github/datasets/S3D/refactor/";
const std::vector<std::string> species = {"H2", "O2", "H2O", "H", "O", "OH", "HO2", "H2O2"};

std::vector<double> Xi_ori;
std::vector<double> Xj_ori;
double * Xi_dec = NULL;
double * Xj_dec = NULL;
double * XiXj_ori = NULL;
std::vector<double> error_XiXj;
std::vector<double> error_est_XiXj;


template<class T>
bool halfing_error_XiXj_uniform(const T * Xi, const T * Xj, size_t n, const double tau, std::vector<double>& ebs){
	double eb_Xi = ebs[0];
	double eb_Xj = ebs[1];
	double max_value = 0;
	int max_index = 0;
	int n_variable = ebs.size();
	// double max_err = -1;
	// int max_err_index = 0;
	for(int i=0; i<n; i++){
		double e_XiXj = compute_bound_multiplication(Xi[i], Xj[i], eb_Xi, eb_Xj);
        double XiXj = Xi[i] * Xj[i];

		error_est_XiXj[i] = e_XiXj;
		error_XiXj[i] = XiXj - XiXj_ori[i];
		if(max_value < error_est_XiXj[i]){
			max_value = error_est_XiXj[i];
			max_index = i;
		}
		// if(max_err < fabs(error_XiXj[i])){
		// 	max_err = fabs(error_XiXj[i]);
		// 	max_err_index = i;
		// }
	}
	// std::cout << "XiXj: max estimated error = " << max_value << ", index = " << max_index << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		auto i = max_index;
		double estimate_error = max_value;
		double eb_Xi = ebs[0];
		double eb_Xj = ebs[1];
		while(estimate_error > tau){
    		// std::cout << "uniform decrease\n";
			eb_Xi = eb_Xi / 1.5;
			eb_Xj = eb_Xj / 1.5;
			estimate_error = compute_bound_multiplication(Xi[i], Xj[i], eb_Xi, eb_Xj);
		}
		ebs[0] = eb_Xi;
		ebs[1] = eb_Xj;
		return false;
	}
	return true;
}


int main(int argc, char **argv){

	using T = double;
	int argv_id = 1;
	int id_i = atoi(argv[argv_id++]);
	int id_j = atoi(argv[argv_id++]);
	double target_rel_eb = atof(argv[argv_id++]);
	std::string data_prefix_path = argv[argv_id++];
	std::string data_file_prefix = data_prefix_path + "/data/";
	std::string rdata_file_prefix = data_prefix_path + "/refactor/";

	std::vector<int> index = {id_i, id_j};

    size_t num_elements = 0;
    Xi_ori = MGARD::readfile<T>((data_file_prefix + species[id_i] + ".dat").c_str(), num_elements);
    Xj_ori = MGARD::readfile<T>((data_file_prefix + species[id_j] + ".dat").c_str(), num_elements);
    std::vector<double> ebs;
    ebs.push_back(compute_value_range(Xi_ori)*target_rel_eb);
    ebs.push_back(compute_value_range(Xj_ori)*target_rel_eb);
	int n_variable = ebs.size();

	std::vector<double> var_range = {compute_value_range(Xi_ori), compute_value_range(Xj_ori)};

    std::vector<T> XiXj(num_elements);
	for(int i=0; i<num_elements; i++){
		XiXj[i] = Xi_ori[i] * Xj_ori[i];
	}
	XiXj_ori = XiXj.data();
	double tau = compute_value_range(XiXj)*target_rel_eb;
  
    std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<double>(num_elements));
	std::vector<size_t> total_retrieved_sizes(n_variable, 0);

	struct timespec start, end;
	int err;
	double elapsed_time;

	err = clock_gettime(CLOCK_REALTIME, &start);

    int iter = 0;
    int max_iter = 5;
	bool tolerance_met = false;
	double max_act_error = 0, max_est_error = 0;
    std::vector<int> current_ind(n_variable, -1);
    T * reconstructed_data = (T *) malloc(num_elements * sizeof(T));
    while((!tolerance_met) && (iter < max_iter)){
    	iter ++;
	    for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + species[index[i]];
            double file_eb = 0.1;
            auto file_ind = find_index(ebs[i]/var_range[i], file_eb);
            // std::cout << "file_ind = " << file_ind << std::endl;
            // std::cout << "Requested relative tolerance = " << ebs[i]/var_range[i] << ", expected relative tolerance = " << file_eb << "\n"; 
            if(file_ind > current_ind[i]){
                for(int j=current_ind[i]+1; j<=file_ind; j++){
                    std::string filename = rdir_prefix + "_refactored/SZ3_delta_eb_" + std::to_string(j) + ".bin";
                    size_t n = 0;
                    auto cmpData = MGARD::readfile<char>(filename.c_str(), n);
                    total_retrieved_sizes[i] += n;
                    SZ3_decompress(cmpData.data(), n, reconstructed_data);
					for(int j=0; j<num_elements; j++){
						reconstructed_vars[i][j] += reconstructed_data[j];
					}
                }
                current_ind[i] = file_ind;            
            }
	    }
	    Xi_dec = reconstructed_vars[0].data();
	    Xj_dec = reconstructed_vars[1].data();
	    // MGARD::print_statistics(Xi_ori.data(), Xi_dec, num_elements);
	    // MGARD::print_statistics(Xj_ori.data(), Xj_dec, num_elements);
	    error_XiXj = std::vector<double>(num_elements);
	    error_est_XiXj = std::vector<double>(num_elements);
		// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
	    // MDR::print_vec(ebs);
	    tolerance_met = halfing_error_XiXj_uniform(Xi_dec, Xj_dec, num_elements, tau, ebs);
		// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
	    // MDR::print_vec(ebs);
	    // std::cout << "requested error = " << tau << std::endl;
	    // print_max_abs("error_est", error_est_XiXj);   	
	    // print_max_abs("actual error", error_XiXj);
		std::string name = "X" + std::to_string(id_i) + "X" + std::to_string(id_j);
	    // std::cout << name << " requested error = " << tau << std::endl;
	    max_est_error = print_max_abs(name + " error_est", error_est_XiXj);   	
	    max_act_error = print_max_abs(name + " actual error", error_XiXj);
    }
	err = clock_gettime(CLOCK_REALTIME, &end);
	elapsed_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000;

	std::cout << "requested_error = " << tau << std::endl;
	std::cout << "max_est_error = " << max_est_error << std::endl;
	std::cout << "max_act_error = " << max_act_error << std::endl;
	std::cout << "iter = " << iter << std::endl;
   
   	size_t total_size = std::accumulate(total_retrieved_sizes.begin(), total_retrieved_sizes.end(), 0);
	double cr = n_variable * num_elements * sizeof(T) * 1.0 / total_size;
	std::cout << "each retrieved size:";
    for(int i=0; i<n_variable; i++){
        std::cout << total_retrieved_sizes[i] << ", ";
    }
    std::cout << std::endl;
	// MDR::print_vec(total_retrieved_size);
	std::cout << "aggregated cr = " << cr << std::endl;
	std::cout << "bitrate = " << ((sizeof(T) * 8) / cr) << std::endl;
	printf("elapsed_time = %.6f\n", elapsed_time);

	return 0;
}