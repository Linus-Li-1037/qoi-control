#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include "utils.hpp"
#include "qoi_utils.hpp"
#include "Reconstructor/Reconstructor.hpp"
#include "Synthesizer4GE.hpp"

using namespace MDR;

std::vector<double> P_ori;
std::vector<double> D_ori;
std::vector<double> Vx_ori;
std::vector<double> Vy_ori;
std::vector<double> Vz_ori;
double * P_dec = NULL;
double * D_dec = NULL;
double * Vx_dec = NULL;
double * Vy_dec = NULL;
double * Vz_dec = NULL;
double * Temp_ori = NULL;
std::vector<double> error_Temp;
std::vector<double> error_est_Temp;


template<class T>
bool halfing_error_T_uniform(const T * P, const T * D, size_t n, const double tau, std::vector<double>& ebs){
	double eb_P = ebs[0];
	double eb_D = ebs[1];
	double R = 287.1;
	double c_1 = 1.0 / R;
	double max_value = 0;;
	int max_index = 0;
	for(int i=0; i<n; i++){
		// error of temperature
		double e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
		double Temp = P[i] / (D[i] * R);
		// print_error("T", Temp, Temp_ori[i], e_T);

		error_est_Temp[i] = e_T;
		error_Temp[i] = Temp - Temp_ori[i];

		if(max_value < error_est_Temp[i]){
			max_value = error_est_Temp[i];
			max_index = i;
		}
	}
	// std::cout << "P = " << P[max_index] << " D = " << D[max_index] << std::endl;
	// std::cout << "eb_P = " << eb_P << " eb_D = " << eb_D << std::endl;
	// std::cout << "coeff_P = " << fabs(P[max_index])*eb_D << " coeff_D = " << fabs(D[max_index])*eb_P << std::endl;
	// std::cout << names[1] << ": max estimated error = " << max_value << ", index = " << max_index << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		auto i = max_index;
		double estimate_error = max_value;
        // double T = c_1 * P[i] / D[i];
        double eb_P = ebs[0];
        double eb_D = ebs[1];
        while(estimate_error > tau){
    		// std::cout << "uniform decrease\n";
            eb_P = eb_P / 1.5;
            eb_D = eb_D / 1.5;
            estimate_error = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
        }
        ebs[0] = eb_P;
        ebs[1] = eb_D;
		return false;
	}
	return true;
}



int main(int argc, char ** argv){

    using T = double;
	int argv_id = 1;
    double target_rel_eb = atof(argv[argv_id++]);
	std::string data_prefix_path = argv[argv_id++];
	std::string data_file_prefix = data_prefix_path + "/data/";
	std::string rdata_file_prefix = data_prefix_path + "/refactor/";

    size_t num_elements = 0;
    P_ori = MGARD::readfile<T>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    D_ori = MGARD::readfile<T>((data_file_prefix + "Density.dat").c_str(), num_elements);
    std::vector<double> ebs;
    ebs.push_back(compute_value_range(P_ori)*target_rel_eb);
    ebs.push_back(compute_value_range(D_ori)*target_rel_eb);
	int n_variable = ebs.size();
    std::vector<std::vector<T>> vars_vec = {P_ori, D_ori};
    std::vector<double> var_range(n_variable);
    for(int i=0; i<n_variable; i++){
        var_range[i] = compute_value_range(vars_vec[i]);
    } 

	struct timespec start, end;
	int err;
	double elapsed_time;

	err = clock_gettime(CLOCK_REALTIME, &start);

    std::vector<T> Temp(num_elements);
    compute_T(P_ori.data(), D_ori.data(), num_elements, Temp.data());
	Temp_ori = Temp.data();
    double tau = compute_value_range(Temp)*target_rel_eb;

    std::vector<MDR::ComposedReconstructor<T, MGARDHierarchicalDecomposer<T>, DirectInterleaver<T>, PerBitBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + varlist[i+3];
        std::string metadata_file = rdir_prefix + "_refactored/metadata.bin";
        std::vector<std::string> files;
        int num_levels = 9;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        auto decomposer = MGARDHierarchicalDecomposer<T>();
        auto interleaver = DirectInterleaver<T>();
        auto encoder = PerBitBPEncoder<T, uint32_t>();
        auto compressor = AdaptiveLevelCompressor(64);
        auto estimator = MaxErrorEstimatorHB<T>();
        auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>(estimator);
        auto retriever = ConcatLevelFileRetriever(metadata_file, files);
        reconstructors.push_back(generateReconstructor<T>(decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever));
        reconstructors.back().load_metadata();
    }    
    std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<double>(num_elements));
	std::vector<size_t> total_retrieved_size(n_variable, 0);

    int iter = 0;
    int max_iter = 5;
	bool tolerance_met = false;
	double max_act_error = 0, max_est_error = 0;
    while((!tolerance_met) && (iter < max_iter)){
    	iter ++;
	    for(int i=0; i<n_variable; i++){
	        auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
	        memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
	    }
	    P_dec = reconstructed_vars[0].data();
	    D_dec = reconstructed_vars[1].data();
	    // MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
	    // MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
	    error_Temp = std::vector<double>(num_elements);
	    error_est_Temp = std::vector<double>(num_elements);
		// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
	    // MDR::print_vec(ebs);
	    tolerance_met = halfing_error_T_uniform(P_dec, D_dec, num_elements, tau, ebs);
		// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
	    // MDR::print_vec(ebs);
	    // std::cout << names[1] << " requested error = " << tau << std::endl;
	    max_act_error = print_max_abs(names[1] + " error", error_Temp);
	    max_est_error = print_max_abs(names[1] + " error_est", error_est_Temp);   	
    }
	err = clock_gettime(CLOCK_REALTIME, &end);
	elapsed_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000;

	std::cout << "requested error = " << tau << std::endl;
	std::cout << "max_est_error = " << max_est_error << std::endl;
	std::cout << "max_act_error = " << max_act_error << std::endl;
	std::cout << "iter = " << iter << std::endl;
   
   	size_t total_size = std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0);
	double cr = n_variable * num_elements * sizeof(T) * 1.0 / total_size;
	std::cout << "each retrieved size:";
    for(int i=0; i<n_variable; i++){
        std::cout << total_retrieved_size[i] << ", ";
    }
    std::cout << std::endl;
	// MDR::print_vec(total_retrieved_size);
	std::cout << "aggregated cr = " << cr << std::endl;
	printf("elapsed_time = %.6f\n", elapsed_time);

    return 0;
}