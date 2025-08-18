#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include <limits>
#include "mpi.h"
#include <sstream>
#include "utils.hpp"
#include "qoi_utils.hpp"
#include "Reconstructor/Reconstructor.hpp"
#include "Synthesizer4GE.hpp"

using namespace MDR;

const std::vector<std::string> var_name_out{"Pressure", "Density"};

std::vector<double> P_ori;
std::vector<double> D_ori;
double * P_dec = NULL;
double * D_dec = NULL;
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

template<class T>
T compute_global_value_range(const std::vector<T> data_vec){
	T global_max = 0, global_min = 0;
	T local_max = -std::numeric_limits<T>::max();
	T local_min = std::numeric_limits<T>::max();
	for(int i=0; i<data_vec.size(); i++){
		if(data_vec[i] > local_max) local_max = data_vec[i];
		if(data_vec[i] < local_min)	local_min = data_vec[i];
	}
	if(std::is_same<T, double>::value){
		MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
	else if(std::is_same<T, float>::value){
		MPI_Allreduce(&local_min, &global_min, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&local_max, &global_max, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
	}
	return global_max - global_min;
}

int main(int argc, char ** argv){
    if(argc != 4){
		std::cout << "usage: para_VTOT [eb] [input path] [output path]" << std::endl;
 	}
	MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	std::ostringstream oss;
	oss << rank;

    using T = double;
	int argv_id = 1;
    double target_rel_eb = atof(argv[argv_id++]);
	std::string data_prefix_path = argv[argv_id++];
    std::string output_path = argv[argv_id++];

	data_prefix_path += oss.str();
	std::string data_file_prefix = data_prefix_path + "/data/";
	std::string rdata_file_prefix = data_prefix_path + "/refactor/";
    // int exp = static_cast<int>(std::round(std::log10(target_rel_eb)));
	// std::string wdata_file_prefix = output_path + "/1e" + std::to_string(exp) + "/";
    
    const int target_level = 4;
	const int num_levels = target_level + 1;

    size_t num_elements = 0;
    P_ori = MGARD::readfile<T>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    D_ori = MGARD::readfile<T>((data_file_prefix + "Density.dat").c_str(), num_elements);
    std::vector<double> ebs;
    ebs.push_back(compute_global_value_range(P_ori)*target_rel_eb);
    ebs.push_back(compute_global_value_range(D_ori)*target_rel_eb);
	int n_variable = ebs.size();

    std::vector<T> Temp(num_elements);
    compute_T(P_ori.data(), D_ori.data(), num_elements, Temp.data());
    
    double tau, global_tau;
	global_tau = compute_global_value_range(Temp) * target_rel_eb;
	// tau = (local_max - local_min) * target_rel_eb;
	tau = global_tau;
    Temp_ori = Temp.data();

    std::vector<MDR::ComposedReconstructor<T, MGARDHierarchicalDecomposer<T>, DirectInterleaver<T>, PerBitBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + var_name_out[i] + "_refactored/";
        std::string metadata_file = rdir_prefix + "metadata.bin";
        std::vector<std::string> files;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "level_" + std::to_string(i) + ".bin";
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

    int iter = 0;
    int max_iter = 30;
	bool tolerance_met = false;
    size_t local_total_size = 0;
	double local_elapsed_time = 0, global_eplapsed_time = 0;
	local_elapsed_time = -MPI_Wtime();
    while((!tolerance_met) && (iter < max_iter)){
    	iter ++;
	    local_total_size = 0;
			for(int i=0; i<n_variable; i++){
				auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
			    local_total_size += reconstructors[i].get_retrieved_size();
                memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			}
	    P_dec = reconstructed_vars[0].data();
	    D_dec = reconstructed_vars[1].data();
	    error_Temp = std::vector<double>(num_elements);
	    error_est_Temp = std::vector<double>(num_elements);
	    tolerance_met = halfing_error_T_uniform(P_dec, D_dec, num_elements, tau, ebs);
    }
    local_elapsed_time += MPI_Wtime();
    MPI_Reduce(&local_elapsed_time, &global_eplapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // std::cout << "rank = " << rank << " act_iter = " << iter << std::endl;
	
    if(!rank) printf("requested_error = %.10f\n", global_tau);

	double local_max_est_error = print_max_abs("", error_est_Temp);
	double global_max_est_error = 0;
	MPI_Reduce(&local_max_est_error, &global_max_est_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) printf("max_est_error = %.10f\n", global_max_est_error);

	double local_max_error = 0;
	local_max_error = print_max_abs("", error_Temp);
	double global_max_act_error = 0;
	MPI_Reduce(&local_max_error, &global_max_act_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) printf("max_act_error = %.10f\n", global_max_act_error);

	unsigned long long int global_total_num = 0;
	MPI_Reduce(&num_elements, &global_total_num, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	unsigned long long int global_total_retrieved = 0;
	MPI_Reduce(&local_total_size, &global_total_retrieved, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	if(!rank) printf("Aggregated bitrate = %.10f, retrieved_size = %ld, total_num_elements = %ld\n", 8*global_total_retrieved * 1.0 / (global_total_num * n_variable), global_total_retrieved, global_total_num);
	if(!rank) printf("elapsed_time = %.6f\n", global_eplapsed_time);
    MPI_Finalize();
    return 0;
}