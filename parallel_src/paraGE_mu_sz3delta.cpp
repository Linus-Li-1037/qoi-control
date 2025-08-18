#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include "mpi.h"
#include <sstream>
#include "adios2.h"
#include "utils.hpp"
#include "qoi_utils.hpp"
#include "SZ3/api/sz.hpp"
#include "Synthesizer4GE.hpp"

using namespace MDR;

const std::vector<std::string> var_name_out{"Pressure", "Density"};

std::vector<double> P_ori;
std::vector<double> D_ori;
double * P_dec = NULL;
double * D_dec = NULL;
double * mu_ori = NULL;
std::vector<double> error_mu;
std::vector<double> error_est_mu;

template<class T>
bool halfing_error_mu_uniform(const T * P, const T * D, size_t n, const double tau, std::vector<double>& ebs){
	double eb_P = ebs[0];
	double eb_D = ebs[1];
	double R = 287.1;
	double gamma = 1.4;
	double mi = 3.5;
	double mu_r = 1.716e-5;
	double T_r = 273.15;
	double S = 110.4;
	double c_1 = 1.0 / R;
	double c_2 = sqrt(gamma * R);
	double c_3 = T_r + S;
	double max_value = 0;
	int max_index = 0;
	int n_variable = ebs.size();
	for(int i=0; i<n; i++){
		double e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
		double Temp = P[i] / (D[i] * R);
		double e_TrS_TS = c_3 * compute_bound_radical(Temp, S, e_T);
		double TrS_TS = c_3 / (Temp + S);
		double e_T_Tr_3 = 3*pow(Temp/T_r, 2)*(e_T/T_r) + 3*Temp/T_r*(e_T/T_r)*(e_T/T_r) + (e_T/T_r)*(e_T/T_r)*(e_T/T_r);
		double T_Tr_3 = pow(Temp/T_r, 3);
		double e_T_Tr_3_sqrt = compute_bound_square_root_x(T_Tr_3, e_T_Tr_3);
		double T_Tr_3_sqrt = sqrt(T_Tr_3);
		double e_mu = mu_r * compute_bound_multiplication(T_Tr_3_sqrt, TrS_TS, e_T_Tr_3_sqrt, e_TrS_TS);
		double mu = mu_r * T_Tr_3_sqrt * TrS_TS;

		error_est_mu[i] = e_mu;
		error_mu[i] = mu - mu_ori[i];

		if(max_value < error_est_mu[i]){
			max_value = error_est_mu[i];
			max_index = i;
		}
	}
	// std::cout << names[3] << ": max estimated error = " << max_value << ", index = " << max_index << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		auto i = max_index;
		double estimate_error = max_value;
		double eb_P = ebs[0];
		double eb_D = ebs[1];
		while(estimate_error > tau){
    		// std::cout << "uniform decrease\n";
			eb_P = eb_P / 1.5;
			eb_D = eb_D / 1.5;
			double e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
			double Temp = P[i] / (D[i] * R);
			double e_TrS_TS = c_3 * compute_bound_radical(Temp, S, e_T);
			double TrS_TS = c_3 / (Temp + S);
			double e_T_Tr_3 = 3*pow(Temp/T_r, 2)*(e_T/T_r) + 3*Temp/T_r*(e_T/T_r)*(e_T/T_r) + (e_T/T_r)*(e_T/T_r)*(e_T/T_r);
			double T_Tr_3 = pow(Temp/T_r, 3);
			double e_T_Tr_3_sqrt = compute_bound_square_root_x(T_Tr_3, e_T_Tr_3);
			double T_Tr_3_sqrt = sqrt(T_Tr_3);
			estimate_error = mu_r * compute_bound_multiplication(T_Tr_3_sqrt, TrS_TS, e_T_Tr_3_sqrt, e_TrS_TS);			
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

    const int target_level = 8;
    // read_file
    size_t num_elements = 0;
    P_ori = MGARD::readfile<T>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    D_ori = MGARD::readfile<T>((data_file_prefix + "Density.dat").c_str(), num_elements);
    std::vector<double> ebs;
    ebs.push_back(compute_global_value_range(P_ori)*target_rel_eb);
    ebs.push_back(compute_global_value_range(D_ori)*target_rel_eb);
	int n_variable = ebs.size();

    std::vector<T> mu(num_elements);
    compute_mu(P_ori.data(), D_ori.data(), num_elements, mu.data());

    double tau, global_tau;
	global_tau = compute_global_value_range(mu) * target_rel_eb;
	// tau = (local_max - local_min) * target_rel_eb;
	tau = global_tau;
    mu_ori = mu.data();

    std::vector<double> value_range(n_variable);
    value_range[0] = compute_value_range(P_ori);
    value_range[1] = compute_value_range(D_ori);

    int iter = 0;
    int max_iter = 30;
	bool tolerance_met = false;
    size_t local_total_size = 0;
	double local_elapsed_time = 0, global_elapsed_time = 0;
    std::vector<int> current_ind(n_variable, -1);
    std::vector<int> prev_ind(n_variable, 0);
	std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<double>(num_elements));
	T * reconstructed_data = (T *) malloc(num_elements * sizeof(T));
	local_elapsed_time = -MPI_Wtime();
    while((!tolerance_met) && (iter < max_iter)){
    	iter ++;
	    for(int i=0; i<n_variable; i++){
	        std::string rdir_prefix = rdata_file_prefix + var_name_out[i] + "_refactored/";
            double rel_eb = ebs[i]/value_range[i];
            double file_eb = 0.1;
            auto file_ind = find_index(rel_eb, file_eb);
	    	if(file_ind > 17) file_ind = 17;
			if(file_ind > current_ind[i]){
                for(int j=current_ind[i]+1; j<=file_ind; j++){
	                std::string filename = rdir_prefix + "SZ3_delta_eb_" + std::to_string(j) + ".bin";
                    size_t n = 0;
                    auto cmpData = MGARD::readfile<char>(filename.c_str(), n);
		    		local_total_size += n;
                    SZ3_decompress(cmpData.data(), n, reconstructed_data);
					// reconstruct with mask
					int index = 0;
					for(int j=0; j<num_elements; j++){
						reconstructed_vars[i][j] += reconstructed_data[j];
					}
                }
                current_ind[i] = file_ind;
            }
	    }
	    P_dec = reconstructed_vars[0].data();
	    D_dec = reconstructed_vars[1].data();
	    error_mu = std::vector<double>(num_elements);
	    error_est_mu = std::vector<double>(num_elements);
	    tolerance_met = halfing_error_mu_uniform(P_dec, D_dec, num_elements, tau, ebs);
    }
    local_elapsed_time += MPI_Wtime();
    MPI_Reduce(&local_elapsed_time, &global_elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    free(reconstructed_data);

	// std::cout << "rank = " << rank << " act_iter = " << iter << std::endl;

    if(!rank) printf("requested_error = %.10f\n", global_tau);

	double local_max_est_error = print_max_abs("", error_est_mu);
	double global_max_est_error = 0;
	MPI_Reduce(&local_max_est_error, &global_max_est_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) printf("max_est_error = %.10f\n", global_max_est_error);

	double local_max_error = 0;
	local_max_error = print_max_abs("", error_mu);
	double global_max_act_error = 0;
	MPI_Reduce(&local_max_error, &global_max_act_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) printf("max_act_error = %.10f\n", global_max_act_error);

	unsigned long long int global_total_num = 0;
	MPI_Reduce(&num_elements, &global_total_num, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	unsigned long long int global_total_retrieved = 0;
	MPI_Reduce(&local_total_size, &global_total_retrieved, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	if(!rank) printf("Aggregated bitrate = %.10f, retrieved_size = %ld, total_num_elements = %ld\n", 8*global_total_retrieved * 1.0 / (global_total_num * n_variable), global_total_retrieved, global_total_num);
	if(!rank) printf("elapsed_time = %.6f\n", global_elapsed_time);
    MPI_Finalize();
    return 0;
}
