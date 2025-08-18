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
double * C_ori = NULL;
std::vector<double> error_C;
std::vector<double> error_est_C;


template<class T>
bool halfing_error_C_uniform(const T * P, const T * D, size_t n, const double tau, std::vector<double>& ebs){
	double eb_P = ebs[0];
	double eb_D = ebs[1];
	double R = 287.1;
	double gamma = 1.4;
	double c_1 = 1.0 / R;
	double c_2 = sqrt(gamma * R);
	double max_value = 0;
	int max_index = 0;
	for(int i=0; i<n; i++){
		double e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
		double Temp = P[i] / (D[i] * R);
		// error of C
		double e_C = c_2*compute_bound_square_root_x(Temp, e_T);
		double C = c_2 * sqrt(Temp);

		error_est_C[i] = e_C;
		error_C[i] = C - C_ori[i];

		if(max_value < error_est_C[i]){
			max_value = error_est_C[i];
			max_index = i;
		}

	}
	// std::cout << names[2] << ": max estimated error = " << max_value << ", index = " << max_index << std::endl;
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
			estimate_error = c_2*compute_bound_square_root_x(Temp, e_T);
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

    std::vector<T> C(num_elements);
    compute_C(P_ori.data(), D_ori.data(), num_elements, C.data());
	C_ori = C.data();
    double tau = compute_value_range(C)*target_rel_eb;

    std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<double>(num_elements));
	std::vector<size_t> total_retrieved_size(n_variable, 0);

    int iter = 0;
    int max_iter = 30;
	bool tolerance_met = false;
	double max_act_error = 0, max_est_error = 0;
    std::vector<int> current_ind(n_variable, -1);
    T * reconstructed_data = (T *) malloc(num_elements * sizeof(T));
    while((!tolerance_met) && (iter < max_iter)){
        iter ++;
        for(int i=0; i<n_variable; i++){
            std::string rdir_prefix = rdata_file_prefix + varlist[i+3];
            double file_eb = 0.1;
            auto file_ind = find_index(ebs[i]/var_range[i], file_eb);
            // std::cout << "file_ind = " << file_ind << std::endl;
            // std::cout << "Requested relative tolerance = " << ebs[i]/var_range[i] << ", expected tolerance = " << file_eb << "\n"; 
            // std::cout << "Requested tolerance = " << ebs[i] << ", expected tolerance = " << file_eb * var_range[i] << "\n";
            if(file_ind > current_ind[i]){
                for(int j=current_ind[i]+1; j<=file_ind; j++){
                    std::string filename = rdir_prefix + "_refactored/SZ3_delta_eb_" + std::to_string(j) + ".bin";
                    size_t n = 0;
                    auto cmpData = MGARD::readfile<char>(filename.c_str(), n);
                    total_retrieved_size[i] += n;
                    SZ3_decompress(cmpData.data(), n, reconstructed_data);
					for(int j=0; j<num_elements; j++){
						reconstructed_vars[i][j] += reconstructed_data[j];
					}
                }
                current_ind[i] = file_ind;            
            }
        }
	    P_dec = reconstructed_vars[0].data();
	    D_dec = reconstructed_vars[1].data();
	    // MGARD::print_statistics(P_ori.data(), P_dec, num_elements);
	    // MGARD::print_statistics(D_ori.data(), D_dec, num_elements);
	    error_C = std::vector<double>(num_elements);
	    error_est_C = std::vector<double>(num_elements);
		// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
	    // MDR::print_vec(ebs);
	    tolerance_met = halfing_error_C_uniform(P_dec, D_dec, num_elements, tau, ebs);
		// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
	    // MDR::print_vec(ebs);
	    // std::cout << names[2] << " requested error = " << tau << std::endl;
	    max_est_error = print_max_abs(names[2] + " error_est", error_est_C);
	    max_act_error = print_max_abs(names[2] + " actual error", error_C);
    }
	err = clock_gettime(CLOCK_REALTIME, &end);
	elapsed_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000;

	std::cout << "requested_error = " << tau << std::endl;
	std::cout << "max_est_error = " << max_est_error << std::endl;
	std::cout << "max_act_error = " << max_act_error << std::endl;
    free(reconstructed_data);
	std::cout << "iter = " << iter << std::endl;
	std::cout << "each retrieved size:";
    for(int i=0; i<n_variable; i++){
        std::cout << total_retrieved_size[i] << ", ";
    }
    std::cout << std::endl;
	// MDR::print_vec(total_retrieved_size);
	size_t total_size = std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), size_t(0));
	double cr = n_variable * num_elements * sizeof(T) * 1.0 / total_size;
	std::cout << "aggregated cr = " << cr << std::endl;
	std::cout << "bitrate = " << ((sizeof(T) * 8) / cr) << std::endl;
	printf("elapsed_time = %.6f\n", elapsed_time);

    return 0;
}