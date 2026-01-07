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

const std::vector<std::string> var_name_out{"VelocityX", "VelocityY", "VelocityZ", "Pressure", "Density"};

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
double * PT_ori = NULL;
std::vector<double> error_PT;
std::vector<double> error_est_PT;

template<class T>
bool halfing_error_PT_uniform(const T * Vx, const T * Vy, const T * Vz, const T * P, const T * D, size_t n, const std::vector<unsigned char>& mask, const double tau, std::vector<double>& ebs){
	double eb_Vx = ebs[0];
	double eb_Vy = ebs[1];
	double eb_Vz = ebs[2];
	double eb_P = ebs[3];
	double eb_D = ebs[4];
	double R = 287.1;
	double gamma = 1.4;
	double mi = 3.5;
	double mu_r = 1.716e-5;
	double T_r = 273.15;
	double S = 110.4;
	double c_1 = 1.0 / R;
	double c_2 = sqrt(gamma * R);
	int C7i[8] = {1, 7, 21, 35, 35, 21, 7, 1};
	double max_value = 0;
	int max_index = 0;
	int n_variable = ebs.size();
	double Mach_tmp_pow[8];
    double e_Mach_tmp_pow[8];
	for(int i=0; i<n; i++){
		double e_V_TOT_2 = 0;
		e_V_TOT_2 = mask[i] ? compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz) : 0;
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		double e_V_TOT = 0;
		e_V_TOT = mask[i] ? compute_bound_square_root_x(V_TOT_2, e_V_TOT_2) : 0;
		double V_TOT = sqrt(V_TOT_2);
		double e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
		double Temp = P[i] / (D[i] * R);
		double e_C = c_2*compute_bound_square_root_x(Temp, e_T);
		double C = c_2 * sqrt(Temp);
		double e_Mach = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
		double Mach = V_TOT / C;
		double e_Mach_tmp = ldexp(gamma - 1, -1) * compute_bound_x_square(Mach, e_Mach);
		double Mach_tmp = 1 + ldexp(gamma - 1, -1) * Mach * Mach;
		double e_Mach_tmp_mi = 0;
        Mach_tmp_pow[0] = 1;
        e_Mach_tmp_pow[0] = 1;
        for (int k = 1; k <= 7; k++) {
            Mach_tmp_pow[k] = Mach_tmp_pow[k - 1] * Mach_tmp;
            e_Mach_tmp_pow[k] = e_Mach_tmp_pow[k - 1] * e_Mach_tmp;
        }
        for (int k = 1; k <= 7; k++) {
            e_Mach_tmp_mi += C7i[k] * Mach_tmp_pow[7 - k] * e_Mach_tmp_pow[k];
        }
		double Mach_tmp_mi = sqrt(Mach_tmp * Mach_tmp * Mach_tmp * Mach_tmp * Mach_tmp * Mach_tmp * Mach_tmp);
		double e_PT = compute_bound_multiplication(P[i], Mach_tmp_mi, eb_P, e_Mach_tmp_mi);
		double PT = P[i] * Mach_tmp_mi;

		error_est_PT[i] = e_PT;
		error_PT[i] = PT - PT_ori[i];
		if(max_value < error_est_PT[i]){
			max_value = error_est_PT[i];
			max_index = i;
		}
	}
	// std::cout << names[3] << ": max estimated error = " << max_value << ", index = " << max_index << std::endl;
	// estimate error bound based on maximal errors
	if(max_value > tau){
		auto i = max_index;
		double estimate_error = max_value;
		double eb_Vx = ebs[0];
		double eb_Vy = ebs[1];
		double eb_Vz = ebs[2];
		double eb_P = ebs[3];
		double eb_D = ebs[4];
		while(estimate_error > tau){
    		// std::cout << "uniform decrease\n";
			eb_Vx = eb_Vx / 1.5;
			eb_Vy = eb_Vy / 1.5;
			eb_Vz = eb_Vz / 1.5; 
			eb_P = eb_P / 1.5;
			eb_D = eb_D / 1.5;
			double e_V_TOT_2 = 0;
			if(mask[i]) e_V_TOT_2 = compute_bound_x_square(Vx[i], eb_Vx) + compute_bound_x_square(Vy[i], eb_Vy) + compute_bound_x_square(Vz[i], eb_Vz);
			double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
			double e_V_TOT = 0;
			if(mask[i]) e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
			double V_TOT = sqrt(V_TOT_2);
			double e_T = c_1 * compute_bound_division(P[i], D[i], eb_P, eb_D);
			double Temp = P[i] / (D[i] * R);
			double e_C = c_2*compute_bound_square_root_x(Temp, e_T);
			double C = c_2 * sqrt(Temp);
			double e_Mach = compute_bound_division(V_TOT, C, e_V_TOT, e_C);
			double Mach = V_TOT / C;
			double e_Mach_tmp = (gamma-1) / 2 * compute_bound_x_square(Mach, e_Mach);
			double Mach_tmp = 1 + (gamma-1)/2 * Mach * Mach;
			double e_Mach_tmp_mi = 0;
			Mach_tmp_pow[0] = 1;
			e_Mach_tmp_pow[0] = 1;
			for (int k = 1; k <= 7; k++) {
				Mach_tmp_pow[k] = Mach_tmp_pow[k - 1] * Mach_tmp;
				e_Mach_tmp_pow[k] = e_Mach_tmp_pow[k - 1] * e_Mach_tmp;
			}
			for (int k = 1; k <= 7; k++) {
				e_Mach_tmp_mi += C7i[k] * Mach_tmp_pow[7 - k] * e_Mach_tmp_pow[k];
			}
			double Mach_tmp_mi = sqrt(Mach_tmp * Mach_tmp * Mach_tmp * Mach_tmp * Mach_tmp * Mach_tmp * Mach_tmp);
			estimate_error = compute_bound_multiplication(P[i], Mach_tmp_mi, eb_P, e_Mach_tmp_mi);
			// if((ebs[0] / eb_Vx) > 10) break;
		}
		ebs[0] = eb_Vx;
		ebs[1] = eb_Vy;
		ebs[2] = eb_Vz;
		ebs[3] = eb_P;
		ebs[4] = eb_D;
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
	int exp = static_cast<int>(std::round(std::log10(target_rel_eb)));
	std::string wdata_file_prefix = output_path + "/1e" + std::to_string(exp) + "/";

    size_t num_elements = 0;
    P_ori = MGARD::readfile<T>((data_file_prefix + "Pressure.dat").c_str(), num_elements);
    D_ori = MGARD::readfile<T>((data_file_prefix + "Density.dat").c_str(), num_elements);
    Vx_ori = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    Vy_ori = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    Vz_ori = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    std::vector<double> ebs;
    ebs.push_back(compute_value_range(Vx_ori)*target_rel_eb);
    ebs.push_back(compute_value_range(Vy_ori)*target_rel_eb);
    ebs.push_back(compute_value_range(Vz_ori)*target_rel_eb);
    ebs.push_back(compute_value_range(P_ori)*target_rel_eb);
    ebs.push_back(compute_value_range(D_ori)*target_rel_eb);
	int n_variable = ebs.size();
    std::vector<std::vector<T>> vars_vec = {Vx_ori, Vy_ori, Vz_ori, P_ori, D_ori};
    std::vector<double> var_range(n_variable);
    for(int i=0; i<n_variable; i++){
        var_range[i] = compute_value_range(vars_vec[i]);
    } 

    std::vector<T> PT(num_elements);
    compute_PT(Vx_ori.data(), Vy_ori.data(), Vz_ori.data(), P_ori.data(), D_ori.data(), num_elements, PT.data());

    double tau, global_tau;
	global_tau = compute_global_value_range(PT)*target_rel_eb;
	tau = global_tau;
	PT_ori = PT.data();

    std::vector<double> value_range(n_variable);
    value_range[0] = compute_value_range(Vx_ori);
    value_range[1] = compute_value_range(Vy_ori);
    value_range[2] = compute_value_range(Vz_ori);
    value_range[3] = compute_value_range(P_ori);
    value_range[4] = compute_value_range(D_ori);
    std::string mask_file = rdata_file_prefix + "/mask.bin";
	uint32_t mask_file_size = 0;
    auto mask = readmask(mask_file.c_str(), mask_file_size);

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
					if(i < 3){
						// reconstruct with mask
						for(int j=0; j<num_elements; j++){
							if(mask[j]){
								reconstructed_vars[i][j] += reconstructed_data[j];
							}
							else reconstructed_vars[i][j] = 0;
						}
					}
					else{
						for(int j=0; j<num_elements; j++){
							reconstructed_vars[i][j] += reconstructed_data[j];
						}						
					}
                }
                current_ind[i] = file_ind;
            }
	    }
	    Vx_dec = reconstructed_vars[0].data();
	    Vy_dec = reconstructed_vars[1].data();
	    Vz_dec = reconstructed_vars[2].data();
        P_dec = reconstructed_vars[3].data();
	    D_dec = reconstructed_vars[4].data();
	    error_PT = std::vector<double>(num_elements);
	    error_est_PT = std::vector<double>(num_elements);
	    tolerance_met = halfing_error_PT_uniform(Vx_dec, Vy_dec, Vz_dec, P_dec, D_dec, num_elements, mask, tau, ebs);
    }
	local_elapsed_time += MPI_Wtime();
	MPI_Reduce(&local_elapsed_time, &global_elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	int global_max_iter = 0;
    MPI_Reduce(&iter, &global_max_iter, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    if(!rank) std::cout << "max_iter = " << global_max_iter << std::endl;

    free(reconstructed_data);

    if(!rank) std::cout << "requested_error = " << global_tau << std::endl;

	double local_max_est_error = print_max_abs("", error_est_PT);
	double global_max_est_error = 0;
	MPI_Reduce(&local_max_est_error, &global_max_est_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) std::cout << "max_est_error = " << global_max_est_error << std::endl;

	double local_max_error = 0;
	local_max_error = print_max_abs("", error_PT);
	double global_max_act_error = 0;
	MPI_Reduce(&local_max_error, &global_max_act_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) std::cout << "max_act_error = " << global_max_act_error << std::endl;

    local_total_size += mask_file_size;

	unsigned long long int global_total_num = 0;
	MPI_Reduce(&num_elements, &global_total_num, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	unsigned long long int global_total_retrieved = 0;
	MPI_Reduce(&local_total_size, &global_total_retrieved, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	if(!rank) printf("Aggregated bitrate = %.10f, retrieved_size = %ld, total_num_elements = %ld\n", 8 * global_total_retrieved * 1.0 / (global_total_num * n_variable), global_total_retrieved, global_total_num);
	if(!rank) printf("elapsed_time = %.6f\n", global_elapsed_time);

	std::vector<char> data_buffer(num_elements*sizeof(double));
	// printf("Processor %d: ???\n", rank);
    for(int i=0; i<n_variable; i++){
	    std::string rdir_prefix = rdata_file_prefix + var_name_out[i] + "_refactored/";;
    	// copy local data
    	size_t local_size = 0;
    	char * data_buffer_pos = data_buffer.data();
    	for(int j=prev_ind[i]; j<=current_ind[i]; j++){
			size_t n = 0;
			std::string filename = rdir_prefix + "SZ3_delta_eb_" + std::to_string(j) + ".bin";
			auto cmpData = MGARD::readfile<char>(filename.c_str(), n); 
			memcpy(data_buffer_pos, cmpData.data(), n);	
			data_buffer_pos += n;
			local_size += n;
    	}
		// printf("Processor %d: local_size = %ld, total_retrieved_size = %ld, total_size = %ld\n", rank, local_size, total_retrieved_sizes[i], total_size);
    	unsigned long long int count = local_size;
		unsigned long long int offset = 0;
		unsigned long long int buffer = 0;
		for(int j=0; j<size; j++){
			if(j == rank){
				if(j != 0) {
					MPI_Recv(&offset, 1, MPI_UNSIGNED_LONG_LONG, j-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
				buffer = offset + count;
				if(j != size - 1) MPI_Send(&buffer, 1, MPI_UNSIGNED_LONG_LONG, j+1, 0, MPI_COMM_WORLD);
			}
		}
		MPI_File file;
		std::string filename = wdata_file_prefix + "/SZ3_delta_" + var_name_out[i] + "_aggregated.dat";
		MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
		MPI_File_write_at(file, offset, data_buffer.data(), local_size, MPI_SIGNED_CHAR, MPI_STATUS_IGNORE);
		MPI_File_close(&file);
		// printf("Processor %d: offset = %ld, count = %ld, offset + count = %ld\n", rank, offset, local_size, local_size + offset);
		if(i == 0){
            // mask file size already known
			unsigned long long int mask_offset = 0;
			unsigned long long int mask_buffer; 
			for(int j=0; j<size; j++){
				if(j == rank){
					if(j != 0) {
						MPI_Recv(&mask_offset, 1, MPI_UNSIGNED_LONG_LONG, j-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					}
					mask_buffer = mask_offset + mask_file_size;
					if(j != size - 1) {
						MPI_Send(&mask_buffer, 1, MPI_UNSIGNED_LONG_LONG, j+1, 0, MPI_COMM_WORLD);
					}
				}
			}
			size_t mask_num_char = 0;
			auto mask_data = MGARD::readfile<unsigned char> (mask_file.c_str(), mask_num_char);
			MPI_File mask_file;
			std::string mask_filename = wdata_file_prefix + "aggregated_mask.bin";
			MPI_File_open(MPI_COMM_WORLD, mask_filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mask_file);
			MPI_File_write_at(mask_file, mask_offset, mask_data.data(), mask_file_size, MPI_SIGNED_CHAR, MPI_STATUS_IGNORE);
			MPI_File_close(&mask_file);
        }
    }
	//printf("Processor %d total size = %d\n", rank, total_2);

    MPI_Finalize();
    return 0;
}
