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

const std::vector<std::string> var_name_out{"VelocityX", "VelocityY", "VelocityZ"};

std::vector<float> Vx_ori;
std::vector<float> Vy_ori;
std::vector<float> Vz_ori;
float * Vx_dec = NULL;
float * Vy_dec = NULL;
float * Vz_dec = NULL;
float * V_TOT_ori = NULL;
std::vector<double> error_V_TOT;
std::vector<double> error_est_V_TOT;

template<class T>
bool halfing_error_V_TOT_uniform(const T * Vx, const T * Vy, const T * Vz, size_t n, const double tau, std::vector<double>& ebs){
	double eb_Vx = ebs[0];
	double eb_Vy = ebs[1];
	double eb_Vz = ebs[2];
	double max_value = 0;
	int max_index = 0;
	for(int i=0; i<n; i++){
		// error of total velocity square
		double e_V_TOT_2 = 0;
		e_V_TOT_2 = compute_bound_x_square((double) Vx[i], eb_Vx) + compute_bound_x_square((double) Vy[i], eb_Vy) + compute_bound_x_square((double) Vz[i], eb_Vz);
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		// error of total velocity
		double e_V_TOT = 0;
		e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
		double V_TOT = sqrt(V_TOT_2);
		// print_error("V_TOT", V_TOT, V_TOT_ori[i], e_V_TOT);

		error_est_V_TOT[i] = e_V_TOT;
		error_V_TOT[i] = V_TOT - V_TOT_ori[i];

		if(max_value < error_est_V_TOT[i]){
			max_value = error_est_V_TOT[i];
			max_index = i;
		}

	}
	// estimate error bound based on maximal errors
	if(max_value > tau){
		// estimate
		auto i = max_index;
		double estimate_error = max_value;
		double V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		double V_TOT = sqrt(V_TOT_2);
		double eb_Vx = ebs[0];
		double eb_Vy = ebs[1];
		double eb_Vz = ebs[2];
		while(estimate_error > tau){
			// change error bound
			eb_Vx = eb_Vx / 1.5;
			eb_Vy = eb_Vy / 1.5;
			eb_Vz = eb_Vz / 1.5; 							        		
			double e_V_TOT_2 = compute_bound_x_square((double) Vx[i], eb_Vx) + compute_bound_x_square((double) Vy[i], eb_Vy) + compute_bound_x_square((double) Vz[i], eb_Vz);
			// double e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
			estimate_error = compute_bound_square_root_x<double>(V_TOT_2, e_V_TOT_2);
		}
		ebs[0] = eb_Vx;
		ebs[1] = eb_Vy;
		ebs[2] = eb_Vz;
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

    using T = float;
	int argv_id = 1;
    double target_rel_eb = atof(argv[argv_id++]);
	std::string data_prefix_path = argv[argv_id++];
	std::string output_path = argv[argv_id++];

	data_prefix_path += oss.str();
	std::string data_file_prefix = data_prefix_path + "/data/";
	std::string rdata_file_prefix = data_prefix_path + "/refactor/";
	int exp = static_cast<int>(std::round(std::log10(target_rel_eb)));
	std::string wdata_file_prefix = output_path + "/1e" + std::to_string(exp) + "/";

    const int target_level = 8;
    // read_file
    size_t num_elements = 0;
	Vx_ori = MGARD::readfile<T>((data_file_prefix + var_name_out[0] + ".dat").c_str(), num_elements);
    Vy_ori = MGARD::readfile<T>((data_file_prefix + var_name_out[1] + ".dat").c_str(), num_elements);
    Vz_ori = MGARD::readfile<T>((data_file_prefix + var_name_out[2] + ".dat").c_str(), num_elements);

    std::vector<double> ebs;
    ebs.push_back(compute_global_value_range(Vx_ori));
    ebs.push_back(compute_global_value_range(Vy_ori));
    ebs.push_back(compute_global_value_range(Vz_ori));
	int n_variable = ebs.size();

    for(int i=0; i<ebs.size(); i++){
    	ebs[i] *= target_rel_eb;
    }

    std::vector<T> V_TOT(num_elements);
    compute_VTOT(Vx_ori.data(), Vy_ori.data(), Vz_ori.data(), num_elements, V_TOT.data());

    double tau, global_tau;
	global_tau = compute_global_value_range(V_TOT)*target_rel_eb;
	tau = global_tau;
	V_TOT_ori = V_TOT.data();

    std::vector<double> value_range(n_variable);
    value_range[0] = compute_value_range(Vx_ori);
    value_range[1] = compute_value_range(Vy_ori);
    value_range[2] = compute_value_range(Vz_ori);

    int iter = 0;
    int max_iter = 30;
	bool tolerance_met = false;
    size_t local_total_size = 0;
	double local_elapsed_time = 0, global_elapsed_time = 0;
    std::vector<int> current_ind(n_variable, -1);
    std::vector<int> prev_ind(n_variable, 0);
	std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
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
					int index = 0;
					for(int j=0; j<num_elements; j++){
						reconstructed_vars[i][j] += reconstructed_data[index ++];
					}
                }
                current_ind[i] = file_ind;
            }
	    }
	    Vx_dec = reconstructed_vars[0].data();
	    Vy_dec = reconstructed_vars[1].data();
	    Vz_dec = reconstructed_vars[2].data();
	    error_V_TOT = std::vector<double>(num_elements);
	    error_est_V_TOT = std::vector<double>(num_elements);
	    tolerance_met = halfing_error_V_TOT_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, tau, ebs);
    }
	local_elapsed_time += MPI_Wtime();
	MPI_Reduce(&local_elapsed_time, &global_elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	// std::cout << "rank = " << rank << " act_iter = " << iter << std::endl;

    free(reconstructed_data);

    if(!rank) printf("requested_error = %.10f\n", global_tau);

	double local_max_est_error = print_max_abs("", error_est_V_TOT);
	double global_max_est_error = 0;
	MPI_Reduce(&local_max_est_error, &global_max_est_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) printf("max_est_error = %.10f\n", global_max_est_error);

	double local_max_error = 0;
	local_max_error = print_max_abs("", error_V_TOT);
	double global_max_act_error = 0;
	MPI_Reduce(&local_max_error, &global_max_act_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) printf("max_act_error = %.10f\n", global_max_act_error);

	unsigned long long int global_total_num = 0;
	MPI_Reduce(&num_elements, &global_total_num, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	unsigned long long int global_total_retrieved = 0;
	MPI_Reduce(&local_total_size, &global_total_retrieved, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	if(!rank) printf("Aggregated bitrate = %.10f, retrieved_size = %ld, total_num_elements = %ld\n", 8*global_total_retrieved * 1.0 / (global_total_num * n_variable), global_total_retrieved, global_total_num);
	if(!rank) printf("elapsed_time = %.6f\n", global_elapsed_time);
	std::vector<char> data_buffer(num_elements*sizeof(double));
	// printf("Processor %d: ???\n", rank);
    for(int i=0; i<n_variable; i++){
	    std::string rdir_prefix = rdata_file_prefix + var_name_out[i] + "_refactored/";
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
		std::string filename = wdata_file_prefix + var_name_out[i] + "_aggregated.dat";
		MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
		MPI_File_write_at(file, offset, data_buffer.data(), local_size, MPI_SIGNED_CHAR, MPI_STATUS_IGNORE);
		MPI_File_close(&file);
		// printf("Processor %d: offset = %ld, count = %ld, offset + count = %ld\n", rank, offset, local_size, local_size + offset);
    }
	//printf("Processor %d total size = %d\n", rank, total_2);
    MPI_Finalize();
    return 0;
}
