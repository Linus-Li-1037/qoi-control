#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include<limits>
#include "mpi.h"
#include <sstream>
#include "utils.hpp"
#include "qoi_utils.hpp"
#include "Reconstructor/Reconstructor.hpp"

const std::vector<std::string> var_name_out{"VelocityX", "VelocityY", "VelocityZ"};
const int n_vars = 3;

using namespace MDR;

std::vector<float> Vx_ori;
std::vector<float> Vy_ori;
std::vector<float> Vz_ori;
float * Vx_dec = NULL;
float * Vy_dec = NULL;
float * Vz_dec = NULL;
float * V_TOT_ori = NULL;
std::vector<double> error_V_TOT;
std::vector<double> error_est_V_TOT;


template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
MDR::ComposedReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever> generateReconstructor(Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::ComposedReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    return reconstructor;
}

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
	// std::cout << names[0] << ": max estimated error = " << max_value << ", index = " << max_index << std::endl;
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
    		// std::cout << "uniform decrease\n";
			eb_Vx = eb_Vx / 1.5;
			eb_Vy = eb_Vy / 1.5;
			eb_Vz = eb_Vz / 1.5; 							        		
			double e_V_TOT_2 = compute_bound_x_square((double) Vx[i], eb_Vx) + compute_bound_x_square((double) Vy[i], eb_Vy) + compute_bound_x_square((double) Vz[i], eb_Vz);
			// float e_V_TOT = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
			estimate_error = compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
		}
		ebs[0] = eb_Vx;
		ebs[1] = eb_Vy;
		ebs[2] = eb_Vz;
		return false;
	}
	return true;
}

template <class T>
T print_max_abs(int rank, const std::string& name, const std::vector<T>& vec){
	T max = fabs(vec[0]);
	for(int i=1; i<vec.size(); i++){
		if(max < fabs(vec[i])) max = fabs(vec[i]);
	}
	// printf("Processor %d var %s: max absolute value =  %.10f\n", rank, name.c_str(), max);
	return max;
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


    const int target_level = 4;
	const int num_levels = target_level + 1;
    // read_file
    size_t num_elements = 0;
	Vx_ori = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    Vy_ori = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    Vz_ori = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);

    std::vector<double> ebs;
    ebs.push_back(compute_value_range(Vx_ori));
    ebs.push_back(compute_value_range(Vy_ori));
    ebs.push_back(compute_value_range(Vz_ori));
	int n_variable = ebs.size();

    for(int i=0; i<ebs.size(); i++){
    	ebs[i] *= target_rel_eb;
    }

    std::vector<T> V_TOT(num_elements);
    compute_VTOT(Vx_ori.data(), Vy_ori.data(), Vz_ori.data(), num_elements, V_TOT.data());

    double tau, global_tau;
	double local_max = -std::numeric_limits<double>::max();
	double local_min = std::numeric_limits<double>::max();
	double global_max = 0, global_min = 0;
	for(int i=0; i<num_elements; i++){
		if(V_TOT[i] > local_max) local_max = V_TOT[i];
		if(V_TOT[i] < local_min) local_min = V_TOT[i];
	}
	MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	double max_value_range = global_max - global_min;
	global_tau = max_value_range * target_rel_eb;
	// tau = (local_max - local_min) * target_rel_eb;
	tau = global_tau;
	V_TOT_ori = V_TOT.data();

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
    std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));
	
    int iter = 0;
    int max_iter = 20;
	bool tolerance_met = false;
    size_t total_size = 0;
	double local_elapsed_time = 0, max_time = 0;
	local_elapsed_time = -MPI_Wtime();
    while((!tolerance_met) && (iter < max_iter)){
			iter ++;
			// std::cout << "*rank " << rank << ", iter " << iter << std::endl;
			total_size = 0;
			for(int i=0; i<n_variable; i++){
				auto reconstructed_data = reconstructors[i].progressive_reconstruct(ebs[i], -1);
			    total_size += reconstructors[i].get_retrieved_size();
				memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
			}
			Vx_dec = reconstructed_vars[0].data();
			Vy_dec = reconstructed_vars[1].data();
			Vz_dec = reconstructed_vars[2].data();
			error_V_TOT = std::vector<double>(num_elements);
			error_est_V_TOT = std::vector<double>(num_elements);
			tolerance_met = halfing_error_V_TOT_uniform(Vx_dec, Vy_dec, Vz_dec, num_elements, tau, ebs);
    }
	for(int i=iter; i<max_iter; i++){
		// std::cout << "#rank " << rank << ", iter " << iter << std::endl;
		for(int i=0; i<n_variable; i++) {
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
	// std::cout << "rank = " << rank << " act_iter = " << iter << std::endl;
	local_elapsed_time += MPI_Wtime();
	MPI_Reduce(&local_elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(!rank) printf("Target V_TOT error = %.10f\n", global_tau);
	double local_IO_time = 0, gloabl_IO_time = 0;
	for(int i=0; i<n_variable; i++){
		local_IO_time += reconstructors[i].get_IO_time();
	}
	MPI_Reduce(&local_IO_time, &gloabl_IO_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) printf("IO_time = %.6f\n", gloabl_IO_time);
	double max_error = 0;
	max_error = print_max_abs(rank, "V_TOT error", error_V_TOT);
	double max_vtot_error = 0;
	MPI_Reduce(&max_error, &max_vtot_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) printf("Max aggregated V_TOT act_error = %.10f\n", max_vtot_error);
	double max_est_error = print_max_abs(rank, "V_TOT error", error_est_V_TOT);
	double max_vtot_error_est = 0;
	MPI_Reduce(&max_est_error, &max_vtot_error_est, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!rank) printf("Max aggregated V_TOT est_error = %.10f\n", max_vtot_error_est);
	unsigned long long int total_num = 0;
	MPI_Reduce(&num_elements, &total_num, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	unsigned long long int total_retrieved = 0;
	MPI_Reduce(&total_size, &total_retrieved, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	if(!rank) printf("Aggregated bitrate = %.10f, retrieved_size = %ld, total_num_elements = %ld\n", 8*total_retrieved * 1.0 / (total_num * n_variable), total_retrieved, total_num);
	if(!rank) printf("elapsed_time = %.6f\n", max_time);
    // size_t total_2 = 0;
    // for(int i=0; i<n_variable; i++){
    //     auto count = reconstructors[i].get_offsets();
    //     auto offsets(count);
    //     for(int j=0; j<offsets.size(); j++) offsets[j] = 0;
    //     auto buffer(offsets);
    //     for(int j=0; j<size; j++){
    //         if(j == rank){
    //             if(j != 0) {
    //                 MPI_Recv(&offsets[0], offsets.size(), MPI_UNSIGNED, j-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //             }
    //             for(int k=0; k<offsets.size(); k++){
    //                 buffer[k] = offsets[k] + count[k];
    //             }
    //             if(j != size - 1) MPI_Send(&buffer[0], offsets.size(), MPI_UNSIGNED, j+1, 0, MPI_COMM_WORLD);
    //         }
    //     }
    //     for(int k=0; k<offsets.size(); k++){
    //         std::string rdir_prefix = rdata_file_prefix + var_name_out[i] + "_refactored/";
    //         std::string file_level = rdir_prefix + "level_" + std::to_string(k) + ".bin";
    //         size_t num_char = 0;
    //         auto level_data = MGARD::readfile<unsigned char>(file_level.c_str(), num_char);
    //         MPI_File file;
    //         std::string filename = wdata_file_prefix + var_name_out[i] + "_aggregated_level_" + std::to_string(k) + ".dat";
    //         MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    //         MPI_File_write_at(file, offsets[k], level_data.data(), count[k], MPI_SIGNED_CHAR, MPI_STATUS_IGNORE);
    //         MPI_File_close(&file);
    //         total_2 += count[k];
    //     }
    // }
    MPI_Finalize();
    return 0;
}
