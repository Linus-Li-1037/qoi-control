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
#define BLOCK_SIZE 256

using namespace MDR;

template <class T>
__device__ inline T compute_bound_x_square(T x, T eb){
	return 2*fabs(x)*eb + eb*eb;
}

// f(x) = sqrt(x)
template <class T>
__device__ inline T compute_bound_square_root_x(T x, T eb){
	if(x == 0) {
		return sqrt(eb);
	}
	if(x > eb){
		return eb / (sqrt(x - eb) + sqrt(x));
	}
	else{
		return eb / sqrt(x);
	}
}

template<class T>
__global__ void V_TOT_error_estimation(const T *Vx, const T * Vy, const T * Vz, size_t n, const unsigned char * mask, const T tau, T * ebs, T *error_est_V_TOT, T *error_V_TOT, T *V_TOT_ori, int *tolerance_exceed_flag){
	T eb_Vx = ebs[0];
	T eb_Vy = ebs[1];
	T eb_Vz = ebs[2];
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for(int i = tid; i < n; i += stride){
		if (atomicAdd(tolerance_exceed_flag, 0)) return;
		T e_V_TOT_2 = 0;
		if(mask[i]) e_V_TOT_2 = ::compute_bound_x_square(Vx[i], eb_Vx) + ::compute_bound_x_square(Vy[i], eb_Vy) + ::compute_bound_x_square(Vz[i], eb_Vz);
		T V_TOT_2 = Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i];
		// error of total velocity
		T e_V_TOT = 0;
		if(mask[i]) e_V_TOT = ::compute_bound_square_root_x(V_TOT_2, e_V_TOT_2);
		T V_TOT = sqrt(V_TOT_2);
		// print_error("V_TOT", V_TOT, V_TOT_ori[i], e_V_TOT);

		error_est_V_TOT[i] = e_V_TOT;
		error_V_TOT[i] = V_TOT - V_TOT_ori[i];
		if (e_V_TOT > tau){
			atomicExch(tolerance_exceed_flag, 1);
			return;
		}
	}
	return;
}

int main(int argc, char ** argv){

    using T = float;
	int argv_id = 1;
    T target_rel_eb = atof(argv[argv_id++]);
	std::string data_prefix_path = argv[argv_id++];
	std::string data_file_prefix = data_prefix_path + "/data/";
	std::string rdata_file_prefix = data_prefix_path + "/refactor/";
	std::vector<T> P_ori;
	std::vector<T> D_ori;
	std::vector<T> Vx_ori;
	std::vector<T> Vy_ori;
	std::vector<T> Vz_ori;
	T * P_dec = NULL;
	T * D_dec = NULL;
	T * Vx_dec = NULL;
	T * Vy_dec = NULL;
	T * Vz_dec = NULL;
	T * V_TOT_ori = NULL;
	std::vector<T> error_V_TOT;
	std::vector<T> error_est_V_TOT;
    size_t num_elements = 0;
    Vx_ori = MGARD::readfile<T>((data_file_prefix + "VelocityX.dat").c_str(), num_elements);
    Vy_ori = MGARD::readfile<T>((data_file_prefix + "VelocityY.dat").c_str(), num_elements);
    Vz_ori = MGARD::readfile<T>((data_file_prefix + "VelocityZ.dat").c_str(), num_elements);
    std::vector<T> ebs(3, 0);
	int n_variable = ebs.size();
    std::vector<std::vector<T>> vars_vec = {Vx_ori, Vy_ori, Vz_ori};

	struct timespec start, end;
	int err;
	double elapsed_time;

	err = clock_gettime(CLOCK_REALTIME, &start);

    std::vector<T> V_TOT(num_elements);
    compute_VTOT(Vx_ori.data(), Vy_ori.data(), Vz_ori.data(), num_elements, V_TOT.data());
	V_TOT_ori = V_TOT.data();
    T tau = compute_value_range(V_TOT)*target_rel_eb;

    std::string mask_file = rdata_file_prefix + "mask.bin";
    size_t num_valid_data = 0;
    auto mask = MGARD::readfile<unsigned char>(mask_file.c_str(), num_valid_data);
    std::vector<MDR::SegmentedReconstructor<T, MGARDHierarchicalDecomposer<T>, DirectInterleaver<T>, PerBitBPEncoder<T, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<T>>, MaxErrorEstimatorHB<T>, ConcatLevelFileRetriever>> reconstructors;
    for(int i=0; i<n_variable; i++){
        std::string rdir_prefix = rdata_file_prefix + varlist[i];
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
        reconstructors.push_back(generateSegmentedReconstructor<T>(decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever));
        reconstructors.back().load_metadata();
    }    
    std::vector<std::vector<T>> reconstructed_vars(n_variable, std::vector<T>(num_elements));

	std::vector<size_t> total_retrieved_size(n_variable, 0);

    int iter = 0;
    int max_iter = 100;
	bool tolerance_met = false;
	T max_act_error = 0, max_est_error = 0;
	int tolerance_exceed_flag_h;

	T *Vx_dec_d, *Vy_dec_d, *Vz_dec_d;
	unsigned char *mask_d;
	T *ebs_d;
	T *error_est_V_TOT_d, *error_V_TOT_d, *V_TOT_ori_d;
	int *tolerance_exceed_flag_d;
	cudaMalloc((void**) &Vx_dec_d, num_elements * sizeof(T));
	cudaMalloc((void**) &Vy_dec_d, num_elements * sizeof(T));
	cudaMalloc((void**) &Vz_dec_d, num_elements * sizeof(T));
	cudaMalloc((void**) &mask_d, num_elements * sizeof(unsigned char));
	cudaMalloc((void**) &ebs_d, ebs.size() * sizeof(T));
	cudaMalloc((void**) &error_est_V_TOT_d, num_elements * sizeof(T));
	cudaMalloc((void**) &error_V_TOT_d, num_elements * sizeof(T));
	cudaMalloc((void**) &V_TOT_ori_d, num_elements * sizeof(T));
	cudaMalloc((void**) &tolerance_exceed_flag_d, sizeof(int));

	cudaMemcpy(mask_d, mask.data(), num_elements * sizeof(unsigned char), cudaMemcpyHostToDevice);
	cudaMemcpy(V_TOT_ori_d, V_TOT_ori, num_elements * sizeof(T), cudaMemcpyHostToDevice);

	dim3 block(BLOCK_SIZE);
	dim3 grid((num_elements + BLOCK_SIZE - 1) / BLOCK_SIZE);
	
    while((!tolerance_met) && (iter < max_iter)){
    	iter ++;
		cudaMemset(tolerance_exceed_flag_d, 0, sizeof(int));
	    for(int i=0; i<n_variable; i++){
	        auto reconstructed_data = reconstructors[i].progressive_reconstruct(50000000, -1);
			total_retrieved_size[i] = reconstructors[i].get_retrieved_size();
			ebs[i] = reconstructors[i].get_error_bound();
	        if(i < 3){
	            // reconstruct with mask
	            int index = 0;
	            for(int j=0; j<num_elements; j++){
	                if(mask[j]){
	                    reconstructed_vars[i][j] = reconstructed_data[index ++];
	                }
	                else reconstructed_vars[i][j] = 0;
	            }
	        }
	        else{
	            memcpy(reconstructed_vars[i].data(), reconstructed_data, num_elements*sizeof(T));
	        }
	    }

	    Vx_dec = reconstructed_vars[0].data();
	    Vy_dec = reconstructed_vars[1].data();
	    Vz_dec = reconstructed_vars[2].data();
	    // MGARD::print_statistics(Vx_ori.data(), Vx_dec, num_elements);
	    // MGARD::print_statistics(Vy_ori.data(), Vy_dec, num_elements);
	    // MGARD::print_statistics(Vz_ori.data(), Vz_dec, num_elements);
	    error_V_TOT = std::vector<T>(num_elements);
	    error_est_V_TOT = std::vector<T>(num_elements);
		// std::cout << "iter" << iter << ": The old ebs are:" << std::endl;
	    // MDR::print_vec(ebs);
		cudaMemcpy(Vx_dec_d, Vx_dec, num_elements * sizeof(T), cudaMemcpyHostToDevice);
		cudaMemcpy(Vy_dec_d, Vy_dec, num_elements * sizeof(T), cudaMemcpyHostToDevice);
		cudaMemcpy(Vz_dec_d, Vz_dec, num_elements * sizeof(T), cudaMemcpyHostToDevice);
		cudaMemcpy(ebs_d, ebs.data(), ebs.size() * sizeof(T), cudaMemcpyHostToDevice);
	    V_TOT_error_estimation<T><<<grid, block>>>(Vx_dec_d, Vy_dec_d, Vz_dec_d, num_elements, mask_d, tau, ebs_d, error_est_V_TOT_d, error_V_TOT_d, V_TOT_ori_d, tolerance_exceed_flag_d);
		cudaError_t err = cudaGetLastError();
		if (err != cudaSuccess) {
			std::cerr << "CUDA kernel failed: " << cudaGetErrorString(err) << std::endl;
			exit(EXIT_FAILURE);
		}
		cudaMemcpy(&tolerance_exceed_flag_h, tolerance_exceed_flag_d, sizeof(int), cudaMemcpyDeviceToHost);
		tolerance_met = (tolerance_exceed_flag_h == 0);
		// std::cout << "iter" << iter << ": The new ebs are:" << std::endl;
	    // MDR::print_vec(ebs);
	    // std::cout << names[0] << " requested error = " << tau << std::endl;
	    if(tolerance_met){
			cudaMemcpy(error_V_TOT.data(), error_V_TOT_d, num_elements * sizeof(T), cudaMemcpyDeviceToHost);
			cudaMemcpy(error_est_V_TOT.data(), error_est_V_TOT_d, num_elements * sizeof(T), cudaMemcpyDeviceToHost);
			max_act_error = print_max_abs(names[0] + " error", error_V_TOT);
	    	max_est_error = print_max_abs(names[0] + " error_est", error_est_V_TOT);  
		}
    }
	err = clock_gettime(CLOCK_REALTIME, &end);
	elapsed_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000;

	std::cout << "requested error = " << tau << std::endl;
	std::cout << "max_est_error = " << max_est_error << std::endl;
	std::cout << "max_act_error = " << max_act_error << std::endl;
	std::cout << "iter = " << iter << std::endl;

	size_t total_size = std::accumulate(total_retrieved_size.begin(), total_retrieved_size.end(), 0);
	// std::cout << "total_size: " << total_size << std::endl;
	// std::cout << "original_size: " << n_variable * num_elements * sizeof(T) << std::endl;
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