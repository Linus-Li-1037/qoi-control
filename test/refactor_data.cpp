#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include "utils.hpp"
#include "Synthesizer4GE.hpp"

using namespace MDR;

int main(int argc, char** argv){

    using T = float;
    int argv_id = 1;
    int mode = atoi(argv[argv_id++]);
    std::string data = argv[argv_id++];
	std::string data_prefix_path = argv[argv_id++];
	std::string data_file_prefix = data_prefix_path + "/data/";
	std::string rdata_file_prefix = data_prefix_path + "/refactor/";

	struct timespec start, end;
	int err;
	double elapsed_time;

	err = clock_gettime(CLOCK_REALTIME, &start);

    if(data == "GE"){
        if (mode == 1) refactor_GE<T>(data_file_prefix, rdata_file_prefix);
        // refactor_GE_SZ3<T>(data_file_prefix, rdata_file_prefix);
        else refactor_GE_SZ3_delta<T>(data_file_prefix, rdata_file_prefix);
    }
    else if(data == "NYX" || data == "Hurricane" || data == "SCALE" || data == "Miranda" || data == "S3D"){
        if (mode == 1){
            refactor_velocities_1D<T>(data_file_prefix, rdata_file_prefix);
        }
        else if (mode == 3){
            if(data == "Hurricane"){
                refactor_Vtot_SZ3_delta<T>(100, 500, 500, data_file_prefix, rdata_file_prefix);
            }
            else if(data == "NYX"){
                refactor_Vtot_SZ3_delta<T>(512, 512, 512, data_file_prefix, rdata_file_prefix);
            }
            else if(data == "SCALE"){
                refactor_Vtot_SZ3_delta<T>(98, 1200, 1200, data_file_prefix, rdata_file_prefix);
            }
            else if(data == "Miranda"){
                refactor_Vtot_SZ3_delta<T>(256, 384, 384, data_file_prefix, rdata_file_prefix);
            }
            else if(data == "S3D"){
                refactor_Vtot_SZ3_delta<T>(500, 500, 500, data_file_prefix, rdata_file_prefix);
            }
        }
        else if (mode == 2){
            if(data == "Hurricane"){
                refactor_velocities_3D<T>(100, 500, 500, data_file_prefix, rdata_file_prefix);
            }
            else if(data == "NYX"){
                refactor_velocities_3D<T>(512, 512, 512, data_file_prefix, rdata_file_prefix);
            }
            else if(data == "SCALE"){
                refactor_velocities_3D<T>(98, 1200, 1200, data_file_prefix, rdata_file_prefix);
            }
            else if(data == "Miranda"){
                refactor_velocities_3D<T>(256, 384, 384, data_file_prefix, rdata_file_prefix);
            }
            else if(data == "S3D"){
                refactor_velocities_3D<T>(500, 500, 500, data_file_prefix, rdata_file_prefix);
            }
        }
    }
    // else{
    //     uint32_t n1 = atoi(argv[argv_id++]);
    //     uint32_t n2 = atoi(argv[argv_id++]);
    //     uint32_t n3 = atoi(argv[argv_id++]);
    //     if(data == "S3D"){
    //         refactor_S3D<T>(n1, n2, n3, data_file_prefix, rdata_file_prefix);
    //         refactor_S3D_SZ3<T>(n1, n2, n3, data_file_prefix, rdata_file_prefix);
    //         refactor_S3D_SZ3_delta<T>(n1, n2, n3, data_file_prefix, rdata_file_prefix);
    //     }
    // }
    
	err = clock_gettime(CLOCK_REALTIME, &end);
	elapsed_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000;
	printf("elapsed_time = %.6f\n", elapsed_time);

    return 0;

}