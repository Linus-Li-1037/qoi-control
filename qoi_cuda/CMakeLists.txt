function(add_my_executable target source)
    add_executable(${target} ${source})
    target_include_directories(${target} PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
    target_link_libraries(${target} ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})
endfunction()

# add_my_executable(halfing_Vtot_cuda halfing_Vtot_cuda.cu)
