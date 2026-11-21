set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(lsquant_core STATIC
    core/algebra/sparse_matrix.cpp
    core/algebra/composite_operators.cpp
    core/algebra/linalg.cpp
    core/algebra/state_generator.cpp
    core/algebra/spectral_bounds.cpp
    core/chebyshev/rescaled_operator.cpp
    core/chebyshev/chebyshev_recursion.cpp
    core/chebyshev/chebyshev_vectors.cpp
    core/chebyshev/chebyshev_moments.cpp
    core/chebyshev/moment_engine.cpp
    io/csr_text.cpp
    io/chebmom_file.cpp
)

target_include_directories(lsquant_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_features(lsquant_core PUBLIC cxx_std_17)
target_compile_options(lsquant_core PUBLIC -Wall -Wextra)
target_link_libraries(lsquant_core PUBLIC Eigen3::Eigen OpenMP::OpenMP_CXX)

add_executable(stage1_moments apps/stage1_moments.cpp)
set_target_properties(stage1_moments PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
target_compile_options(stage1_moments PRIVATE -Wall -Wextra)
target_link_libraries(stage1_moments PRIVATE lsquant_core)

add_executable(stage_1_check "${PROJECT_SOURCE_DIR}/tests/dev/stage_1_check.cpp")
set_target_properties(stage_1_check PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
target_compile_options(stage_1_check PRIVATE -Wall -Wextra)
target_link_libraries(stage_1_check PRIVATE lsquant_core)

# Unit tests of the stage-1 core (tests/unit/test_core_stage1.cpp) and the dense
# self-check (tests/dev/stage_1_check.cpp), registered with CTest:  ctest -R core
add_executable(test_core_stage1 "${PROJECT_SOURCE_DIR}/tests/unit/test_core_stage1.cpp")
set_target_properties(test_core_stage1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
target_compile_options(test_core_stage1 PRIVATE -Wall -Wextra)
target_link_libraries(test_core_stage1 PRIVATE lsquant_core)

enable_testing()
add_test(NAME core_stage1_unit COMMAND test_core_stage1)
add_test(NAME core_stage1_dense_check COMMAND stage_1_check)
# The unit test uses matrices of dimension <= 40: a few threads are enough, and the
# default (one thread per core) only adds contention on a loaded node.
set_tests_properties(core_stage1_unit PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")
