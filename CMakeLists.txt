# LSQUANT — Linear Quantum Transport (KPM core)
#
# Reflects the build as it actually is: a single Eigen numeric backend plus
# FFTW3, GSL and OpenMP. The former INTEL_MKL / EIGEN / BLAS_LIB switches and the
# (nonexistent) icc+mkl toolchain are gone; dependencies are resolved with
# find_package and consumed as imported targets, so no hand-passed toolchain file
# is required:  cmake -B build && cmake --build build -j
cmake_minimum_required(VERSION 3.21)

project(LSQUANT VERSION 2.0.0 LANGUAGES CXX)

include(GNUInstallDirs)

# --- Language / build type ---------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# Current Boost.Math headers need at least C++14. Use C++17 throughout, matching
# the new core, while preserving the legacy Boost Bessel implementation.
# Warning: changing the numerical implementation can invalidate byte-exact goldens.
if(NOT "cxx_std_17" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    message(FATAL_ERROR
        "LSQuant requires a compiler with C++17 support. Select one with -DCMAKE_CXX_COMPILER=...")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# --- Options -----------------------------------------------------------------
option(LSQUANT_BUILD_UTILITIES   "Build the Python helper utilities" OFF)
option(LSQUANT_BUILD_NEW_CORE    "Build the new C++ core (lsquant/) next to the legacy drivers" ON)
option(LSQUANT_BUILD_INPUT_TOOLS "Fetch + build the external wannier2sparse for tests" ON)
# OFF by default: normal builds and the test suite never touch Doxygen/Graphviz.
option(BUILD_DOCS                "Build the Doxygen documentation (needs doxygen + graphviz)" OFF)

# New core: numerical build profile and components.
set(LSQUANT_BUILD_PROFILE "golden" CACHE STRING
    "Numerical build profile: 'golden' (portable IEEE, no -march=native or -ffast-math) or 'native' (-O3 -march=native -mtune=native -fno-math-errno, defines LSQUANT_NATIVE_BUILD=1)")
set_property(CACHE LSQUANT_BUILD_PROFILE PROPERTY STRINGS golden native)
if(NOT LSQUANT_BUILD_PROFILE MATCHES "^(golden|native)$")
    message(FATAL_ERROR "LSQUANT_BUILD_PROFILE must be 'golden' or 'native' (got '${LSQUANT_BUILD_PROFILE}')")
endif()
message(STATUS "LSQUANT build profile: ${LSQUANT_BUILD_PROFILE}")
option(LSQUANT_BUILD_PYTHON     "Build the pybind11 extension module lsquant._core" ON)
option(LSQUANT_BUILD_BENCHMARKS "Build the C++ kernel benchmark executable stage8_kernels" OFF)
option(LSQUANT_BUILD_DEV_CHECKS "Build the C++ developer self-check executable stage_1_check" ON)
option(LSQUANT_WITH_GINKGO      "Build the optional Ginkgo LinearOperator backend (native profile only)" OFF)
option(LSQUANT_WITH_CUDA        "Build the CUDA atom-block kernel; also enables the Ginkgo CUDA executor when Ginkgo is on (native profile only)" OFF)
option(LSQUANT_WITH_MPI         "Build the MPI distributed operator (native profile only)" OFF)

if(LSQUANT_WITH_GINKGO AND LSQUANT_BUILD_PROFILE STREQUAL "golden")
    message(FATAL_ERROR "LSQUANT_WITH_GINKGO requires LSQUANT_BUILD_PROFILE=native "
                        "(the golden profile must stay on the native SparseMatrix kernels)")
endif()

if((LSQUANT_WITH_CUDA OR LSQUANT_WITH_MPI) AND LSQUANT_BUILD_PROFILE STREQUAL "golden")
    message(FATAL_ERROR "LSQUANT_WITH_CUDA/LSQUANT_WITH_MPI require LSQUANT_BUILD_PROFILE=native "
                        "(the golden profile is OpenMP-only by construction)")
endif()

# Shared warning flags of the new core.
add_library(lsquant_warnings INTERFACE)
target_compile_options(lsquant_warnings INTERFACE
    $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>)

# --- Dependencies ------------------------------------------------------------
# Boost.Math is header-only here; no compiled Boost libraries are needed.
find_package(Boost QUIET)
if(NOT Boost_FOUND)
    message(FATAL_ERROR
        "Boost headers are required for the legacy Boost.Math Bessel functions. "
        "Install libboost-dev or boost-cpp, or set Boost_ROOT to the Boost installation prefix.")
endif()

# Eigen3: prefer a real config-package (Eigen3::Eigen); fall back to a header dir.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
    find_path(EIGEN3_INCLUDE_DIR
        NAMES signature_of_eigen3_matrix_library
        HINTS $ENV{CONDA_PREFIX}/include/eigen3 /usr/include/eigen3 /usr/local/include/eigen3
        DOC "Directory containing the Eigen/ folder")
    if(NOT EIGEN3_INCLUDE_DIR)
        message(FATAL_ERROR
            "Eigen3 not found. Install libeigen3-dev or pass -DEIGEN3_INCLUDE_DIR=/path/to/eigen3")
    endif()
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
    message(STATUS "Eigen3 headers: ${EIGEN3_INCLUDE_DIR}")
endif()

find_package(OpenMP REQUIRED)
find_package(GSL    REQUIRED)   # -> GSL::gsl (carries gslcblas)
find_package(Threads REQUIRED)

# FFTW3 (single precision not used). The OpenMP-threaded variant lives in a
# separate library (libfftw3_omp); link it when present.
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)
find_library(FFTW3_OMP_LIBRARY fftw3_omp)
if(NOT FFTW3_OMP_LIBRARY)
    message(STATUS "libfftw3_omp not found; building without threaded FFTW")
    set(FFTW3_OMP_LIBRARY "")
endif()

# pybind11 for the new core's Python module (found from the running Python when
# pybind11_DIR is not given).
if(LSQUANT_BUILD_NEW_CORE AND LSQUANT_BUILD_PYTHON)
    # Modern Python discovery avoids pybind11's removed FindPythonInterp/FindPythonLibs
    # modules on recent CMake.
    if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.18")
        set(_lsquant_legacy_python "${PYTHON_EXECUTABLE}")
        if(NOT Python_EXECUTABLE AND _lsquant_legacy_python)
            set(Python_EXECUTABLE "${_lsquant_legacy_python}")
        endif()
        find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
        set(PYBIND11_FINDPYTHON ON)
        set(LSQUANT_PYTHON_EXECUTABLE "${Python_EXECUTABLE}")
    endif()
    if(NOT pybind11_DIR)
        if(NOT LSQUANT_PYTHON_EXECUTABLE)
            find_program(LSQUANT_PYTHON_EXECUTABLE NAMES python3 python)
        endif()
        if(LSQUANT_PYTHON_EXECUTABLE)
            execute_process(
                COMMAND "${LSQUANT_PYTHON_EXECUTABLE}" -m pybind11 --cmakedir
                OUTPUT_VARIABLE _pybind11_cmakedir
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
                RESULT_VARIABLE _pybind11_rc
            )
            if(_pybind11_rc EQUAL 0 AND EXISTS "${_pybind11_cmakedir}")
                list(APPEND CMAKE_PREFIX_PATH "${_pybind11_cmakedir}")
            endif()
        endif()
    endif()
    find_package(pybind11 CONFIG QUIET)
    if(NOT pybind11_FOUND)
        message(FATAL_ERROR
            "pybind11 not found (required when LSQUANT_BUILD_PYTHON=ON). "
            "Install pybind11 ('pip install pybind11' or 'apt install pybind11-dev') "
            "or pass -Dpybind11_DIR=$(python3 -m pybind11 --cmakedir).")
    endif()
endif()

include_directories(${PROJECT_SOURCE_DIR}/include)

# --- Targets -----------------------------------------------------------------
add_subdirectory(src)

if(LSQUANT_BUILD_NEW_CORE)
    add_subdirectory(lsquant lsquant_core)
endif()

if(LSQUANT_BUILD_UTILITIES)
    add_subdirectory(utilities)
endif()

# Golden-input generation depends on the external wannier2sparse, which is NOT
# vendored or distributed. It is fetched + built only on demand for tests.
# (cmake/FetchWannier2Sparse.cmake wraps scripts/fetch_wannier2sparse.sh, pinned
#  to a commit of github.com/adamecius/wannier2sparse.)
if(LSQUANT_BUILD_INPUT_TOOLS)
    include(cmake/FetchWannier2Sparse.cmake OPTIONAL)
endif()

include(CTest)
if(BUILD_TESTING AND EXISTS "${PROJECT_SOURCE_DIR}/test/CMakeLists.txt")
    add_subdirectory(test)
endif()

# Documentation (opt-in; OFF by default so it never affects normal builds/tests).
if(BUILD_DOCS)
    add_subdirectory(docs)
endif()
