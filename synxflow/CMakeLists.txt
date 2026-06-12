find_package(CUDA)

IF(WIN32)
    set(PLATFORM_CXX_FLAGS /Zi)
    set(CUDA_SEPARABLE_COMPILATION ON)
    # arch patch (transport T1): arch overridable via -DCUDA_NVCC_FLAGS or env SYNXFLOW_CUDA_ARCH;
    # fallback sm_50. --expt-extended-lambda always kept (required by device lambdas).
    if(NOT CUDA_NVCC_FLAGS)
        if(DEFINED ENV{SYNXFLOW_CUDA_ARCH})
            set(CUDA_NVCC_FLAGS -arch=$ENV{SYNXFLOW_CUDA_ARCH};--expt-extended-lambda)
        else()
            set(CUDA_NVCC_FLAGS -arch=sm_50;--expt-extended-lambda)
        endif()
    endif()
ENDIF(WIN32)

IF(UNIX)
    set(PLATFORM_CXX_FLAGS "-pthread -std=c++11 -fPIC -O3")
    #set(CUDA_SEPARABLE_COMPILATION ON)
    # arch patch (transport T1): arch overridable via -DCUDA_NVCC_FLAGS or env SYNXFLOW_CUDA_ARCH
    # (setup.py build_ext does NOT forward -D cmake args, so env is the working injection path);
    # SYNXFLOW_CUDA_ARCH=sm_120 for Blackwell RTX 5070. fallback sm_50 + --expt-extended-lambda.
    if(NOT CUDA_NVCC_FLAGS)
        if(DEFINED ENV{SYNXFLOW_CUDA_ARCH})
            set(CUDA_NVCC_FLAGS -arch=$ENV{SYNXFLOW_CUDA_ARCH};--expt-extended-lambda)
        else()
            set(CUDA_NVCC_FLAGS  -arch=sm_50;--expt-extended-lambda)
        endif()
    endif()
ENDIF(UNIX)

list(APPEND CMAKE_CXX_FLAGS ${PLATFORM_CXX_FLAGS})

message("CXX  Flags: " ${CMAKE_CXX_FLAGS})
message("NVCC Flags: " ${CUDA_NVCC_FLAGS})

find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)

include(FetchContent)

FetchContent_Declare(
  pybind11_sources
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.12.0
)

FetchContent_GetProperties(pybind11_sources)

if(NOT pybind11_sources_POPULATED)
  FetchContent_Populate(pybind11_sources)

  add_subdirectory(
    ${pybind11_sources_SOURCE_DIR}
    ${pybind11_sources_BINARY_DIR}
    )
endif()

message("pybind11 : "${pybind11_sources_SOURCE_DIR})

add_subdirectory(lib)
add_subdirectory(apps/cudaFloodSolversPybind)
add_subdirectory(apps/cudaDebrisFlowSolverPybind)
add_subdirectory(apps/cudaLandslideSolverPybind)
