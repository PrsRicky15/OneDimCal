# Minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Project name and language
project(OneDimension LANGUAGES Fortran)

# Set Fortran standard (optional)
set(CMAKE_Fortran_STANDARD 2008)

# Set build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -Wall -Wextra -fcheck=all")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")

# Find source files
file(GLOB_RECURSE SOURCES "*.f08" "*.F08" "*.f" "*.F")

# Create executable
add_executable(${PROJECT_NAME} potential.f08 gridinfo.f08
        LinearAlgebra.f08 derivative_operator.f08
        test.f08 main.f08)

# Find BLAS and LAPACK libraries
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(FFTW3 REQUIRED)

# Display found libraries (optional, for debugging)
message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")
message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")
message(STATUS "FFTW3 libraries: ${FFTW3_LIBRARIES}")

# Optional: Set output directory
set_target_properties(${PROJECT_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Link BLAS and LAPACK libraries
target_link_libraries(${PROJECT_NAME}
        ${LAPACK_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${FFTW3_LIBRARIES}
)

install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin
)