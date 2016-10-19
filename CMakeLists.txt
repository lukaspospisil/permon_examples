#*******************************************************************************
# PERMON examples
# Lukas Pospisil
# USI Lugano, 2016
# lukas.pospisil@usi.ch
#
# *******************************************************************************/
#
# cmake -DFIND_PETSC=ON -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx ..
#

project(PERMON_EXAMPLES)
cmake_minimum_required(VERSION 2.8)

# set the root to the examples code
set(EXAMPLE_ROOT "${CMAKE_SOURCE_DIR}")

# include cmake functions
set(CMAKE_MODULE_PATH "${EXAMPLE_ROOT}/cmake/" ${CMAKE_MODULE_PATH})

# default settings
option(USE_CUDA "USE_CUDA" OFF)
option(USE_PETSC "USE_PETSC" ON)
option(USE_PERMON "USE_PERMON" ON)
option(FIND_PETSC "FIND_PETSC" ON)

# include cmake functions
include(common)

# print information about settings
include(printsetting)

# --------- include test cmake functions --------
# load variables from console
option(TEST_01 "TEST_01" OFF)

# print info
print("\nPERMON tests")
printinfo_onoff(" TEST_01            (hello_permon)      " "${TEST_01}")
message("\n----------------------------------------------------------\n")


if(${TEST_01})
	if(${USE_CUDA})
		add_executable_ext("test_01/hello_permon.cu" "test_01")
	else()
		add_executable_ext("test_01/hello_permon.cpp" "test_01")
	endif()
endif()

