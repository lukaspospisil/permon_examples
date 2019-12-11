#*******************************************************************************
# PERMON examples
# Lukas Pospisil
# VSB-TU Ostrava
# lukas.pospisil@vsb.cz
#
# *******************************************************************************/
#
# cmake -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx ..
# cmake -DCMAKE_CXX_COMPILER=/home/lukas/Soft/petsc/arch-linux2-c-debug/bin/mpicxx -DTEST_01=ON ..
#

project(PERMON_EXAMPLES)
cmake_minimum_required(VERSION 2.8)

# set the root to the examples code
set(EXAMPLE_ROOT "${CMAKE_SOURCE_DIR}")

# include cmake functions
set(CMAKE_MODULE_PATH "${EXAMPLE_ROOT}/cmake/" ${CMAKE_MODULE_PATH})

# default settings
option(USE_PETSC "USE_PETSC" ON)
option(USE_PERMON "USE_PERMON" ON)
option(FIND_PETSC "FIND_PETSC" ON)

# include cmake functions
include(common)

# print information about settings
include(printsetting)

# we will store the output with results into this folder
make_directory("results")

# --------- include test cmake functions --------
# load variables from console
option(TEST_01 "TEST_01" OFF)
option(TEST_02 "TEST_02" OFF)
option(TEST_03 "TEST_03" OFF)

# print info
print("\nPERMON tests")
printinfo_onoff(" TEST_01            (hello_permon)      " "${TEST_01}")
printinfo_onoff(" TEST_02            (string)            " "${TEST_02}")
printinfo_onoff(" TEST_03            (polytope)          " "${TEST_03}")
message("\n----------------------------------------------------------\n")


if(${TEST_01})
	add_executable_ext("test_01/hello_permon.cpp" "test_01")
endif()

if(${TEST_02})
	add_executable_ext("test_02/string.cpp" "test_02")
endif()

if(${TEST_03})
	add_executable_ext("test_03/polytope.cpp" "test_03")
endif()
