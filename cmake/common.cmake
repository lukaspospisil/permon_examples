# define general functions and variables used in other scripts

# set default build type
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Debug")
endif()

# here will be loaded options for compilers
set(FLAGS_DEF "")
set(FLAGS_DEF_D "")
set(LIBRARIES_DEF "")
set(COMPILE_FIRST "")

# add debug definitions to compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

# define some colors for funny cmake messages
option(CMAKE_USE_COLOR "CMAKE_USE_COLOR" ON)
if(NOT WIN32 AND ${CMAKE_USE_COLOR})
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
endif()

# define macros for printing messages from cmake
macro(PRINT value)
 message("${value}")
endmacro()

macro(PRINTINFO name value)
 message(" ${name} : ${Yellow}${value}${ColourReset}")
endmacro()

macro(PRINTINFO_RED name value)
 message(" ${name} : ${Red}${value}${ColourReset}")
endmacro()

macro(PRINTINFO_GREEN name value)
 message(" ${name} : ${Green}${value}${ColourReset}")
endmacro()

macro(PRINTINFO_ONOFF name value)
	if(${value})
		message(" ${ColourReset}${name} : ${Green}${value}${Yellow}${ColourReset}")
	else()
		message(" ${ColourReset}${name} : ${Red}${value}${Yellow}${ColourReset}")
	endif()
endmacro()


# addition of executable file based on file extension
macro(ADD_EXECUTABLE_EXT filename outname)
	# choose compiler subject to file extension
	get_filename_component(FILE_EXT ${filename} EXT)

	if(${FILE_EXT} MATCHES ".cpp")
		# compile with g++

		# add executable file
		add_executable(${outname} ${filename})

		# link external libraries	
		target_link_libraries(${outname} ${LIBRARIES_DEF})

		# set the name of output file
		set_source_files_properties(${filename}
				COMPILE_FLAGS "${FLAGS_DEF_D}")

		set_target_properties(${outname} PROPERTIES
			OUTPUT_NAME ${outname}
#			DEBUG ${CMAKE_CXX_FLAGS_DEBUG}
		)
	endif()	
	
endmacro()

include(load_petsc) # PETSC
include(load_permon) # PERMON


