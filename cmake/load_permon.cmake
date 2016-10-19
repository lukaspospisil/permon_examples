# permon stuff

if(${USE_PERMON})
	message(STATUS "${Blue}loading permon library${ColourReset}")

	# TODO: check ENV variable PERMON_DIR

	# define where can be PERMON libraries and includes found
	set(PERMON_INCLUDE "$ENV{PERMON_DIR}/include")
	set(PERMON_LIBRARIES "$ENV{PERMON_DIR}/${PETSC_ARCH}/lib/libpermon.so")

	# include header files
	include_directories(${PERMON_INCLUDE})
	if(${USE_CUDA})
		cuda_include_directories(${PERMON_INCLUDE})
	endif()

	# link library
	set(LIBRARIES_DEF "${PERMON_LIBRARIES};${LIBRARIES_DEF}")

	# append to flags definitions
	set(FLAGS_DEF "-USE_PERMON ${FLAGS_DEF}")
	set(FLAGS_DEF_D "-DUSE_PERMON ${FLAGS_DEF_D}")	
endif()



# define print info (will be called in printsetting.cmake)
macro(PRINTSETTING_PERMON)
	printinfo_onoff("USE_PERMON\t\t\t" "${USE_PERMON}")
	if(${USE_PERMON})
		printinfo(" - PERMON_DIR\t\t\t" "$ENV{PERMON_DIR}")
		printinfo(" - PERMON_INCLUDE\t\t" "${PERMON_INCLUDE}")
		printinfo(" - PERMON_LIBRARIES\t\t" "${PERMON_LIBRARIES}")
	endif()
endmacro()

